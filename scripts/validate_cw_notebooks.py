"""Snapshot inputs, test diagnostics, and execute both CW notebooks on N522."""

import argparse
import ast
import base64
import hashlib
import json
import os
from pathlib import Path
import platform
import time

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "validation" / "cw_residual_20260903"
NOTEBOOKS = ["02_CW批次注入.ipynb", "04_CW泄漏_PFOS与天空图.ipynb"]


def hashes():
    paths = [p for folder in ("demo/data", "demo/cw_batch_data")
             for p in (ROOT / folder).rglob("*") if p.is_file()]
    paths += [ROOT / "demo" / name for name in
              ("01_CW注入_手写原稿.ipynb", "03_CW本征频率泄漏.ipynb", "cw_intrinsic.py")]
    return {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(paths)}


def save(name, value):
    RUN.mkdir(parents=True, exist_ok=True)
    (RUN / name).write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


def snapshot():
    path = RUN / "input_sha256.json"
    if path.exists():
        raise RuntimeError(f"Refusing to replace original snapshot: {path}")
    save(path.name, hashes())
    print(f"Saved {len(hashes())} input hashes to {path}", flush=True)


def verify():
    before = json.loads((RUN / "input_sha256.json").read_text())
    after = hashes()
    changed = [p for p, digest in before.items() if after.get(p) != digest]
    audit = RUN / "concurrent_input_changes.json"
    audit_data = json.loads(audit.read_text()) if audit.exists() else {}
    observed = audit_data.get("observed_changes", {})
    active = audit_data.get("user_confirmed_active_paths", [])
    unexpected = [p for p in changed
                  if not any(p == a or (a.endswith("/") and p.startswith(a)) for a in active)
                  and (p not in observed or after.get(p) != observed[p]["observed_after"])]
    assert not unexpected, f"Protected inputs changed: {unexpected}"
    result = {"unchanged_original_files": len(before) - len(changed),
              "separately_recorded_concurrent_changes": changed,
              "new_files": sorted(set(after) - set(before))}
    save("input_integrity.json", result)
    print(result, flush=True)


def sources():
    notebooks = [json.loads((ROOT / "demo" / name).read_text()) for name in NOTEBOOKS]
    code = []
    for name, nb in zip(NOTEBOOKS, notebooks):
        for cell in nb["cells"]:
            if cell["cell_type"] == "code":
                ast.parse("".join(cell["source"]), filename=f"{name}:{cell['id']}")
        code.append(next("".join(c["source"]) for c in nb["cells"]
                         if c["id"] == "residual-plot"))
    assert code[0] == code[1], "Diagnostic code differs between notebooks"
    return code[0]


def test():
    import numpy as np
    import nbformat

    code = sources()
    for name in NOTEBOOKS:
        nbformat.validate(nbformat.read(ROOT / "demo" / name, as_version=4))
    definitions = ast.Module(body=[n for n in ast.parse(code).body
                                  if isinstance(n, ast.FunctionDef)], type_ignores=[])
    ns = {"np": np}
    exec(compile(definitions, "notebook-diagnostic-functions", "exec"), ns)
    center, stats = ns["cw_weighted_center"], ns["cw_error_stats"]
    t = np.linspace(0, 1, 261, dtype=np.longdouble)
    baseline = 0.3 * np.cos(2 * np.pi * 1.3 * t) + 0.2
    cases = []
    for cycles in (0, 5, 5.5):
        signal = 29 * np.sin(2 * np.pi * cycles * t)
        for unequal in (False, True):
            weights = (1 + t)**-2 if unequal else np.ones_like(t)
            delta = center(baseline + signal, weights) - center(baseline, weights)
            error = delta - center(signal, weights)
            assert np.max(np.abs(error)) < 1e-10
            assert abs(np.average(center(signal, weights), weights=weights)) < 1e-10
            assert np.max(np.abs((baseline + signal) - baseline - signal)) < 1e-10
            cases.append({"cycles": cycles, "unequal_weights": unequal,
                          "max_error_ns": float(np.max(np.abs(error)))})
    constant = stats(np.array([2.0, 2.0, 2.0]))
    assert constant["std_ns"] == 0 and constant["rms_ns"] == 2
    assert constant["mean_ns"] == 2 and constant["max_abs_ns"] == 2
    save("synthetic_tests.json", {"cases": cases, "constant_error": constant,
                                   "notebook_sources_identical": True})
    print("PASS: notebook format, syntax, matching diagnostics, 6 synthetic cases, RMS test", flush=True)


def integration():
    """Test real PINT objects without saving any simulated TOAs."""
    import numpy as np
    import astropy.units as u
    from copy import deepcopy

    ns = {}
    nb = json.loads((ROOT / "demo" / NOTEBOOKS[0]).read_text())
    os.chdir(ROOT / "demo")
    for cell in nb["cells"]:
        if cell["id"] in ("imports", "parameters", "paths", "ideal-loader", "inject-function"):
            exec("".join(cell["source"]), ns)
    definitions = ast.Module(body=[n for n in ast.parse(sources()).body
                                  if isinstance(n, ast.FunctionDef)], type_ignores=[])
    exec(compile(definitions, "notebook-diagnostic-functions", "exec"), ns)
    ideal = ns["load_ideal_psr"](ns["ideal_dir"] / "J00.par", ns["ideal_dir"] / "J00.tim")
    cw_key = "J00_cw"
    cases = []
    for unequal in (False, True):
        baseline = deepcopy(ideal)
        if unequal:
            baseline.toas.table["error"] *= np.linspace(0.5, 2.0, len(baseline.toas))
        zero = deepcopy(baseline)
        zero.added_signals_time[cw_key] = np.zeros(len(baseline.toas)) * u.s
        zero_result = ns["diagnose_cw_injection"](baseline, zero)
        assert zero_result["raw"]["max_abs_ns"] == zero_result["projected"]["max_abs_ns"] == 0
        cases.append({"cycles": 0, "unequal_weights": unequal,
                      "raw": zero_result["raw"], "projected": zero_result["projected"]})
        for cycles in (5.0, 5.5):
            injected = ns["inject_cw"]([baseline], cycles * ns["f"])[0]
            result = ns["diagnose_cw_injection"](baseline, injected)
            for convention in ("raw", "projected"):
                assert result[convention]["rms_ns"] < 1
                assert result[convention]["max_abs_ns"] < 3
            cases.append({"cycles": cycles, "unequal_weights": unequal,
                          "raw": result["raw"], "projected": result["projected"]})
            if not unequal and cycles == 5.5:
                # Astropy Time subtraction retains two-part JD precision.
                toa_shift = np.array([(a - b).to_value(u.ns) for a, b in
                                      zip(injected.toas.table["mjd"], baseline.toas.table["mjd"])])
                tdb = np.asarray(baseline.toas.table["tdbld"], dtype=np.longdouble)
                tdb_shift = (np.asarray(injected.toas.table["tdbld"], dtype=np.longdouble) - tdb) * 86400 * 1e9
                residuals = ns["Residuals"]
                pint_shift = (residuals(injected.toas, injected.model, subtract_mean=False).time_resids
                              - residuals(baseline.toas, baseline.model, subtract_mean=False).time_resids).to_value(u.ns)
                precision = {"toa_shift_minus_cw": ns["cw_error_stats"](toa_shift - result["cw_raw"]),
                             "tdbld_shift_minus_cw": ns["cw_error_stats"](tdb_shift - result["cw_raw"]),
                             "pint_shift_minus_tdbld_shift": ns["cw_error_stats"](pint_shift - tdb_shift),
                             "tdbld_spacing_ns": float(np.spacing(tdb[0]) * 86400 * 1e9),
                             "longdouble_mantissa_bits": np.finfo(np.longdouble).nmant}
    save("pint_integration_tests.json", {"cases": cases, "precision": precision})
    verify()
    print(json.dumps({"cases": cases, "precision": precision}, indent=2), flush=True)


def execute(name):
    import nbformat
    from nbclient import NotebookClient

    path = ROOT / "demo" / name
    source_digest = hashlib.sha256(path.read_bytes()).hexdigest()
    nb = nbformat.read(path, as_version=4)
    for cell in nb.cells:
        if cell.cell_type == "code":
            cell.outputs = []
            cell.execution_count = None
            cell.metadata.pop("execution", None)
    nb.metadata.kernelspec = {"display_name": "Leakage N522 (Python 3.11)",
                              "language": "python", "name": "leakage-n522-py311"}
    log = []

    def progress(cell, cell_index, **kwargs):
        entry = {"cell_index": cell_index, "cell_id": cell.get("id"),
                 "time": time.time()}
        log.append(entry)
        print(f"{name} cell {cell_index + 1}/{len(nb.cells)}: {cell.get('id')}", flush=True)

    started = time.time()
    client = NotebookClient(nb, timeout=7200, kernel_name="leakage-n522-py311",
                            resources={"metadata": {"path": str(ROOT / "demo")}},
                            allow_errors=False, on_cell_start=progress)
    try:
        client.execute()
    except BaseException:
        nbformat.write(nb, RUN / name.replace(".ipynb", ".failed.ipynb"))
        save(name + ".progress.json", log)
        verify()
        raise
    nbformat.validate(nb)
    assert all(c.execution_count is not None for c in nb.cells if c.cell_type == "code")
    verify()
    assert hashlib.sha256(path.read_bytes()).hexdigest() == source_digest, "Notebook changed during execution; refusing to overwrite"
    nbformat.write(nb, path)
    stats_cell = next(c for c in nb.cells if c.get("id") == "residual-plot")
    for output in stats_cell.outputs:
        if "image/png" in output.get("data", {}):
            (RUN / name.replace(".ipynb", ".residual.png")).write_bytes(
                base64.b64decode(output.data["image/png"]))
    output = "".join(o.get("text", "") for o in stats_cell.outputs if o.output_type == "stream")
    scan_messages = "".join(o.get("text", "") for c in nb.cells
                            for o in c.get("outputs", []) if o.output_type == "stream"
                            and ("扫描" in o.get("text", "") or "scan data" in o.get("text", "")))
    summary = {"notebook": name, "host": platform.node(), "python": platform.python_version(),
               "elapsed_seconds": time.time() - started,
               "executed_code_cells": sum(c.cell_type == "code" for c in nb.cells),
               "diagnostic_output": output, "scan_messages": scan_messages, "cell_log": log}
    save(name + ".run.json", summary)
    print(json.dumps({k: v for k, v in summary.items() if k != "cell_log"},
                     ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("snapshot", "verify", "test", "integration", "execute"))
    parser.add_argument("--notebook", choices=NOTEBOOKS)
    args = parser.parse_args()
    RUN.mkdir(parents=True, exist_ok=True)
    if args.action == "execute":
        if not args.notebook:
            parser.error("--notebook is required for execute")
        os.environ.setdefault("MPLCONFIGDIR", str(RUN / "mpl-cache"))
        execute(args.notebook)
    else:
        globals()[args.action]()
