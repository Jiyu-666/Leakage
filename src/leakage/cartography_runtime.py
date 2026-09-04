"""Load audited upstream code with a recorded compatibility patch, never a map reimplementation."""
from pathlib import Path
import hashlib
import importlib
import json
import subprocess
import sys

COMMIT = "67beabce445b3391bce488c8162581542d6608b6"


def load_cartography(root):
    root = Path(root)
    upstream = root / "external/cartography"
    revision = subprocess.check_output(["git", "-C", str(upstream), "rev-parse", "HEAD"], text=True).strip()
    if revision != COMMIT:
        raise RuntimeError(f"Wrong cartography revision: {revision}")
    # Check all upstream tracked files, including the MATLAB implementation.
    subprocess.run(["git", "-C", str(upstream), "diff", "--exit-code", "HEAD"], check=True, capture_output=True)
    patch = root / "external/patches/cartography-runtime.patch"
    key = hashlib.sha256(patch.read_bytes()).hexdigest()
    cache = root / "data/cache/cartography" / key
    package = cache / "leakage_cartography_maps"
    manifest = json.loads((root / "external/patches/upstream-sha256.json").read_text())
    package.mkdir(parents=True, exist_ok=True)
    for name, expected in manifest.items():
        original = upstream / "PTA_AnalysisUtils/maps_mod" / name
        if hashlib.sha256(original.read_bytes()).hexdigest() != expected:
            raise RuntimeError(f"Upstream source changed: {original}")
        (package / name).write_bytes(original.read_bytes())
    (package / "__init__.py").touch()
    subprocess.run(["patch", "--batch", "--fuzz=0", "-p1", "-i", str(patch)],
                   cwd=package, check=True, capture_output=True)
    if "leakage_cartography_maps.anis_pta" in sys.modules:
        loaded = Path(sys.modules["leakage_cartography_maps.anis_pta"].__file__).parent
        if loaded != package:
            raise RuntimeError("Compatibility patch changed: restart the notebook kernel.")
    sys.path.insert(0, str(cache))
    return importlib.import_module("leakage_cartography_maps.anis_pta"), {
        "commit": revision, "compatibility_patch_sha256": key,
        "basis_coordinates": "propagation; apply explicit antipodal coordinate conversion after makemap",
    }
