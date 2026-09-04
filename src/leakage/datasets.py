"""Small, explicit I/O helpers shared by injection validation and the notebook."""
from pathlib import Path
import hashlib
import json


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify_dataset(folder, manifest_name):
    folder = Path(folder)
    metadata = json.loads((folder / manifest_name).read_text())
    names = [f"J{i:02d}" for i in range(25)]
    for suffix in ("par", "tim"):
        assert sorted(p.stem for p in (folder / suffix).glob(f"*.{suffix}")) == names
    assert len(metadata["files"]) == 50
    for relative, digest in metadata["files"].items():
        assert sha256(folder / relative) == digest, relative
    return metadata


def load_pulsar(parfile, timfile):
    # Supplying the model preserves CLK=TT(BIPM2016), essential for a round trip.
    from pint.models import get_model
    from pint.toa import get_TOAs
    from pint.residuals import Residuals
    from pta_replicator.simulate import SimulatedPulsar

    model = get_model(parfile)
    toas = get_TOAs(timfile, model=model, ephem="DE440", planets=True, usepickle=False)
    return SimulatedPulsar(
        ephem="DE440", model=model, toas=toas, residuals=Residuals(toas, model),
        name=model.PSR.value, loc={"RAJ": model.RAJ.value, "DECJ": model.DECJ.value},
        added_signals={}, added_signals_time={},
    )
