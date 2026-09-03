"""Cache the unchanged DE440 kernel using the NANOGrav ephemeris mirror."""

import hashlib
import json
from pathlib import Path

from astropy.utils.data import download_file, import_file_to_cache
from jplephem.spk import SPK

CANONICAL = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp"
MIRROR = "https://data.nanograv.org/static/data/ephem/de440.bsp"
PINT_KEY = "ftp://ssd.jpl.nasa.gov/pub/eph/planets/bsp/de440.bsp"

if __name__ == "__main__":
    path = Path(download_file(CANONICAL, sources=[MIRROR, CANONICAL],
                              cache=True, timeout=120, show_progress=True))
    with SPK.open(str(path)) as kernel:
        # This is the full DE440 kernel, not DE440s or a different ephemeris.
        assert path.stat().st_size == 119799808, "Unexpected full DE440 file size"
        assert kernel.segments and all(s.start_jd < 2453000 and s.end_jd > 2456700
                                       for s in kernel.segments)
        description = str(kernel)
    import_file_to_cache(PINT_KEY, str(path))
    record = {"ephemeris": "DE440", "canonical_url": CANONICAL,
              "download_mirror": MIRROR, "cache_path": str(path),
              "bytes": path.stat().st_size,
              "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
              "segments": description}
    root = Path(__file__).resolve().parents[1]
    output = root / "validation" / "cw_residual_20260903" / "ephemeris.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record, indent=2))
