# Validation records

`cw_benchmark/` verifies the current canonical parent and the single RA=6h, Dec=-45°, fCW=1.5/T injection. `scripts/validate_cw_notebooks.py` retains the old raw/weighted-centered PINT checks, now imports their extracted definitions and reads the new data layout; it no longer executes obsolete batch or map notebooks.

`cartography/` verifies the upstream reconstruction, covariance use, source-coordinate convention and MATLAB pixel conversion. Its independent reference algebra is validation only.

`restructure/migration.json` records preserved inputs and removed paths.

`cw_residual_20260903/` contains immutable **historical injection-validation evidence for the old 5.5/T experiment**. Its path/hash references describe that historical run, not current inputs. The obsolete failed notebook has been replaced by its compact error record; no archive notebooks are retained. Historical JSON logs and the two residual diagnostic figures remain because they document the PINT round-trip repair.
