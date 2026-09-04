# Current research workflow

**Scientific question:** off-bin CW frequency leakage and its manifestation in per-frequency PTA sky maps.

**Current benchmark:** 25 pulsars; RA = 6 h, Dec = −45°; fCW = 2.5/T; Earth term only (`psrTerm=False`, `evolve=False`). Chirp mass 10⁹ M☉, distance 100 Mpc, phase/polarisation 0, inclination π/3.

**Data:** `data/baseline/ideal_25psr/` is the byte-preserved canonical `baseline_v1`. `data/injections/cw/ra6h_dec-45_f2p5T/` records its parent SHA256 and all injection parameters. T = 314496938.28331375 s is the baseline's global barycentric span; fCW = 7.949202983171435 nHz. Nominal 14-day cadence acquires small inherited barycentric deviations; these are recorded, not explored here.

**Mapping:** `mpta-gw/cartography` submodule @ `67beabce445b3391bce488c8162581542d6608b6`; DEFIANT 0.4.4 @ `25d38082bbf2b9b33829f59438ac66a5dfc112fa`; narrowband PFOS with full pair-covariance computation → modified MAPS Fisher/dirty map → official regularised clean reconstruction → official MATLAB `makemap` run by Octave 9.4.0.

**Manual decision:** measured Fisher spectra support the recorded exploratory choice lmax = 6, 29 retained modes, retained condition ≈13.23. The spectrum does not select a unique cutoff. Every bin repeats the same checks.

**Current notebook:** `notebooks/01_MeerKAT_clean_power_hotspot.ipynb`. Products are separated into `fbin_01/` … `fbin_04/`. Saved coefficients retain upstream Sk normalisation; the shared-scale figure restores Sk and shows signed clean correlation power in s². At the exact RA = 6 h, Dec = −45° grid point, the notebook evaluates one radiometer scalar per bin. It does not construct a radiometer map. No division by pixel uncertainty, no monopole removal, no clipping.

**Scope:** one injection, four bins. Compatibility differences and physical limits: [external/COMPATIBILITY.md](external/COMPATIBILITY.md). `src/leakage/cw_response.py` is independent frequency-response code and is not used for mapping.
