# Audited upstream implementation and compatibility

The submodule is pinned to `mpta-gw/cartography@67beabce445b3391bce488c8162581542d6608b6` (2025-09-18). Its tracked files are unchanged. Reviewed: README, `OS_Anisotropy_Analysis_defiant.py`, the example data/noise/model configuration, `PTA_AnalysisUtils`, `maps_mod`, both conversion drivers, their MATLAB call tree, PSF code and sbatch examples.

## DEFIANT and PTA model

DEFIANT **0.4.4 at `25d38082bbf2b9b33829f59438ac66a5dfc112fa`** (2025-03-04) has the driver's seven-output `compute_PFOS` interface, `narrowband=True`, `return_pair_vals=True`, selected frequencies and full pair covariance. This replaces the previously installed 1.0.1 development commit, whose PFOS returns a dictionary. The tuple interface, 300-pair ordering, shapes and numerical results were actually exercised with the existing enterprise 3.5.0 / PINT 1.1.5 environment.

`select_freq` is a **zero-based array index**, including in this historical commit. Physical Fourier k=1 uses `select_freq=0`. The upstream CLI forwards its `fbin` argument without this conversion; blindly copying that would mislabel frequencies.

The notebook constructs the minimal enterprise equivalent of the official WN example: `TimingModel(use_svd=True)`, `MeasurementNoise(efac=1)`, and ten `gw` Fourier components with gamma=0. The official full `PTA_model.py` imports unpinned `enterprise_extensions` and passes custom `orf_matrix`/`psrs_pos` keywords. Its unrelated model branches are not required or copied here. The canonical par files have only the implicit Offset design column. No extra timing parameters are fitted.

DEFIANT requires a positive `gw` spectrum. The explicit numerical floor `gw_log10_A=-20` is about 1.1e-12 of white-noise Fourier variance. It does not inject a stochastic background. Full pair covariance is evaluated, with off-diagonal norm fraction about 8.7e-13 in this white-noise limit. It is **not a physical covariance model for a deterministic CW**; these matrices are not used for significance or confidence claims. `sigk` and the full `covk` are separately retained as returned.

## Runtime patch (no alternative reconstruction)

`patches/cartography-runtime.patch` is applied only to hash-verified copies of the five upstream `maps_mod` Python files under ignored `data/cache/cartography/<patch-sha>/`. `src/leakage/cartography_runtime.py` refuses a different submodule commit, modified upstream files, or a changed source hash. It never patches the installed enterprise or SciPy packages.

The explicit changes are:

1. Import the **bundled** `anis_coefficients.py` rather than the unmodified enterprise module, which lacks `sph_complex`.
2. Remove the unsupported `to_source=True` keyword. The bundled function has no such argument. Its response uses the propagation-coordinate convention. The final pixel coordinates therefore receive the explicit antipodal relabelling described below. No guessed replacement of that missing function is supplied.
3. Replace removed `np.math.factorial` with `math.factorial` and the removed NumPy string dtype `'complex_'` with `np.complex128`.
4. Adapt removed SciPy `sph_harm(m,l,phi,theta)` to `sph_harm_y(l,m,theta,phi)`, preserving the argument and phase conventions.

Fisher construction, dirty-map construction, SVD and `max_lkl_clm` remain the upstream implementation. The ordinary regularisation branch calls removed `scipy.linalg.pinvh(cond=...)` and supplies `sv[cutoff]`; historical SciPy `cond` has version-dependent/deprecated relative-threshold semantics. We do **not** silently replace it with an absolute threshold.

Instead the notebook supplies the current M's SVD factors to the **existing upstream `U_0/Vh_0` branch**. That branch inverts the selected upper-left subspace, and `cutoff=N_KEEP` is its subspace size. This is distinct from the driver's `cutoff-1` convention. Comparison with a Hermitian eigendecomposition gives relative inverse errors below 1e-13. The original branch allocates a real intermediate inverse; imaginary roundoff in this projected diagonal block is negligible and is covered by those checks.

The CLI also contains independent defects: `for pp, psr in PSRs` instead of enumeration, use of undefined `args.incMonopole`, `args.data`, and `save_complex_matrix`; the sbatch options do not all match the parser. The notebook calls the actual library methods rather than claiming this CLI runs unchanged.

## Official pixel conversion and coordinates

`scripts/cartography_clean_map.m` contains only input/output and the explicit coefficient permutation from `complex_analysis/maps_conversion_nocleandistribution.m`. It calls the upstream `getLMvec`, then **`makemap(pOpt,1)` → `plm2plmreal` → `spherelib/plm2xyz`**. These functions execute unchanged in pinned **Octave 9.4.0**. Native MATLAB itself was not available/tested.

The result is a 360×181 regular RA/Dec grid, not a new HEALPix synthesis. Python merely relabels propagation coordinates as source coordinates: RA → (RA+12 h) mod 24 h, Dec → −Dec. A separate antenna-product check against the CW injection convention validates this transformation at relative error ~1.2e-15. Sparse complex-harmonic evaluations independently check the official MATLAB permutation, phase and normalisation at errors below 1e-13; those evaluations are tests, never production maps.

A single angular-integration check at fixed lmax=6 compares the official NSIDE=16 response with NSIDE=32: the relative Frobenius difference in `Gamma_lm` is 0.004993 (about 0.50%). This is a quadrature diagnostic, not a bound on every map pixel. It is separate from the much smaller algebra/conversion errors above. The production NSIDE=16 follows the upstream example. Octave's deprecated-function warnings and shutdown diagnostic are retained in each conversion log; conversion exits successfully and the pixel products pass the independent checks.

Upstream `makemap` uses its default `norm=0`; output is pointwise angular power, not power multiplied by pixel solid angle. The official all-in-one conversion drivers also generate S/N, radiometer, sensitivity and Monte Carlo products. Those branches are deliberately not executed. There is another upstream inconsistency: `map_conversion_nomonopole.m` assigns `pOpt(1)=0` **after LIGO reordering**. For lmax>0, `getLMvec` puts (lmax,−lmax) first, not (0,0), so this assignment does not actually remove the monopole as its README claims. This benchmark never executes that assignment and preserves its clean reconstruction's monopole.

## Power definition and scientific limits

`anis_pta(os=Sk)` divides rho and sig by Sk and covariance by Sk². Saved **`clms`, M, X, Mprime_inv and covariance are in that native normalised convention**. `official_pixel_map.mat` and `clean_normalized` retain exactly `makemap(clms)`.

The final four-panel figure displays **`Sk * clean_normalized`**, restoring the PFOS correlated residual-power scale in s² (PSD/T). All panels have the same linear colour scale. This restoration is essential for displaying the leakage amplitude across frequency bins; it is not division by uncertainty, and it is not a strain² calibration or a claim of physical point-source luminosity. Both normalisations are saved. Signed values are kept, without clipping or removing the monopole.

This is the official unpolarised correlation-power reconstruction applied to a coherent fixed-polarisation, Earth-term CW. Its broadened, displaced hotspot and signed sidelobes are not a replacement for a coherent CW likelihood. Its per-bin values need not reproduce a pure Fejér curve: real signals include both frequency signs, and PFOS includes the timing projection and angular reconstruction.

The preserved canonical data are uniform in **nominal observing MJD**, not exactly in enterprise barycentric TOAs. The actual baseline span is 314496938.28331375 s versus 314496000 nominal seconds; maximum inherited cadence departure is ~120.84 s over 14 days. T is fixed to the former for both injection and Fourier basis. No new irregular-sampling experiment or timing-model modification was performed. PINT round-trip error is below 0.78 ns, much smaller than the canonical 100 ns TOA uncertainty.

Reproduction records: `run_metadata.json`, per-bin products, `validation/cartography/numerical_checks.json`, `validation/cw_benchmark/`, and the two Linux environment lock files. The prior stand-alone `maps` package is not an analysis dependency and was removed from the local environment.
