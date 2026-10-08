# UDCT v7.11 - Preregistration: Constant-Anisotropy Jeans and Outer-Boundary Sensitivity in Hydrus I

## Registration and repository

v7.11 reserved Zenodo DOI: https://doi.org/10.5281/zenodo.23244522

v7.11 GitHub repository: https://github.com/wspaik5/UDCT_v7.11_Preregistration_Anisotropic_Jeans

The DOI is reserved for this deposit; it becomes registered when the upload is published. DOI reservation alone is not a public preregistration timestamp.

## Status and question

Won Shik Paik | Independent Researcher | Auckland, New Zealand | 9 October 2026

Prospective follow-up protocol, prepared after the v7.10 results were known. This is not an independent blind prediction of the earlier negative cases. No v7.11 galaxy cases have been executed during preparation. Synthetic solver and decision tests only have passed. Public preregistration is not complete until this package is deposited with a public timestamp.

The primary question is whether the negative radial Jeans pressure of the specified GQUMOND implementation persists when constant velocity anisotropy replaces isotropic closure. The secondary diagnostic asks what nonnegative outer radial pressure would remove sampled negative pressure.

This protocol studies a spherical-average, stationary, non-streaming second-moment closure. It does not solve a full three-dimensional equilibrium, introduce LMC/SMC fields, fit observed velocities, or establish distribution-function positivity.

## Known evidence and prospective expectations

Known before this protocol: v7.10 reports G2/J1/K2 under its frozen assumptions. Its global work sign reversal includes the f1/2, 10 pc Hydrus I cases. Negative Jeans pressure and central-boundary sensitivity are already known; these facts motivate this follow-up.

Directional expectation: changing beta can change the local pressure, but is not assumed to remove or preserve its sign. Both outcomes and unresolved cases will be reported. No probability or claim of a decisive MOND exclusion is preregistered.

Fixed density and gravity imply that the previously defined global work does not depend on beta. A local nonnegative Jeans profile therefore cannot by itself rescue a negative global-work case. Finite inner and outer pressure terms must be retained in the virial audit.

## Frozen design

Target: Hydrus I only. Stellar inputs, Galactic position, mass-to-light ratio, Milky Way acceleration and Hessian, density normalization, force conventions and action-consistent gravity remain those of the frozen v7.10 runner and its v7.06 dependency.

Primary density: unsoftened Tapered cusp. Comparison density: Plummer. Functions: BASE, A-f1/2, A-f2, A-Qs, A-Qd. Screening lengths: 1, 3, 10 pc. QUMOND controls: nu0, nu1, nu2, with screening absent; controls are not duplicated by length.

Constant beta = -1, -0.5, 0, 0.25, 0.5. Negative beta is tangential bias; positive beta is radial bias. The general definition is beta = 1 - (sigma_theta^2 + sigma_phi^2)/(2 sigma_r^2). This is a finite family, not all possible anisotropy profiles.

Both profiles and all functions/controls run on smoke/default/refined v7.10 grids. Lower limits: 0.002 and 0.0005 baseline r_h. Outer radii: 40 r_h for cusp and 80 r_h for Plummer, unchanged. Three grids x two profiles x (15 screened fields + three controls) x two inner limits = 216 field evaluations; five betas give 1080 pressure records.

Primary classification uses common radii 0.002 to 10 r_h (1200 logarithmic points), on default/refined grids and both inner limits. Native-grid minima, negative intervals, and minima at the inner endpoint are retained separately. Smoke results are audits, not the deciding grid.

Central-slope eligibility: positive-beta Plummer cases are diagnostic only, excluded from physical-candidate claims under the regular-center spherical interpretation. Cusp beta <= 0.5 cases are only necessary-condition candidates. Actual applicability of slope constraints depends on the potential; a positive distribution function and regular center are not established by this filter.

## Jeans solver and boundary diagnostic

Use P_r = rho_star sigma_r^2 and positive inward acceleration g_in. Solve dP_r/dr + 2 beta P_r/r = -rho_star g_in. Set q(r)=(r/R_out)^(2 beta), so q(R_out)=1. Then P_r(r)=[P_out + I(r)]/q(r), where I(r)=integral from r to R_out of q(s) rho_star(s) g_in(s) ds.

Primary runs set P_out=0. The signed cumulative trapezoidal integral is evaluated inward. Negative pressure and variance remain signed; no clipping, absolute-value conversion, or square-root replacement is permitted. Tail cells with rho <= 1e-280 are excluded, matching the inherited density-underflow convention.

For valid native radii up to 10 r_h, P_out,crit = max(0, -min I(r)). Report this threshold in SI pressure units and divided by P_ref=rho_star(r_h) r_h g_N,star(r_h). This is the exact threshold on the sampled domain for the discrete solver, not an all-radius proof.

Nonzero pressure at a far, near-vacuum boundary may be physically implausible. The threshold is a diagnostic of boundary dependence, not evidence that collisionless stars receive gas pressure or that an externally confined equilibrium exists. Truncation, a different outer radius, a variable beta, and softened profiles require a separately labelled extension.

## Gates, tolerances and decision rules

Gate R: beta=0 and P_out=0 must reproduce the inherited isotropic sigma_r^2 array on every field evaluation to maximum absolute error / v_ref^2 <= 1e-10, where v_ref^2 = r_h g_N,star(r_h). Missing or modified frozen dependencies abort execution. Any Gate R failure suppresses all outcome claims.

At each case and beta, form a=default/0.002, b=refined/0.002, c=default/0.0005, d=refined/0.0005 on the common radii, divided by v_ref^2. Define E=max(max|a-b|,max|c-d|,max|b-d|) and tau=max(1e-8,5 E). This is a conservative numerical diagnostic, not a rigorous error bound.

N (persistent negative band): at least one common radius has all four values below -tau. P (positive second-moment candidate): every common value exceeds tau; inspect native minima as well. I (indeterminate): all other cases, invalid/nonfinite data, missing comparisons, or Gate R failure. Exact zero is I. Report every case and all thresholds, including exclusions.

A numerical P is never called a physical equilibrium or observational success. Report which tested, eligible beta values give N/P/I for each function and length. Do not claim that all anisotropies are excluded. If QUMOND controls also give negative or indeterminate results, flag the closure/numerical/control problem before attributing a result specifically to screening.

Finite-shell audit: 2K_shell = V_shell + 4 pi [R_out^3 P_r(R_out) - r_in^3 P_r(r_in)], with 2K_shell=4 pi integral (3-2 beta)P_r r^2 dr and V_shell=4 pi integral rho_star g_in r^3 dr. Report the signed normalized residual. It is a diagnostic alongside the original action-consistent work, not a full 3D virial proof.

## Execution, provenance and interpretation

Before public registration, run only --selftest or the default plan display. The self-test checks constant-force analytic solutions (including beta=-0.5), signed negative pressure, critical boundary pressure, invalid data, exact-zero decisions, hashes, and the beta=0 adapter on synthetic arrays. It does not evaluate a galaxy.

After public registration, use --run --registration-url with the actual https registration URL and --output with a new JSON filename. The URL is operator-supplied and not verified by the script. Archive the public timestamp, package hashes, full command, environment, complete output, failures and deviations. The script refuses to overwrite a prior output.

No optional stopping or outcome-based deletion. Run the complete matrix. If resource limits prevent completion, archive the partial attempt and make no complete-set claim. Fixes that can affect registered quantities require a timestamped amendment before the affected cases are rerun.

This preparation includes no post-registration results. A later report must distinguish known v7.10 evidence, registered v7.11 outcomes, and exploratory additions. LMC/SMC environment, non-sphericity, distribution functions and observed stellar velocities are future tests.

## Sources and frozen files

v7.10 post-registration results: https://doi.org/10.5281/zenodo.23233837

v7.10 preregistration: https://doi.org/10.5281/zenodo.23200693

Results repository: https://github.com/wspaik5/UDCT_v7.10_Postregistration_Results

Preregistration repository: https://github.com/wspaik5/UDCT_v7.10_Preregistration_GQUMOND_Sensitivity

Mamon & Boue (2010), Kinematic deprojection and mass inversion of spherical systems of known velocity anisotropy. MNRAS 401, 2433-2450. https://arxiv.org/abs/0906.4971

An & Evans (2006), A Cusp Slope - Central Anisotropy Theorem. ApJ 642, 752-758. https://arxiv.org/abs/astro-ph/0511686

Milgrom (2023), Generalizations of quasilinear MOND (QUMOND). Phys. Rev. D 108, 084005. https://arxiv.org/abs/2305.01589

Frozen v7.10 SHA-256: a4d1a07ebcaa06348be3022d0f1f6ca463347e00194f7fde59f2bd9dbcd521b7

Frozen v7.06 SHA-256: a64440e1350ae8bf10eb3188ec34829b6d586abf275de5cdc841a50ad3478284

v7.11 runner SHA-256: 3f00abae251d7be64a4ca5abb48b4f550644b21c946206571fd448464b13ade2

The reserved v7.11 DOI is 10.5281/zenodo.23244522. Record the public registration URL and publication timestamp in execution provenance after publication.

## Files and commands

Keep these files in the same directory:

- UDCT_v7_11_Preregistration_Anisotropic_Jeans.py: prospective runner.
- UDCT_v7_10_Postregistration_run.py: frozen gravity implementation.
- UDCT_v7_06_Hessian_Term_Sign_Audit_reproduce.py: frozen dependency.
- requirements.txt: inherited pinned environment.
- UDCT_v7_11_Preregistration_Anisotropic_Jeans_2026-10-09.pdf: protocol.

Install and validate (synthetic tests only):

```bash
python -m pip install -r requirements.txt
python UDCT_v7_11_Preregistration_Anisotropic_Jeans.py
python UDCT_v7_11_Preregistration_Anisotropic_Jeans.py --selftest
```

Only AFTER this Zenodo deposit has been publicly published:

```bash
python UDCT_v7_11_Preregistration_Anisotropic_Jeans.py --run --registration-url https://zenodo.org/records/23244522 --output v7_11_results.json
```

The reserved record URL is not evidence of publication; verify that the record is public before execution. The package has no newly computed galaxy results. Do not run --run before public registration if prospective status is intended.
