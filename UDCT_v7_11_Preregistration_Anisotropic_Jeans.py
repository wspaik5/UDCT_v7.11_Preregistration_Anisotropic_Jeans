#!/usr/bin/env python3
"""UDCT v7.11 prospective constant-beta Jeans sensitivity protocol.
Default invocation prints the plan. --selftest uses synthetic data only.
--run requires an externally published registration URL; it evaluates new cases.
A URL is provenance supplied by the operator, not independently verified here.
"""
import argparse
import hashlib
import importlib.util
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import scipy
from scipy.integrate import cumulative_trapezoid, trapezoid

ROOT = Path(__file__).resolve().parent
DEPENDENCIES = {
    'UDCT_v7_10_Postregistration_run.py': 'a4d1a07ebcaa06348be3022d0f1f6ca463347e00194f7fde59f2bd9dbcd521b7',
    'UDCT_v7_06_Hessian_Term_Sign_Audit_reproduce.py': 'a64440e1350ae8bf10eb3188ec34829b6d586abf275de5cdc841a50ad3478284',
}
BETAS = (-1., -.5, 0., .25, .5)
ELLS = (1., 3., 10.)
INNER = (.002, .0005)
GRIDS = ('smoke', 'default', 'refined')
PROFILES = ('Tapered cusp', 'Plummer')
# Fixed comparison radii inside all integration domains; boundary-sensitive inner
# diagnostics are retained separately rather than silently discarded.
COMMON_X = np.geomspace(.002, 10., 1200)
PLAN = dict(version='v7.11', author='Won Shik Paik', preparation_date='2026-10-09',
    status='prospective protocol; no new galaxy cases executed during preparation',
    target='Hydrus I', betas=BETAS, screening_lengths_pc=ELLS,
    profiles=PROFILES, grids=GRIDS, inner_radii_rh=INNER,
    laws=('BASE', 'A-f1/2', 'A-f2', 'A-Qs', 'A-Qd'),
    controls=('QUMOND-nu0', 'QUMOND-nu1', 'QUMOND-nu2'),
    field_evaluations=216, beta_records=1080,
    primary='Tapered cusp; BASE plus four variants; zero outer pressure',
    comparison='Plummer and QUMOND controls; positive-beta Plummer excluded from physical-candidate claims',
    boundary='minimum nonnegative outer radial pressure required on the evaluated domain; diagnostic only',
    outcome_scope='tested spherical-average closures only; not MOND exclusion or observational rejection')


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_baseline():
    for name, expected in DEPENDENCIES.items():
        path = ROOT / name
        if not path.is_file() or digest(path) != expected:
            raise RuntimeError('Missing or modified frozen dependency: ' + name)
    spec = importlib.util.spec_from_file_location('udct711_baseline', ROOT / 'UDCT_v7_10_Postregistration_run.py')
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def solve(r, rho, inward, beta, pout=0.):
    """Signed pressure. Positive inward force. Never clip negative values.
    Integrating factor q=(r/Rout)^(2 beta), so q(Rout)=1.
    Outer pressure is in the same SI pressure units as rho*velocity^2.
    """
    r, rho, inward = [np.asarray(v, float) for v in (r, rho, inward)]
    if r.ndim != 1 or len(r) < 3 or rho.shape != r.shape or inward.shape != r.shape:
        raise ValueError('Invalid array dimensions')
    if not all(np.all(np.isfinite(v)) for v in (r, rho, inward)):
        raise ValueError('Nonfinite input')
    if np.any(r <= 0) or np.any(np.diff(r) <= 0) or np.any(rho < 0):
        raise ValueError('Invalid radius/density')
    if not np.isfinite(beta) or beta >= 1 or not np.isfinite(pout) or pout < 0:
        raise ValueError('Invalid beta/boundary pressure')
    q = (r/r[-1]) ** (2*beta)
    integ = -np.r_[0., cumulative_trapezoid((q*rho*inward)[::-1], r[::-1])][::-1]
    return (integ + pout)/q, integ, q


def critical_boundary(integ, mask):
    return float(max(0., -np.min(integ[mask])))


def summary(r, rho, inward, beta, rh, gn):
    pressure, integ, q = solve(r, rho, inward, beta)
    mask = (rho > 1e-280)
    sigma = np.divide(pressure, rho, out=np.full_like(r, np.nan), where=mask)
    x = r/rh
    reporting = mask & (x <= 10.)
    common = np.interp(COMMON_X, x[mask], sigma[mask])
    vref2 = float(np.interp(1., x, gn))*rh
    pref = float(np.interp(1., x, rho))*vref2
    idx = np.flatnonzero(reporting)[np.argmin(sigma[reporting])]
    required = critical_boundary(integ, reporting)
    # Finite-shell virial identity, including BOTH boundaries, for the same
    # signed inward force used by Jeans. This is a diagnostic, not an assertion
    # that the spherical closure supplies a full three-dimensional equilibrium.
    kinetic2 = 4*np.pi*trapezoid((3-2*beta)*pressure*r*r, r)
    work = 4*np.pi*trapezoid(rho*inward*r**3, r)
    boundary = 4*np.pi*(r[-1]**3*pressure[-1]-r[0]**3*pressure[0])
    norm = max(abs(work), abs(kinetic2), abs(boundary), 1e-300)
    return dict(beta=beta, common_sigma2_over_vref2=(common/vref2).tolist(),
        common_sigma2_kms2=(common/1e6).tolist(),
        min_sigma2_kms2=float(sigma[idx]/1e6), min_at_r_over_rh=float(x[idx]),
        minimum_at_inner_endpoint=bool(idx == 0),
        negative_band_r_over_rh=([float(x[reporting & (sigma<0)].min()), float(x[reporting & (sigma<0)].max())]
            if np.any(reporting & (sigma<0)) else None),
        required_pout_si=required, required_pout_over_pref=required/pref,
        boundary_scope='removes negative pressure only on valid sampled radii up to 10 rh; not physical confinement evidence',
        finite_shell_2K_si=float(kinetic2), finite_shell_work_si=float(work),
        finite_shell_boundary_si=float(boundary),
        finite_shell_virial_relative_residual=float((kinetic2-work-boundary)/norm),
        isotropic_reproduction_error=None, vref2_si=vref2)


def decide(arrays):
    """Require both lower-limit runs and default/refined grids. Smoke is audit.
    One conservative error bound over the whole fixed comparison domain.
    Exact-zero or tolerance-scale values are indeterminate, never a rescue.
    """
    a, b, c, d = [np.asarray(v, float) for v in arrays]
    if any(v.shape != COMMON_X.shape or not np.all(np.isfinite(v)) for v in (a,b,c,d)):
        return dict(status='I', reason='invalid or incomplete comparison arrays')
    error = max(float(np.max(abs(a-b))), float(np.max(abs(c-d))), float(np.max(abs(b-d))))
    tau = max(1e-8, 5*error)
    lows = [float(np.min(v)) for v in (a,b,c,d)]
    # Require an overlapping negative radius in every comparison array.
    if np.any(np.maximum.reduce((a,b,c,d)) < -tau):
        status = 'N'
    elif min(lows) > tau:
        status = 'P'
    else:
        status = 'I'
    return dict(status=status, tolerance_over_vref2=tau,
                estimated_error_over_vref2=error, minima_over_vref2=lows)


def selftest():
    r = np.geomspace(.1, 3., 20000)
    rho = np.ones_like(r)
    g = np.ones_like(r)*2
    for beta in BETAS:
        p, integ, q = solve(r, rho, g, beta, .7)
        exponent = 2*beta+1
        expected = (2*np.log(r[-1]/r)/r[-1]**(2*beta) if exponent == 0
            else 2*(r[-1]**exponent-r**exponent)/exponent/r[-1]**(2*beta))
        expected = (expected+.7)/q
        assert np.max(abs(p-expected))/np.max(abs(expected)) < 2e-6
    p, integ, q = solve(r, rho, -g, -.5)
    need = critical_boundary(integ, np.ones_like(r, bool))
    corrected, _, _ = solve(r, rho, -g, -.5, need)
    assert np.min(corrected) >= -1e-10
    assert np.min(p) < 0
    assert decide([np.full(1200, -1.)]*4)['status'] == 'N'
    assert decide([np.full(1200, 1.)]*4)['status'] == 'P'
    assert decide([np.zeros(1200)]*4)['status'] == 'I'
    bad = np.ones(1200); bad[0] = np.nan
    assert decide([bad]*4)['status'] == 'I'
    try:
        solve(r, rho, g*np.nan, 0.)
        raise AssertionError('Nonfinite input accepted')
    except ValueError:
        pass
    base = load_baseline()  # hash check and import; does NOT calculate a galaxy
    st = dict(r=r, rho=rho, rh=1.)
    old, sigma = base.V706.aperture_jeans(st, g, 1.)
    p, _, _ = solve(r, rho, g, 0.)
    assert np.allclose(p, sigma, rtol=1e-12, atol=1e-12)
    print('ALL PASSED: synthetic solver, boundary threshold, decisions, dependency hashes, beta=0 adapter. No galaxy cases executed.')


def run(registration_url):
    base = load_baseline()
    records = {}
    gates = []
    laws = (base.BASE_LAW,)+base.PART_A_LAWS
    obj = base.DWARFS['Hydrus I']
    for label in GRIDS:
        grid = base.make_grid(label)
        e, h, _ = base.V706.galactic_field(obj, grid)
        for profile in PROFILES:
            specs = [(law, ell) for law in laws for ell in ELLS]
            specs += [(base.Law('QUMOND-'+kind, 1., kind), None) for kind in ('nu0','nu1','nu2')]
            for law, ell in specs:
                for inner in INNER:
                    ratios = base.mesh(profile, grid, rmin=inner)
                    st = base.model_state(obj, profile, ratios, e, h, grid, law, ell, 0.)
                    work, inward = base.V706.work_from_state(st)
                    _, old_sigma = base.V706.aperture_jeans(st, inward, obj.re_pc)
                    for beta in BETAS:
                        key = f'{profile}|{law.name}|ell={ell}|beta={beta:g}'
                        rec = summary(st['r'], st['rho'], inward, beta, st['rh'], st['gn'])
                        rec['v710_work_ratio'] = work['virial']
                        if beta == 0:
                            p, _, _ = solve(st['r'], st['rho'], inward, 0.)
                            sigma = np.divide(p, st['rho'], out=np.zeros_like(p), where=st['rho']>1e-280)
                            err = float(np.max(abs(sigma-old_sigma))/rec['vref2_si'])
                            rec['isotropic_reproduction_error'] = err
                            gates.append(err <= 1e-10)
                        records.setdefault(key, {})[f'{label}|inner={inner:g}'] = rec
                    print(f'completed {label} {profile} {law.name} ell={ell} inner={inner}', file=sys.stderr)
    decisions = {}
    for key, entries in records.items():
        arrays = [entries[f'{label}|inner={inner:g}']['common_sigma2_over_vref2']
                  for inner in INNER for label in ('default','refined')]
        dec = decide(arrays)
        beta = entries['refined|inner=0.0005']['beta']
        if dec['status'] == 'P':
            tau = dec['tolerance_over_vref2']
            relevant = [entries[f'{lab}|inner={inn:g}'] for inn in INNER for lab in ('default','refined')]
            if any(rec['min_sigma2_kms2']*1e6/rec['vref2_si'] <= tau for rec in relevant):
                dec['status'] = 'I'
                dec['reason'] = 'native-grid minimum fails positive-candidate margin'
        dec['physical_candidate_eligibility'] = not (key.startswith('Plummer|') and beta > 0)
        dec['eligibility_note'] = 'necessary central-slope screen only; no positive distribution function established'
        if not all(gates):
            dec['provisional_status'] = dec['status']; dec['status'] = 'I'
            dec['reason'] = 'beta=0 reproduction gate failed'
        decisions[key] = dec
    return dict(protocol=PLAN, registration_url=registration_url,
        registration_url_verified=False, started_or_completed_utc=datetime.now(timezone.utc).isoformat(),
        environment=dict(python=sys.version, numpy=np.__version__, scipy=scipy.__version__, platform=platform.platform()),
        runner_sha256=digest(Path(__file__)), frozen_dependencies=DEPENDENCIES,
        isotropic_gate_passed=all(gates), common_r_over_rh=COMMON_X.tolist(),
        records=records, decisions=decisions,
        note='N: converged negative band; P: positive second-moment candidate only; I: indeterminate. No full-equilibrium, observational or universal-MOND claim.')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    modes = p.add_mutually_exclusive_group()
    modes.add_argument('--selftest', action='store_true')
    modes.add_argument('--run', action='store_true')
    p.add_argument('--registration-url')
    p.add_argument('--output', type=Path)
    args = p.parse_args()
    if args.selftest:
        selftest(); return
    if not args.run:
        print(json.dumps(PLAN, indent=2)); return
    if not args.registration_url or not args.registration_url.startswith('https://'):
        p.error('--run requires --registration-url https://... after public registration')
    if not args.output:
        p.error('--run requires --output; archive the complete signed results')
    if args.output.exists():
        p.error('Output exists; choose a new path to preserve prior results')
    selftest()
    result = run(args.registration_url)
    with args.output.open('x') as f:
        json.dump(result, f, indent=2, allow_nan=False)
    print('Results saved to '+str(args.output))

if __name__ == '__main__':
    main()
