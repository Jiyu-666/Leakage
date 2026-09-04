"""Validate notebook syntax, canonical lineage and the one CW injection (not mapping)."""
import argparse
from copy import deepcopy
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from leakage.datasets import load_pulsar, sha256, verify_dataset

RUN = ROOT / 'validation/cw_benchmark_f2p5T'
BASE = ROOT / 'data/baseline/ideal_25psr'
INJ = ROOT / 'data/injections/cw/ra6h_dec-45_f2p5T'


def save(name, value):
    RUN.mkdir(parents=True, exist_ok=True)
    (RUN / name).write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n')


def hashes():
    paths = [p for folder in (BASE, INJ) for p in folder.rglob('*') if p.is_file()]
    return {str(p.relative_to(ROOT)): sha256(p) for p in sorted(paths)}


def snapshot():
    if (RUN / 'input_sha256.json').exists():
        raise RuntimeError('Existing input snapshot is immutable; run verify instead.')
    save('input_sha256.json', hashes())


def verify():
    meta = verify_dataset(BASE, 'metadata.json')
    injection = verify_dataset(INJ, 'injection.json')
    assert injection['parent_metadata_sha256'] == sha256(BASE / 'metadata.json')
    assert injection['parameters']['fCW_T'] == 2.5
    assert injection['parameters']['ra_deg'] == 90
    assert injection['parameters']['dec_deg'] == -45
    assert not injection['parameters']['psrTerm'] and not injection['parameters']['evolve']
    before = RUN / 'input_sha256.json'
    if before.exists():
        assert json.loads(before.read_text()) == hashes(), 'Canonical/injection inputs changed'
    save('input_integrity.json', {'files_verified':len(hashes()), 'parent_hash_matches':True,
                                'canonical_source':meta['source_path']})
    print('PASS: baseline, injection and lineage')


def test():
    import ast
    import nbformat
    import numpy as np
    from leakage.injection_validation import cw_weighted_center, cw_error_stats
    from leakage.cw_response import exact_window_projection, fejer_power, timing_projectors

    notebooks = sorted((ROOT / 'notebooks').glob('*.ipynb'))
    assert [n.name for n in notebooks] == [
        '00_CW注入_手写原稿.ipynb',
        '01_MeerKAT_clean_power_hotspot.ipynb',
        '03_CW本征频率泄漏.ipynb',
    ]
    for p in notebooks:
        nb = nbformat.read(p,4)
        nbformat.validate(nb)
        for c in nb.cells:
            if c.cell_type == 'code': ast.parse(c.source)
    t = np.arange(261, dtype=float)
    signal = 3*np.sin(2*np.pi*1.5*t/260) + 2*np.cos(2*np.pi*1.5*t/260)
    centered = cw_weighted_center(signal,np.ones(len(t)))
    assert abs(centered.mean()) < 1e-14
    assert cw_error_stats(np.array([2.,2.,2.]))['rms_ns'] == 2
    direct, predicted, mono_error, window_error = exact_window_projection(signal,t,np.arange(1,5)/260,1.5/260)
    assert mono_error < 1e-14 and window_error < 1e-13
    assert fejer_power(261,1.,0.) == 261**2
    # Rank-deficient design: project only the actual span, not surplus QR columns.
    from types import SimpleNamespace
    design = np.column_stack((np.ones(10),np.ones(10)))
    model = SimpleNamespace(designmatrix=lambda _: (design,['a','duplicate'],None))
    projectors,diagnostics = timing_projectors([SimpleNamespace(model=model,toas=None)],1.0)
    np.testing.assert_allclose(projectors[0],np.eye(10)-np.ones((10,10))/10,atol=1e-14)
    save('synthetic_tests.json',{'exact_window_error':float(window_error),
        'monochromatic_reconstruction_error':float(mono_error), 'projector_rank':int(diagnostics.ranks[0]),
        'notebook_format_and_syntax':True, 'rms_includes_constant_bias':True})
    verify()
    print('PASS: notebook structure, retained response functions, RMS and rank-deficient projection')


def integration():
    import numpy as np
    from pint.logging import setup
    from pta_replicator.deterministic import add_cgw
    from leakage.injection_validation import diagnose_cw_injection
    setup(level='ERROR')
    meta=verify_dataset(INJ,'injection.json');cfg=meta['parameters']
    results=[]
    for name in [f'J{i:02d}' for i in range(25)]:
        baseline=load_pulsar(BASE/'par'/f'{name}.par',BASE/'tim'/f'{name}.tim')
        loaded=load_pulsar(INJ/'par'/f'{name}.par',INJ/'tim'/f'{name}.tim')
        fresh=deepcopy(baseline)
        add_cgw(fresh,gwtheta=np.pi/2-np.deg2rad(cfg['dec_deg']),gwphi=np.deg2rad(cfg['ra_deg']),
            mc=cfg['chirp_mass_solar'],dist=cfg['distance_mpc'],fgw=cfg['fCW_hz'],
            phase0=cfg['phase_rad'],psi=cfg['polarisation_rad'],inc=cfg['inclination_rad'],
            psrTerm=cfg['psrTerm'],evolve=cfg['evolve'],phase_approx=cfg['phase_approx'],
            tref=cfg['tref_seconds'],signal_name='cw')
        loaded.added_signals_time=fresh.added_signals_time
        d=diagnose_cw_injection(baseline,loaded)
        for convention in ('raw','projected'):
            assert d[convention]['max_abs_ns']<3 and d[convention]['rms_ns']<1
        results.append({k:v for k,v in d.items() if k in ('pulsar','raw','projected')})
    save('pint_roundtrip.json',{'parameters':cfg,'results':results})
    verify()
    print('PASS: all 25 CW raw and centered round trips')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=('snapshot','verify','test','integration'))
    args=parser.parse_args()
    globals()[args.action]()
