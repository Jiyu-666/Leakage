"""Audit official clean-map and true-position radiometer products."""
import json
from pathlib import Path
import sys
import numpy as np
from scipy.io import loadmat
from scipy.special import sph_harm_y
import scipy.linalg as sla

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from leakage.cartography_runtime import load_cartography
from leakage.datasets import sha256


def main():
    module, provenance=load_cartography(ROOT)
    products=ROOT/'data/products/meerkat_maps/cw_ra6h_dec-45_f2p5T'
    run=json.loads((products/'run_metadata.json').read_text())
    assert run['commit']==provenance['commit']
    assert run['compatibility_patch_sha256']==provenance['compatibility_patch_sha256']
    assert run['injection_metadata_sha256']==sha256(ROOT/'data/injections/cw/ra6h_dec-45_f2p5T/injection.json')
    assert run['fCW_T']==2.5 and run['true_cw_ra_hours']==6 and run['true_cw_dec_deg']==-45
    ell=run['lmax'];nkeep=run['n_keep']
    lm=[(l,m) for l in range(ell+1) for m in range(-l,l+1)]
    gamma=np.load(products/'Gamma_lm.npy')
    results=[]
    for k in (1,2,3,4):
        p=products/f'fbin_{k:02d}'
        arrays={n:np.load(p/f'{n}.npy') for n in ('M','Mprime_inv','covariance','X','clms')}
        pairs=np.load(p/'pfos_pairs.npz');sk=float(pairs['Sk'])
        R=gamma.T;N=pairs['covk']/sk**2;rho=pairs['rhok']/sk
        # Check normalisation, pair order and the actual covariance passed to official MAPS.
        np.testing.assert_allclose(arrays['M'],R.conj().T@sla.solve(N,R),rtol=1e-10,atol=1e-11)
        np.testing.assert_allclose(arrays['X'],R.conj().T@sla.solve(N,rho),rtol=1e-10,atol=1e-10)
        eig,Q=sla.eigh(arrays['M']);reference=(Q[:,-nkeep:]/eig[-nkeep:])@Q[:,-nkeep:].conj().T
        inverse_error=np.linalg.norm(reference-arrays['Mprime_inv'])/np.linalg.norm(reference)
        assert inverse_error<1e-11
        np.testing.assert_allclose(arrays['clms'],reference@arrays['X'],rtol=1e-10,atol=1e-10)
        mat=loadmat(p/'official_pixel_maps.mat',squeeze_me=True)
        # Only sparse evaluation points, not an alternative production map generator.
        indices=np.arange(113,65160,997)
        phi=mat['RA_hours'][indices]*np.pi/12;theta=np.deg2rad(90-mat['DEC_deg'][indices])
        expected=sum(c*sph_harm_y(l,m,theta,phi) for c,(l,m) in zip(arrays['clms'],lm))
        conversion_error=np.max(np.abs(expected-mat['clean_normalized'][indices]))/max(np.max(np.abs(expected)),1e-300)
        assert conversion_error<1e-11
        pixel=np.load(p/'pixel_power_maps.npz')
        np.testing.assert_array_equal(pixel['ra_hours'],(mat['RA_hours']+12)%24)
        np.testing.assert_array_equal(pixel['dec_deg'],-mat['DEC_deg'])
        np.testing.assert_allclose(pixel['clean_power_s2'],sk*mat['clean_normalized'],rtol=1e-14,atol=0)
        target=int(pixel['true_cw_pixel'])
        assert pixel['ra_hours'][target]==6 and pixel['dec_deg'][target]==-45
        u=mat['target_U']
        target_fisher=float(np.real(u@mat['M']@u.conj()))
        target_radio=float(sk*np.real((u@mat['X'])/target_fisher))
        np.testing.assert_allclose(target_fisher,mat['target_diag_fisher'],rtol=1e-13,atol=0)
        np.testing.assert_allclose(pixel['radiometer_power_at_true_cw_s2'],target_radio,rtol=1e-13,atol=0)
        assert pixel['clean_power_s2'].min()<0
        results.append({'k':k,'inverse_relative_error':float(inverse_error),
                        'matlab_vs_scipy_sparse_relative_error':float(conversion_error),
                        'radiometer_formula_matches':True,
                        'negative_power_preserved':True,'pair_covariance_positive':bool(sla.eigvalsh(pairs['covk']).min()>0)})
    # Check the coordinate issue independently of map peak location.
    from leakage.datasets import load_pulsar
    from leakage.cw_response import antenna_patterns
    base=ROOT/'data/baseline/ideal_25psr'
    psrs=[load_pulsar(base/'par'/f'J{i:02d}.par',base/'tim'/f'J{i:02d}.tim') for i in range(25)]
    pt=np.array([np.pi/2-p.model.DECJ.quantity.to_value('rad') for p in psrs])
    pp=np.array([p.model.RAJ.quantity.to_value('rad') for p in psrs])
    fplus,fcross=antenna_patterns(psrs,3*np.pi/4,np.pi/2)
    # Evaluate upstream at the antipode of the injection source.
    official=module.ac.signalResponse_fast(pt,pp,np.array([np.pi/4]),np.array([3*np.pi/2]))
    correlation=official@official.T
    expected=1.5*(np.outer(fplus,fplus)+np.outer(fcross,fcross))
    direction_error=np.linalg.norm(correlation-expected)/np.linalg.norm(expected)
    assert direction_error<1e-12
    # One quadrature check for this fixed dataset, not an injection parameter scan.
    finer=module.ac.anis_basis(np.column_stack((pp,pt)),lmax=ell,nside=32,sph_complex=True)
    a,b=np.triu_indices(25,1)
    finer=finer[:,a,b]
    quadrature_error=np.linalg.norm(finer-gamma)/np.linalg.norm(finer)
    assert quadrature_error<0.01, 'NSIDE=16 angular integration differs by over 1%; review the notebook'
    result={**provenance,'checks':results,'source_convention_relative_error':float(direction_error),
            'Gamma_lm_nside16_vs32_relative_error':float(quadrature_error)}
    out=ROOT/'validation/cartography_f2p5T';out.mkdir(exist_ok=True,parents=True)
    (out/'numerical_checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':
    from pint.logging import setup
    setup(level='ERROR')
    main()
