import numpy as np
import matplotlib as mpl
mpl.use('agg')

import histlite as hl
import csky as cy

timer = cy.timing.Timer()
time = timer.time

one_year = True
if one_year:
    time_str = 'one_year'
    class GFUOneYearDataSpecs(object):
        class GFUOneYearDataSpec(cy.selections.TrackSpec):
            _bins_sindec = np.unique(np.concatenate([
                np.linspace(-1, -0.93, 4 + 1),
                np.linspace(-0.93, -0.3, 10 + 1),
                np.linspace(-0.3, 0.05, 9 + 1),
                np.linspace(0.05, 1, 18 + 1) ]))
            _bins_logenergy = np.arange(1, 9.5 + 0.01, 0.125)
            def dataset_modifications(self, ds):
                max_sigma = np.radians(15)
                ds.data = ds.data[ds.data.sigma < max_sigma]
                ds.sig = ds.sig[ds.sig.sigma < max_sigma]
                #ds.data.sigma = np.minimum(ds.data.sigma, max_sigma)
                #ds.sig.sigma = np.minimum(ds.sig.sigma, max_sigma)

        class GFUOneYear (GFUOneYearDataSpec):
            _path_sig = 'gfu/{version}/IC86_2011_MC.npy'
            _path_data = ['gfu/{{version}}/IC86_201{}_data.npy'.format(i) for i in range(7,8)]
            _path_grls = ['gfu/{{version}}/GRL/IC86_201{}_data.npy'.format(i) for i in range(7,8)]
            def __init__(self, years=list(map(str, 2010 + np.arange(7,8)))):
                self.path_data = ['gfu/{{version}}/IC86_{}_data.npy'.format(y) for y in years]
                self.path_grls = ['gfu/{{version}}/GRL/IC86_{}_data.npy'.format(y) for y in years]
                self._key = 'GFU_for_' + '_'.join(years)

        gfu_oneyear = [GFUOneYear]
    with time('ana setup'):
        ana = cy.analysis.Analysis(cy.selections.repo, [GFUOneYearDataSpecs.GFUOneYear])
else:
    time_str = 'three_years'
    ana_dir = cy.utils.ensure_dir('/data/user/apizzuto/csky_cache/')
    with time('ana setup'):
        ana = cy.get_analysis(
            cy.selections.repo,
            'version-002-p05', cy.selections.GFUDataSpecs.GFU_2015_2017,
            dir=ana_dir)

#gfu_3yr = [GFU_2015_2017]

cy.CONF['ana'] = ana
cy.CONF['mp_cpus'] = 10

sinDecs = np.r_[-0.95:0.95:30j]
inj_gamma = [2.0]
all_gamma_results = dict()

for gamma in inj_gamma:
    results = dict(sens_n = [],
                  sens_e2dnde = [],
                  disc_n = [],
                  disc_e2dnde = [],
                  sindec = [])

    for sd in sinDecs:
        print(f"\nWorking on source at sindec {sd}")

        src = cy.sources(0.0, np.arcsin(sd))
        cy.CONF['src'] = src
        tr = cy.get_trial_runner(flux=cy.hyp.PowerLawFlux(gamma=gamma))
        bg = cy.dists.Chi2TSD(tr.get_many_fits(20000, seed=1))
        sens = tr.find_n_sig(
            # ts, threshold
            bg.median(),
            # beta, fraction of trials which should exceed the threshold
            0.9,
            # n_inj step size for initial scan
            n_sig_step=3,
            # this many trials at a time
            batch_size=800,
            # tolerance, as estimated relative error
            tol=.02
        )

        disc = tr.find_n_sig(bg.isf_nsigma(5), 0.5, n_sig_step=5, batch_size=800, tol=.02)
        results['sens_n'].append(sens['n_sig'])
        results['disc_n'].append(disc['n_sig'])
        results['sens_e2dnde'].append(tr.to_E2dNdE(sens, E0=100, unit=1e3))
        results['disc_e2dnde'].append(tr.to_E2dNdE(disc['n_sig'], E0=100, unit=1e3))
        results['sindec'].append(sd)
    all_gamma_results[gamma] = results

def conv_ref_en(fl, e1, e0, gamma):
    return np.asarray(fl) * (e1/e0)**(2.-gamma)

import pickle
with open('./dump/gfu_offline_sens_100tev_ovro_proposal_{}.pkl'.format(time_str), 'wb') as f:
    pickle.dump(all_gamma_results, f)

