import numpy as np

import histlite as hl
import csky as cy

timer = cy.timing.Timer()
time = timer.time

ana_dir = cy.utils.ensure_dir('/data/user/apizzuto/csky_cache/')
with time('ana setup'):
    ana = cy.get_analysis(
        cy.selections.repo,
        'version-002-p05', cy.selections.GFUDataSpecs.GFU_IC86,
        dir=ana_dir)

conf = {'extended': True,
       'space': "ps",
        'time': "transient",
        'sig': 'transient',
       }
cy.CONF['mp_cpus'] = 5
sinDecs = [-0.5, 0.0, 0.5]
inj_gamma = [2.0, 2.5, 3.0]
delta_ts = np.logspace(3., 7., 9)[::]

all_gamma_results = dict()

for gamma in inj_gamma[:]:
    all_gamma_results[gamma] = dict()
    for sd in sinDecs[:]:
        results = dict(sens_n = [],
                  sens_e2dnde = [],
                  disc_n = [],
                  disc_e2dnde = [],
                  delta_t = [])
        dec = np.arcsin(sd)
        for delta_t in delta_ts:
            src = cy.utils.Sources(ra=0.0,
                                   dec=dec,
                                   mjd=57000.,
                                   sigma_t=0.,
                                   t_100=delta_t/86400.)
            cy.CONF['src'] = src

            tr = cy.get_trial_runner(conf, ana=ana, src=src)
            bg = cy.dists.Chi2TSD(tr.get_many_fits(10000))
            
            tr = cy.get_trial_runner(conf, ana=ana, src=src,
                                     inj_conf={'flux': cy.hyp.PowerLawFlux(gamma)})

            sensitivity = tr.find_n_sig(bg.median(), 0.9,
                                   batch_size=1000,
                                   n_sig_step=0.25,
                                   #max_batch_size=0,
                                   logging=False,
                                   n_bootstrap=1, tol=.05)

            thresh_ts = bg.isf_nsigma(5.)
            discovery = tr.find_n_sig(thresh_ts, 0.5,
                                   batch_size=1000,
                                   n_sig_step=0.25,
                                   #max_batch_size=0,
                                   logging=False,
                                   n_bootstrap=1, tol=.05)
            
            results['sens_n'].append(sensitivity['n_sig'])
            results['disc_n'].append(discovery['n_sig'])
            results['sens_e2dnde'].append(tr.to_E2dNdE(sensitivity, E0=100, unit=1e3))
            results['disc_e2dnde'].append(tr.to_E2dNdE(discovery['n_sig'], E0=100, unit=1e3))
            results['delta_t'].append(delta_t)
            
        all_gamma_results[gamma][sd] = results

def conv_ref_en(fl, e1, e0, gamma):
    return np.asarray(fl) * (e1/e0)**(2.-gamma)

import pickle
with open('./gfu_offline_transient_sens_100tev.pkl', 'wb') as f:
    pickle.dump(all_gamma_results, f)

class GFUOnlineDataSpecs(object):
    class GFUOnlineDataSpec(cy.selections.TrackSpec):
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

    class GFUOnline_IC86 (GFUOnlineDataSpec):
        _path_sig = 'gfu_online/{version}/IC86_2011_MC.npy'
        _path_data = ['gfu_online/{{version}}/IC86_201{}_data.npy'.format(i) for i in range(1,9)]
        _path_grls = ['gfu_online/{{version}}/GRL/IC86_201{}_data.npy'.format(i) for i in range(1,9)]
        def __init__(self, years=list(map(str, 2010 + np.arange(1, 9)))):
            self.path_data = ['gfu_online/{{version}}/IC86_{}_data.npy'.format(y) for y in years]
            self.path_grls = ['gfu_online/{{version}}/GRL/IC86_{}_data.npy'.format(y) for y in years]
            self._key = 'GFUOnline_for_' + '_'.join(years)

    gfuonline_IC86 = [GFUOnline_IC86]

ana_dir = cy.utils.ensure_dir('/data/user/apizzuto/csky_cache/')
with time('ana setup'):
    ana = cy.analysis.Analysis(cy.selections.repo, [GFUOnlineDataSpecs.GFUOnline_IC86])

conf = {'extended': True,
       'space': "ps",
        'time': "transient",
        'sig': 'transient',
       }
cy.CONF['mp_cpus'] = 5
sinDecs = [-0.5, 0.0, 0.5]
inj_gamma = [2.0, 2.5, 3.0]
delta_ts = np.logspace(3., 7., 9)

all_gamma_results_online = dict()

for gamma in inj_gamma:
    all_gamma_results_online[gamma] = dict()
    for sd in sinDecs:
        results = dict(sens_n = [],
                  sens_e2dnde = [],
                  disc_n = [],
                  disc_e2dnde = [],
                  delta_t = [])
        dec = np.arcsin(sd)
        for delta_t in delta_ts:
            src = cy.utils.Sources(ra=0.0,
                                   dec=dec,
                                   mjd=57000.,
                                   sigma_t=0.,
                                   t_100=delta_t/86400.)
            cy.CONF['src'] = src

            tr = cy.get_trial_runner(conf, ana=ana, src=src)
            bg = cy.dists.Chi2TSD(tr.get_many_fits(10000))

            tr = cy.get_trial_runner(conf, ana=ana, src=src,
                                     inj_conf={'flux': cy.hyp.PowerLawFlux(gamma)})

            sensitivity = tr.find_n_sig(bg.median(), 0.9,
                                   batch_size=1000,
                                   n_sig_step=0.25,
                                   #max_batch_size=0,
                                   logging=False,
                                   n_bootstrap=1, tol=0.05)

            thresh_ts = bg.isf_nsigma(5.)
            discovery = tr.find_n_sig(thresh_ts, 0.5,
                                   batch_size=1000,
                                   n_sig_step=0.25,
                                   #max_batch_size=0,
                                   logging=False,
                                   n_bootstrap=1, tol=0.05)
            
            results['sens_n'].append(sensitivity['n_sig'])
            results['disc_n'].append(discovery['n_sig'])
            results['sens_e2dnde'].append(tr.to_E2dNdE(sensitivity, E0=100, unit=1e3))
            results['disc_e2dnde'].append(tr.to_E2dNdE(discovery['n_sig'], E0=100, unit=1e3))
            results['delta_t'].append(delta_t)
            
        all_gamma_results_online[gamma][sd] = results

with open('./gfu_online_transient_sens_100tev.pkl', 'wb') as f:
    pickle.dump(all_gamma_results_online, f)