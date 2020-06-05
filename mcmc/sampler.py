"""
Running and output manipulation for tlemcee MCMC algorithm
"""

import sys
import emcee
import pickle
import warnings
import numpy as np
from time import time

from st.tle import TLE
from mcmc.config import Config
from mcmc.utils import (
    lnprob,
    getBestModel
)
from mcmc.diagnostics import (
    plotTimeSeries,
    plotCorner
)

class Sampler:
    """
    Container for MCMC sampler
    """
    def __init__(self, config):
        """
        Initialise Sampler

        Parameters
        ----------
        config : mcmc.config.Config
            Config containing prior and data information
        """
        if self._validate(config):
            self._setup()
        else:
            sys.exit()

    def _validate(self, config):
        """
        Validate input

        Parameters
        ----------
        config : mcmc.config.Config
            Config containing prior and data information
        """
        self.config = config
        if not isinstance(self.config, Config):
            print('SamplerError: Input must be Config')
            return 0

        if self.config.out_dir is not None:
            self.out_dir = self.config.out_dir
        else:
            self.out_dir = '.'
        return 1

    def _setup(self):
        """
        Set up sampler
        """
        self.sampler = emcee.EnsembleSampler(self.config.n_walkers,
                                             self.config.n_dim,
                                             lnprob,
                                             args=(self.config,))
        self._run = 0

    def run(self, save=False, verbose=False):
        """
        Run MCMC

        Parameters
        ----------
        save : bool, optional
            Toggle to save sampler chain to file
            Default = False
        verbose : bool, optional
            Toggle to print run time
            Default = False
        """
        t0 = time()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            self.sampler.run_mcmc(self.config.init_pos,
                                  self.config.n_steps,
                                  rstate0=np.random.get_state())
        t1 = time()
        self._run = 1

        if save:
            out_path = '{}/sampler_chain_{}_{}.pkl'.format(self.out_dir,
                                                           self.config.n_steps,
                                                           self.config.n_walkers)
            with open(out_path, 'wb') as f:
                pickle.dump(self.sampler.chain, f)

        if verbose:
            print('MCMC completed in {} seconds'.format(t1 - t0))

    def get_best_model(self, burnin=None, save=False):
        """
        Obtain best model from sampler chain

        Parameters
        ----------
        burnin : int, optional
            Burn-in phase for sampler chain (user-assessed)
            If given, use median method on sampler chain with burn-in phase removed
            If None, use argmax method on entire sampler chain
            Default = None
        save : bool, optional
            Toggle to save best model TLE to file
            Default = False

        Returns
        -------
        best_model : st.tle.TLE
            Best model TLE
        """
        return getBestModel(self,
                            burnin=burnin,
                            save=save)

    def plot_time_series(self, save=False):
        """
        Plot time series for parameters

        Parameters
        ----------
        save : bool, optional
            Toggle to save time series plot to file
            Default = False
        """
        plotTimeSeries(self,
                       save=save)

    def plot_corner(self, burnin=0, save=False):
        """
        Plot corner plot for sampler

        Parameters
        ----------
        burnin : int, optional
            Burn-in phase for sampler chain (user-assessed)
            Default = 0 (no burn-in)
        save : bool, optional
            Toggle to save corner plot to file
            Default = False
        """
        plotCorner(self,
                   burnin=burnin,
                   save=save)
