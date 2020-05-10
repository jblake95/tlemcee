"""
Running and output manipulation for tlemcee MCMC algorithm
"""

import sys
import emcee
import pickle
import numpy as np
from time import time

from mcmc.utils import lnprob
from mcmc.config import Config

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
        return 1

    def _setup(self):
        """
        Set up sampler
        """
        self.sampler = emcee.EnsembleSampler(self.config.n_walkers,
                                             self.config.n_dim,
                                             lnprob,
                                             args=(self.config,))

    def _save(self):
        """
        Save sampler chain to file
        """
        if self.config.out_dir is not None:
            out_dir = self.config.out_dir
        else:
            out_dir = '.'
        out_path = '{}/sampler_chain_{}_{}.pkl'.format(out_dir,
                                                       self.config.n_steps,
                                                       self.config.n_walkers)
        with open(out_path, 'wb') as f:
            pickle.dump(self.sampler.chain, f)

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
        self.sampler.run_mcmc(self.config.init_pos,
                              self.config.n_steps,
                              rstate0=np.random.get_state())
        t1 = time()
        if verbose:
            print('MCMC completed in {} seconds'.format(t1 - t0))
        if save:
            self._save()
