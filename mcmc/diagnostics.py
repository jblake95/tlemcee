"""
Plotting for tlemcee MCMC algorithm
"""

import corner
import matplotlib.pyplot as plt

def plotTimeSeries(sampler, save=False):
    """
    Plot time series for parameters

    Parameters
    ----------
    sampler : mcmc.sampler.Sampler
        Sampler containing chain information
    save : bool, optional
        Toggle to save time series plot to file
        Default = False
    """
    fig, ax = plt.subplots(sampler.config.n_dim,
                           figsize=(6, 3 * sampler.config.n_dim),
                           sharex=True)
    for i, (name, init) in enumerate(zip(sampler.config.priors['name'],
                                         sampler.config.priors['init'])):
        ax[i].plot(sampler.sampler.chain[:, :, i].T, color='k', alpha=0.4)
        ax[i].axhline(init, color='#888888', lw=2)

        ax[i].set_ylabel(name)

        if i == sampler.config.n_dim - 1:
            ax[i].set_xlabel('Step number')

    if save:
        out_path = '{}/time_series_{}_{}.pdf'.format(sampler.out_dir,
                                                     sampler.config.n_steps,
                                                     sampler.config.n_walkers)
        plt.savefig(out_path,
                    overwrite=True,
                    bbox_inches='tight',
                    pad_inches=0.01)
    plt.show()

def plotCorner(sampler, burnin=0, save=False):
    """
    Plot corner plot for sampler

    Parameters
    ----------
    sampler : mcmc.sampler.Sampler
        Sampler containing chain information
    burnin : int, optional
        Burn-in phase for sampler chain (user-assessed)
        Default = 0 (no burn-in)
    save : bool, optional
        Toggle to save corner plot to file
        Default = False
    """
    samples = sampler.sampler.chain[:, burnin:, :].reshape((-1,
                                                            sampler.config.n_dim))
    fig = corner.corner(samples,
                        labels=sampler.config.priors['name'],
                        plot_contours=True)
    if save:
        out_path = '{}/corner_{}_{}.pdf'.format(sampler.out_dir,
                                                sampler.config.n_steps,
                                                sampler.config.n_walkers)
        plt.savefig(out_path,
                    overwrite=True,
                    bbox_inches='tight',
                    pad_inches=0.01)
    fig.show()
