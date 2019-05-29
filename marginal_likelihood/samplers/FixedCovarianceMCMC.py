from marginal_likelihood.samplers.MetropolisHastings import \
        MetropolisHastings
from marginal_likelihood.samplers.AdaptingCovarianceMCMC import \
        AdaptingCovarianceMCMC

class FixedCovarianceMCMC (AdaptingCovarianceMCMC):
    """ Objects of this class are able to return a sample of theta using
        the Metropolis-Hastings algorithm. The proposal distribution is 
        Multivariate Lognormal and it's shape is defined on the 
        constructor through a covariance matrix."""

    def __init__ (self, theta, model, experiments, covar, t=1, 
            verbose=False):
        """ Default constructor. """
        super ().__init__ (theta, model, experiments, t=t, 
                verbose=verbose)
        self._jump_S = covar 


    def get_sample (self, N):
        return MetropolisHastings.get_sample (self, N)


    def _iteration_update (self):
        """ There's no structure to be updated on each iteration. """
        pass
