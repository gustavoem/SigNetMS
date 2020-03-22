from marginal_likelihood.samplers.MetropolisHastings import \
        MetropolisHastings
from marginal_likelihood.samplers.AdaptingCovarianceMCMC import \
        AdaptingCovarianceMCMC
from utils import get_current_datetime

class FixedCovarianceMCMC (AdaptingCovarianceMCMC):
    """ Objects of this class are able to return a sample of theta using
        the Metropolis-Hastings algorithm. The proposal distribution is 
        Multivariate Lognormal and it's shape is defined on the 
        constructor through a covariance matrix."""

    def __init__ (self, theta, model, experiments, covar, t=1, 
            verbose=False):
        """ Default constructor. """
        super ().__init__ (theta, model, experiments, 1e10, t=t, 
                verbose=verbose)
        self._jump_S = covar 


    def _open_trace_file (self):
        """ Open a file to write trace. """
        file_name = "trace/" + get_current_datetime () + "_" \
                + str (self._t) + "_" + "3rd_phase"
        self._trace_file = open (file_name, 'w')


    def get_sample (self, N):
        return MetropolisHastings.get_sample (self, N)


    def _iteration_update (self):
        """ There's no structure to be updated on each iteration. """
        pass
