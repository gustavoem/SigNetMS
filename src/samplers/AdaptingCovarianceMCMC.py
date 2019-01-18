import numpy as np
from MetropolisHastings import MetropolisHastings
from LikelihoodFunction import LikelihoodFunction
from MultivariateLognormal import MultivariateLognormal
from CovarianceMatrix import calc_covariance
from utils import safe_power

class AdaptingCovarianceMCMC (MetropolisHastings):
    """ This class is able to return a sample of theta using an adaptive
        metropolis-hastings algorithm. The proposal distribution of this
        algorithm is adaptive and has as covariance matrix an estimate
        of the covariance according to the current sample. To get a 
        sample with an instance of this class, you will need a starting 
        sample to have a first estimate of the covariance matrix. """

    def __init__ (self, theta, model, experiments, t=1, verbose=False):
        """ Default constructor. """
        super ().__init__ (theta, verbose=verbose)
        self.__model = model
        self.__experiments = experiments
        self._jump_S = None
        self.__l_f = LikelihoodFunction (model)
        self._t = 1


    def __calc_jump_S (self):
        """ Calculates jump_S, an estimate of the covariance of 
            parameters. """
        sample_values = [] 
        for t in self._sample:
            sample_values.append (t.get_values ())
        self._jump_S = calc_covariance (sample_values) 


    def _create_jump_dist (self, theta_t):
        """ The jump distribution is Multivariate Lognormal. """
        t_vals = theta_t.get_values ()
        mu = np.log (t_vals) - self._jump_S.diagonal () / 2
        jump_dist = MultivariateLognormal (mu, self._jump_S)
        return jump_dist

    def get_sample (self, N):
        """ Get a sample of size N. """
        if len (self._sample) == 0:
            raise ValueError ("The current sample can't be empty. " \
                    + "Try using the start_sample_from_prior () " \
                    + "method.")

        self.__calc_jump_S ()
        return super ().get_sample (N)

    def _calc_mh_ratio (self, new_t, new_l, old_t, old_l):
        """ In this case, the MH ratio should be:
            [p (y | t*) / p (y | t)] * [p (t*) / p(t)] * 
            [J (t | t*) / J (t* | t)] 
            where t is the current parameter, t* is the proposed 
            parameter, and y is the observations from the experiment. 
        """
        j_gv_old = self._create_jump_dist (old_t)
        j_gv_new = self._create_jump_dist (new_t)
        try:
            new_gv_old = j_gv_old.pdf (new_t.get_values ())
            old_gv_new = j_gv_new.pdf (old_t.get_values ())
        except:
            raise ValueError ("The covariance matrix is not positive " \
                    + "definite. Try using a bigger starting sample.")

        l_ratio = self._calc_likeli_ratio (new_l, old_l)
        prior_ratio = new_t.get_p () / old_t.get_p ()
        jump_ratio = old_gv_new / new_gv_old
        return l_ratio * prior_ratio * jump_ratio


    def _calc_likeli_ratio (self, log_new_l, log_old_l):
        """ Calcultes the ratio (new_l / old_l) ^ t. """
        return safe_power (np.exp (log_new_l - log_old_l), self._t)


    def _calc_log_likelihood (self, theta):
        """ Calculates the log of p (experiments | theta, model). """
        return self.__l_f.get_log_likelihood (self.__experiments, theta)


    def _iteration_update (self):
        """ At the end of each sampling iteration, we should update the
            Covariance Matrix. """
        self.__calc_jump_S ()
