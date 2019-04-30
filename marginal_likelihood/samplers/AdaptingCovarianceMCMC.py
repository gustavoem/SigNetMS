import numpy as np
from marginal_likelihood.samplers.MetropolisHastings import \
        MetropolisHastings
from marginal_likelihood.LikelihoodFunction import LikelihoodFunction
from distributions.MultivariateLognormal import MultivariateLognormal
from distributions.MultivariateNormal import MultivariateNormal
from distributions.MultivariatePositiveNormal import MultivariatePositiveNormal
from marginal_likelihood.CovarianceMatrix import calc_covariance
from utils import safe_power
from utils import safe_exp
from utils import safe_exp_ratio
from utils import safe_pow_exp_ratio

class AdaptingCovarianceMCMC (MetropolisHastings):
    """ Objects of this class are able to return a sample of theta using 
        an adaptive metropolis-hastings algorithm. The proposal 
        distribution of the algorithm is adaptive and has as covariance 
        matrix an estimate of the covariance according to the current 
        sample. To get a sample with an instance of this class, you will 
        need a starting sample to provide a first estimate of the 
        covariance matrix. """

    def __init__ (self, theta, model, experiments, t=1, verbose=False):
        """ Default constructor. """
        super ().__init__ (theta, verbose=verbose)
        self.__model = model
        self.__experiments = experiments
        self._jump_S = None
        self.__l_f = LikelihoodFunction (model)
        self._t = t


    def set_temperature (self, t):
        """ Defines a temperature parameter. """
        self._t = t


    def __calc_jump_S (self):
        """ Calculates jump_S, an estimate of the covariance of 
            parameters. """
        sample_values = [] 
        for t in self._sample:
            sample_values.append (t.get_values ())
        self._jump_S = calc_covariance (sample_values) 
    
    
    def get_jump_covariance (self):
        """ Returns the jump_S matrix. """
        return np.array (self._jump_S)


    def _create_jump_dist (self, theta_t):
        """ The jump distribution is Multivariate Lognormal. """
        t_vals = theta_t.get_values ()
        mu = np.array (t_vals)
        print ("Creating jump distribution with mean: ", end='')
        print (mu)
        print ("and variance: ")
        print (self._jump_S)
        jump_dist = MultivariatePositiveNormal (mu, self._jump_S)
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
        prior_ratio = safe_exp_ratio (new_t.get_log_p (), 
                old_t.get_log_p ())
        jump_ratio = old_gv_new / new_gv_old
        if self._is_verbose:
            print ("\tnew given old: " + str (new_gv_old))
            print ("\told given new: " + str (old_gv_new))
            print ("\tprior ratio: " + str (prior_ratio))
            print ("\tlikelihood ratio: " + str (l_ratio))
            print ("\tjump ratio: " + str (jump_ratio))
        return l_ratio * prior_ratio * jump_ratio


    def _calc_likeli_ratio (self, log_new_l, log_old_l):
        """ Calcultes the ratio (new_l / old_l) ^ t. """
        return safe_pow_exp_ratio (log_new_l, log_old_l, self._t)


    def _calc_log_likelihood (self, theta):
        """ Calculates the log of p (experiments | theta, model). """
        return self.__l_f.get_log_likelihood (self.__experiments, theta)


    def _iteration_update (self):
        """ At the end of each sampling iteration, we should update the
            Covariance Matrix. """
        self.__calc_jump_S ()
