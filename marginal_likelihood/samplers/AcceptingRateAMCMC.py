import numpy as np
from marginal_likelihood.samplers.MetropolisHastings import \
        MetropolisHastings
from marginal_likelihood.LikelihoodFunction import LikelihoodFunction
from utils import safe_pow_exp_ratio
from utils import safe_exp
from utils import safe_exp_ratio
from distributions.MultivariateLognormal import MultivariateLognormal
from distributions.MultivariateNormal import MultivariateNormal
from distributions.MultivariatePositiveNormal import MultivariatePositiveNormal

class AcceptingRateAMCMC (MetropolisHastings):
    """ This class is able to return a sample of theta using an adaptive
        algorithm that samples each parameter independently and with
        a log-normal proposal distribution with adapting variance. The
        variance is updated according to the acceptance ratio of 
        proposals. This method is inspired in the one described on 
        "Supplementary Materials for Inferring Signaling Pathway 
        Topologies from Multiple Perturbation Measurements of Specific 
        Biochemical Species", Tian-Rui Xu et. al. """

    def __init__ (self, theta, model, experiments, sigma_update_n, 
            verbose=False, t=1):
        """ Default constructor. """
        super ().__init__ (theta, verbose=verbose)
        self.__model = model
        self.__experiments = experiments
        self.__sigma_update_n = sigma_update_n
        self._jump_S = self.__init_jump_S ()
        self.__l_f = LikelihoodFunction (model)
        self.__t = t


    def __init_jump_S (self):
        """ Gives an initial value for the proposal distribution sigma
            on phase one of inferece. """
        theta = self._theta
        jump_S = []
        for p in theta:
            param_dist = p.get_distribution ()
            prior_variance = param_dist.variance ()
            sigma2 = prior_variance
            jump_S.append (sigma2)
        print ("Starting jump S: ")
        print (jump_S)
        return jump_S

    
    def _create_jump_dist (self, theta_t):
        """ The jump distribution is Multivariate Lognormal with a 
            diagonal covariance matrix, i.e the jumps on each parameter
            are independent. """
        n = theta_t.get_size ()
        t_vals = np.array (theta_t.get_values ())
        S = np.eye (n)
        for i in range (n):
            S[i, i] = self._jump_S[i]
        mu = t_vals
        print ("Creating jump dist with mean: ")
        print (t_vals)
        print ("and variance: ")
        print (S)
        jump_dist = MultivariatePositiveNormal (mu, S)
        return jump_dist


    def set_temperature (self, t):
        """ Defines a temperature parameter for this sampler. """
        self.__t = t


    def _calc_mh_ratio (self, new_t, new_l, old_t, old_l):
        """ In this case, the MH ratio should be:
            [p (y | t*) / p (y | t)] * [p (t*) / p(t)] * 
            [J (t | t*) / J (t* | t)] 
            where t is the current parameter, t* is the proposed 
            parameter, and y is the observations from the experiment. 
        """
        j_gv_old = self._create_jump_dist (old_t)
        j_gv_new = self._create_jump_dist (new_t)
        new_gv_old = j_gv_old.pdf (new_t.get_values ())
        old_gv_new = j_gv_new.pdf (old_t.get_values ())
        l_ratio = safe_pow_exp_ratio (new_l, old_l, self.__t)
        prior_ratio = safe_exp_ratio (new_t.get_log_p (), 
                old_t.get_log_p ())
        jump_ratio = old_gv_new / new_gv_old
        if self._is_verbose:
            print ("\told log prior: " + str (old_t.get_log_p ()))
            print ("\tnew log prior: " + str (new_t.get_log_p ()))
            print ("\tnew given old: " + str (new_gv_old))
            print ("\told given new: " + str (old_gv_new))
            print ("\tprior ratio: " + str (prior_ratio))
            print ("\tlikelihood ratio: " + str (l_ratio))
            print ("\tjump ratio: " + str (jump_ratio))
        return l_ratio * prior_ratio * jump_ratio
        

    def _calc_log_likelihood (self, theta):
        """ Calculates the log of p (experiments | theta, model). """
        return self.__l_f.get_log_likelihood (self.__experiments, theta)


    def _iteration_update (self):
        if self._n_jumps % self.__sigma_update_n == 0:
            self.__update_Sigma ()
        
    
    def __update_Sigma (self):
        """ Updates jump_S according to current acceptance rate. We do
            so as recommended in section "Efficient Metopolis jumping 
            rules" from Bayesian Data Analysis (Third Edition), Gelman. 
            """
        print ("Sigma started with: ")
        acceptance_rate = self.get_acceptance_ratio ()
        jump_S = self._jump_S
        print (jump_S)
        for i in range (len (jump_S)):
            if acceptance_rate > .4 and jump_S[i] < 10:
                jump_S[i] += jump_S[i] * .5
            if acceptance_rate < .25 and jump_S[i] > 1e-4:
                jump_S[i] -= jump_S[i] * .5
        print ("Updated sigma to: ")
        print (jump_S)
