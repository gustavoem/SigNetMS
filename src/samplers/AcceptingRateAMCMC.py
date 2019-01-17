import numpy as np
from MetropolisHastings import MetropolisHastings
from LikelihoodFunction import LikelihoodFunction
from MultivariateLognormal import MultivariateLognormal

class AcceptingRateAMCMC (MetropolisHastings):
    """ This class is able to return a sample of theta using an adaptive
        algorithm that samples each parameter independently and with
        a log-normal proposal distribution with adapting variance. The
        variance is updated according to the acceptance ratio of 
        proposals. This method is inspired in the one described on 
        "Supplementary Materials for Inferring Signaling Pathway 
        Topologies from Multiple Perturbation Measurements of Specific 
        Biochemical Species", Tian-Rui Xu et. al. """

    def __init__ (self, theta, model, experiments, sigma_update_n):
        """ Default constructor. """
        super ().__init__ (theta, verbose=False)
        self.__model = model
        self.__experiments = experiments
        self.__sigma_update_n = sigma_update_n
        self._jump_S = self.__init_jump_S ()
        self.__l_f = LikelihoodFunction (model)


    def __init_jump_S (self):
        """ Gives an initial value for the proposal distribution sigma
            on phase one of inferece. """
        theta = self._theta
        jump_S = []
        for p in theta:
            param_dist = p.get_distribution ()
            prior_variance = param_dist.variance ()
            sigma2 = np.log (np.sqrt (prior_variance) + 1)
            jump_S.append (np.sqrt (sigma2))
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
        mu = np.log (t_vals) - S.diagonal () / 2
        jump_dist = MultivariateLognormal (mu, S)
        return jump_dist


    def _calc_mh_ratio (self, new_t, new_l, old_t, old_l):
        """ In this case, the MH ratio should be:
            [p (y | t*) / p (y | t)] * [p (t*) / p(t)] * 
            [J (t | t*) / J (t* | t)] 
            where t is the current parameter, t* is the proposed 
            parameter, and y is the observations from the experiment. 
        """
        j_gv_old = self._create_jump_dist (old_t)
        j_gv_new = self._create_jump_dist (new_t)
        new_gv_old = j_gv_old.pdf (new_t)
        old_gv_new = j_gv_new.pdf (old_t)

        l_ratio = np.exp (new_l - old_l)
        prior_ratio = new_t.get_p () / old_t.get_p ()
        jump_ratio = old_gv_new / new_gv_old
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
        acceptance_rate = self.get_acceptance_ratio ()
        jump_S = self._jump_S
        print ("alo")
        for i in range (len (jump_S)):
            if acceptance_rate > .4 and jump_S[i] < 10:
                jump_S[i] += jump_S[i] * .5
            if acceptance_rate < .25 and jump_S[i] > 1e-4:
                jump_S[i] -= jump_S[i] * .5
