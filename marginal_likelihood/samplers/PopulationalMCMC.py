import numpy as np
from marginal_likelihood.LikelihoodFunction import LikelihoodFunction
from utils import safe_exp_ratio
from utils import safe_pow_exp_ratio
from distributions.DiscreteLaplacian import DiscreteLaplacian


class PopulationalMCMC:

    """ Performs a tempered sampling of theta. The sampling chooses 
        temperatures and for each of them conduces a 
        FixedCovarianceMCMC. The algorithm also mixes samples of 
        different temperatures. """


    __SCHEDULE_POWER = 4


    def __init__ (self, n_strata, strata_size, fc_mcmcs, betas=None, 
            verbose=False):
        """ Default constructor. """
        if n_strata * strata_size != len (fc_mcmcs):
            raise ValueError ("The list of covariances and starts " \
                    + "should have the same size as n_strata * " \
                    + "strata_size")

        self.__n_strata = n_strata
        self.__strata_size = strata_size
        if betas == None:
            self.__betas = PopulationalMCMC.sample_scheduled_betas \
                (n_strata * strata_size)
        else:
            self.__betas = betas
        self.__fc_mcmcs = fc_mcmcs
        self.__define_samplers_temp (self.__betas, self.__fc_mcmcs)
        self.__verbose = verbose


    @staticmethod
    def get_sched_power ():
        """ Returns the beta schedule power. """
        return PopulationalMCMC.__SCHEDULE_POWER

    
    @staticmethod
    def sample_betas (n_strata, strata_size):
        """ Samples the temperatures array. """
        betas = [0]
        sched_power = PopulationalMCMC.__SCHEDULE_POWER
        for i in range (n_strata):
            strata_start = (i /  n_strata) ** sched_power
            strata_end = ((i + 1) / n_strata) ** sched_power
            for j in range (strata_size):
                x = np.random.uniform (strata_start, strata_end)
                betas.append (x)
            betas.sort ()
        betas.append (1)
        return betas


    @staticmethod
    def sample_scheduled_betas (n):
        """ Samples betas from a fixed schedule. """
        betas = []
        sched_power = PopulationalMCMC.__SCHEDULE_POWER
        for i in range (n):
            x = (i / (n - 1)) ** sched_power
            betas.append (x)
        betas.sort ()
        return betas


    def __define_samplers_temp (self, betas, fc_mcmcs):
        for i in range (len (betas)):
            fc_mcmcs[i].set_temperature (betas[i])
    

    def get_sample (self, N):
        """ Get a sample of size N. """
        betas = self.__betas
        fc_mcmcs = self.__fc_mcmcs
        for i in range (N):
            if self.__verbose:
                print (str (i) + "-th iteration of PopulationalMCMC.")
            for j in range (len (self.__betas)):
                fc_mcmcs[j].get_sample (1)

            j = np.random.choice (range (len (betas)))
            temp_jump_dist = DiscreteLaplacian (len (betas), j + 1)
            k = temp_jump_dist.rvs () - 1
            inv_temp_jump_dist = DiscreteLaplacian (len (betas), k + 1)
            
            sample, likelihoods = fc_mcmcs[j].get_last_sampled (1)
            thetaj, thetaj_l = sample[0], likelihoods[0]
            sample, likelihoods = fc_mcmcs[k].get_last_sampled (1)
            thetak, thetak_l = sample[0], likelihoods[0]
            tjotk = safe_exp_ratio (thetaj_l, thetak_l)
            tkotj = safe_exp_ratio (thetak_l, thetaj_l)
            j_gv_k = inv_temp_jump_dist.pdf (j + 1)
            k_gv_j = temp_jump_dist.pdf (k + 1)
            
            tktj_pw = safe_pow_exp_ratio (thetak_l, thetaj_l, betas[j])
            tjtk_pw = safe_pow_exp_ratio (thetaj_l, thetak_l, betas[k])
            r = tktj_pw * tjtk_pw * (j_gv_k / k_gv_j)
            
            if self.__verbose:
                print ("\ttheta_j over theta_k: " + str (tjotk))
                print ("\ttheta_k over theta_j: " + str (tkotj))
                print ("\tr: " + str (r))

            if np.random.uniform () <= r:
                if self.__verbose:
                    print ("Inverted j and k.")
                fc_mcmcs[j].manual_jump (thetak.get_copy (), thetak_l)
                fc_mcmcs[k].manual_jump (thetaj.get_copy (), thetaj_l)
        
        sample = []
        likls  = []
        for i in range (len (betas)):
            sub_sample, likeli = fc_mcmcs[i].get_last_sampled (1)
            sample.append (sub_sample[0])
            likls.append (likeli[0])
        return (betas, sample, likls)


    def get_last_sampled (self, N):
        betas = self.__betas
        fc_mcmcs = self.__fc_mcmcs
        sample = []
        likelihoods = []
        for i in range (len (betas)):
            sub_sample, likeli = fc_mcmcs[i].get_last_sampled (N)
            sample.append (sub_sample)
            likelihoods.append (likeli)
        return (betas, sample, likelihoods)
