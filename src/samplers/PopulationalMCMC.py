import numpy as np
from LikelihoodFunction import LikelihoodFunction
from DiscreteLaplacian import DiscreteLaplacian
from utils import safe_power


class PopulationalMCMC:

    """ Performs a tempered sampling of theta. The sampling chooses 
        temperatures and for each of them conduces a 
        FixedCovarianceMCMC. The algorithm also mixes samples of 
        different temperatures. """


    __SCHEDULE_POWER = 4


    def __init__ (self, n_strata, strata_size, fc_mcmcs):
        """ Default constructor. """
        if n_strata * strata_size != len (fc_mcmcs):
            raise ValueError ("The list of covariances and starts " \
                    + "should have the same size as n_strata * " \
                    + "strata_size")

        self.__n_strata = n_strata
        self.__strata_size = strata_size
        self.__sample_betas (n_strata, strata_size)
        self.__fc_mcmcs = fc_mcmcs
        self.__define_samplers_temp (self.__betas, self.__fc_mcmcs)

    
    def __sample_betas (self, n_strata, strata_size):
        """ Samples the temperatures array. """
        betas = []
        sched_power = PopulationalMCMC.__SCHEDULE_POWER
        for i in range (n_strata):
            strata_start = (i /  n_strata) ** sched_power
            strata_end = ((i + 1) / n_strata) ** sched_power
            for j in range (strata_size):
                x = np.random.uniform (strata_start, strata_end)
                betas.append (x)
            betas.sort ()
        self.__betas = betas


    def __define_samplers_temp (self, betas, fc_mcmcs):
        for i in range (len (betas)):
            fc_mcmcs[i].set_temperature (betas[i])


    def get_sample (self, N):
        """ Get a sample of size N. """
        betas = self.__betas
        fc_mcmcs = self.__fc_mcmcs
        for i in range (N):
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
            tjotk = np.exp (thetaj_l - thetak_l)
            tkotj = np.exp (thetak_l - thetaj_l)
            j_gv_k = inv_temp_jump_dist.pdf (j + 1)
            k_gv_j = temp_jump_dist.pdf (k + 1)
            
            r = safe_power (tkotj, betas[j]) * \
                safe_power (tjotk, betas[k]) * \
                (j_gv_k / k_gv_j)

            if np.random.uniform () <= r:
                fc_mcmcs[j].manual_jump (thetak.get_copy (), thetak_l)
                fc_mcmcs[k].manual_jump (thetaj.get_copy (), thetaj_l)
        
        sample = []
        likls  = []
        for i in range (len (betas)):
            sub_sample, likeli = fc_mcmcs[i].get_last_sampled (1)
            sample.append (sub_sample[0])
            likls.append (likeli[0])
        return (sample, likls)
