# In this module we calculate the marginal likelihood of models given
# observed data using the methodology presented in "Inferring Signaling 
# Pathway Topologies from Multiple Perturbation Measurements of Specific 
# Biochemical Species", Tian-Rui Xu, et. al.
from marginal_likelihood.samplers.AcceptingRateAMCMC import \
        AcceptingRateAMCMC
from marginal_likelihood.samplers.AdaptingCovarianceMCMC import \
        AdaptingCovarianceMCMC
from marginal_likelihood.samplers.FixedCovarianceMCMC import \
        FixedCovarianceMCMC
from marginal_likelihood.samplers.PopulationalMCMC import \
        PopulationalMCMC
import multiprocessing

from parallel_map import parallel_map
import numpy as np
class MarginalLikelihood:
    """ This class is able to perform an adaptive MCMC sampling to 
        estimate the likelihood of a model given experimental data. """
    
    def __init__ (self, phase1_iterations, sigma_update_n, 
            phase2_iterations, phase3_iterations, n_strata, 
            strata_size, verbose=False, n_process=0):
        """ Default constructor. phase1_iterations is the number of 
            iterations performed by the AcceptingRateAMCMC, which is
            an adaptive sampler that performs independent MCMC on each
            system variable and adapts the jump covariance with the goal
            of having acceptance rate between 0.25 and 0.75. 
            sigma_update_n is the number of iterations before updating 
            sigma in intial phase. phase2_iterations is the number of 
            iterations performed in the second phase of sampling with 
            an adaptive sampling algorithm that has jumps distributed 
            with a probability function that has as covariance matrix
            a matrix of covariance estimated from the current sample.
            phase3_iterations is the number of iterations performed by
            the PopulationalMCMC algorithm. n_strata is the number of 
            strata used in the populational algorithm, and strata_size 
            is the number of individuals per strata."""
        self.__phase1_iterations = phase1_iterations
        self.__phase2_iterations = phase2_iterations
        self.__phase3_iterations = phase3_iterations
        self.__sigma_update_n = sigma_update_n
        self.__n_strata = n_strata
        self.__strata_size = strata_size
        self.__verbose = verbose
        self.__sample = None

        if n_process == 0:
            self.__n_process = max (1, \
                    multiprocessing.cpu_count () // 2)
            if verbose:
                print ("Using automatic number of process.")
        else:
            self.__n_process = n_process

    
    @staticmethod
    def __run_phase_one_and_two (temp, experiments, model, theta_prior,
            n_acc, n_adap_cov, n_sigma_update, verbose):
        """ Map function to run phase 2 and 3 for each temperature. """
        np.random.RandomState ()
        # Phase 1
        acc_mcmc = AcceptingRateAMCMC (theta_prior, model, experiments, 
                n_sigma_update, verbose=verbose)
        acc_mcmc.set_temperature (temp)
        acc_mcmc.start_sample_from_prior ()
        sample, likelis = acc_mcmc.get_sample (n_acc)

        # Phase 2
        adap_cov_mcmc = AdaptingCovarianceMCMC (theta_prior, model, 
                experiments, n_sigma_update, verbose=verbose)
        adap_cov_mcmc.set_temperature (temp)
        adap_cov_mcmc.define_start_sample (sample, likelis)
        sample, likelis = adap_cov_mcmc.get_sample (n_adap_cov)

        # Construct phase 3 local temperature sampler
        S = adap_cov_mcmc.get_jump_covariance ()
        fc_mcmc = FixedCovarianceMCMC (theta_prior, model, 
                experiments, S, t=temp, verbose=verbose)
        theta = sample[-1]
        log_likeli = likelis[-1]
        fc_mcmc.define_start_sample ([theta], [log_likeli])
        return fc_mcmc


    def __set_sample (self, betas, thetas, log_ls):
        """ Defines the sample created in estimate_marginal_likelihood.
        
        Parameters
            beta: a list of temperatures for which the sample was 
                created.
            thetas: the list of parameters sampled, 
                for each temperature. This should be represented as a 
                list of list; the index of the first list is associated
                to a temperature of same index in beta, and the second
                list represent the sample for that temperature.
            log_ls: has the same indexing as thetas, but insted of 
                storing parameters, stores the log-likelihood of that
                parameter.

        Notes
            This method sets the private attribute __sample, which 
            stores the created sample as a list tuples. Each tuple
            is in the format of (temperature, sample for temperature).
            The sample part of the tuple is also a list of tuples, 
            containing each tuple the format (parameter values, 
            log-likelihood).
        """
        self.__sample = []
        for i in range (len (betas)):
            temp_sample = []
            temp = betas[i]
            for j in range (len (thetas[i])):
                theta = thetas[i][j]
                log_l = log_ls[i][j]
                temp_sample.append ((theta, log_l))
            self.__sample.append ((temp, temp_sample))


    def print_sample (self, output_file=None):
        """ Prints the current self.__sample. 
            
            Parameters
                output_file: a string with a filename. If this is 
                    provided, the sample is also printed this file.
            Note
                This method should only be called after a succesful call
                of estimate_marginal_likelihood.
        """
        if self.__sample == None:
            raise ValueError ("Sample is undefined, before printing" \
                    + " sample, you should first call the method" \
                    + " estimate_marginal_likelihood.")

        if output_file:
            file_obj = open (output_file, "w")

        for temp_tuple in self.__sample:
            temp = temp_tuple[0]
            temp_sample = temp_tuple[1]
            print_str = "Sample for t = " + str (temp)
            print (print_str)
            if output_file:
                file_obj.write (print_str + "\n")
            for sample_tuple in temp_sample:
                theta = sample_tuple[0]
                log_l = sample_tuple[1]
                print_str = "parameter = " + str (theta.get_values ()) \
                        + " likelihood = " + str (log_l)
                print (print_str)
                if output_file:
                    file_obj.write (print_str + "\n")

        
    def estimate_marginal_likelihood (self, experiments, model, 
            theta_prior):
        """ This function estimates the marginal likelihood of a  model.
        """
        n_pop = self.__phase3_iterations
        n_strata = self.__n_strata
        strata_size = self.__strata_size
        verbose = self.__verbose
        n_process = self.__n_process
        self.__sample = None

        # initialize ODEs function and jacobian
        model.evaluate_on ([experiments[0].times[0]])
        model.get_system_jacobian ()
        print ("Initialized ODEs")

        betas = PopulationalMCMC.sample_scheduled_betas (n_strata * 
                strata_size)
        print ("Phase 1 and 2 starts.")
        phase_1_n_2_f = lambda temp : \
                MarginalLikelihood.__run_phase_one_and_two (temp, \
                experiments, model, theta_prior, 
                self.__phase1_iterations, self.__phase2_iterations,
                self.__sigma_update_n, self.__verbose) 
        fc_mcmcs = parallel_map (phase_1_n_2_f, betas, n_process)
                       
        print ("Phase 3 starts.")
        # Phase 3
        pop_mcmc = PopulationalMCMC (n_strata, strata_size, fc_mcmcs,
                betas=betas, verbose=verbose)
        pop_mcmc.get_sample (n_pop)
        betas, thetas, log_ls = pop_mcmc.get_last_sampled (n_pop // 4)
        self.__set_sample (betas, thetas, log_ls)
        print ("Sampling ended.")
        
        if self.__verbose:
            print ("Here are the sampled parameters separated by" \
                    " temperature.")
            self.print_sample ()
            
        ml = self.__calculate_marginal_likelihood (betas, thetas, 
                log_ls)
        return ml


    def __create_fcmcmc_samplers (self, adap_cov_mcmc, m, experiments,
            model, theta_prior):
        """ Creates a list of FixedCovarianceMCMC objects that has the
            same jump covariance matrix as AdaptiveCovarianceMatrix. """
        S = adap_cov_mcmc.get_jump_covariance ()
        start_t, start_l = adap_cov_mcmc.get_last_sampled (1)
        fc_mcmcs = []
        for _ in range (m):
            sampler = FixedCovarianceMCMC (theta_prior, model, 
                    experiments, S)
            sampler.define_start_sample (start_t, start_l)
            fc_mcmcs.append (sampler)
        return fc_mcmcs


    def __calculate_marginal_likelihood (self, betas, thetas, \
            likelihoods):
        """ Given a list with samples, calculates the marginal 
            likelihood. """
        ml = 0
        
        if self.__verbose:
            print ("Estimating marginal likelihood")

        exp_gv_tip1 = 0
        power_sample_ip1 = thetas[0]
        for j in range (len (power_sample_ip1)):
            log_l = likelihoods[0][j]
            exp_gv_tip1 += log_l
        exp_gv_tip1 /= len (power_sample_ip1)

        for i in range (len (betas) - 1):
            tip1 = betas[i + 1]
            ti = betas[i]
                        
            exp_gv_ti = exp_gv_tip1
            exp_gv_tip1 = 0
            power_sample_ip1 = thetas[i + 1]
            for j in range (len (power_sample_ip1)):
                log_l = likelihoods[i + 1][j]
                exp_gv_tip1 += log_l
            exp_gv_tip1 /= len (power_sample_ip1)

            ml += (tip1 - ti) * (exp_gv_tip1 + exp_gv_ti)
            
            if self.__verbose:
                print ("first temperature = " + str (ti))
                print ("second temperature = " + str (tip1))
                print ("log ml given first temperature = " + \
                        str (exp_gv_ti))
                print ("log ml given second temperature = " + \
                        str (exp_gv_tip1))

        ml /= 2
        print ("Calculated log marginal likelihood: " + str (ml))
        return ml
