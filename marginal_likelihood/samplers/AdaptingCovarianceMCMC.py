import numpy as np
from marginal_likelihood.samplers.MetropolisHastings import \
        MetropolisHastings
from marginal_likelihood.LikelihoodFunction import LikelihoodFunction
from distributions.MultivariateLognormal import MultivariateLognormal
from covariance_estimate import calc_covariance
from utils import safe_log
from utils import safe_exp_ratio
from utils import safe_pow_exp_ratio
from utils import get_current_datetime
from pathlib import Path 

class AdaptingCovarianceMCMC (MetropolisHastings):
    """ Objects of this class are able to return a sample of theta using 
        an adaptive metropolis-hastings algorithm. The proposal 
        distribution of the algorithm is adaptive and has as covariance 
        matrix an estimate of the covariance according to the current 
        sample. To get a sample with an instance of this class, you will 
        need a starting sample to provide a first estimate of the 
        covariance matrix. """

    def __init__ (self, theta, model, experiments, covariance_rescale_n,
            t=1, verbose=False):
        """ Default constructor. 
        
            Parameters
                theta: a RandomParameterList with the parameter priors.
                model: an SBML object with the model.
                experiments: an ExperimentSet object with the 
                    experimental data.
                covariance_rescale_n: the number of iterations before
                    a rescale on the proposal distribution covariance.
                t: a float indicating the tempering parameter of the
                    target distribution.
                verbose: boolean indicating if verbosity is wanted.
            
            Returns
                an AdaptingCovarianceMCMC object.
        """
        super ().__init__ (theta, verbose=verbose)
        self.__model = model
        self.__experiments = experiments
        self._jump_S = None
        self._jump_scale = 1
        self._covariance_rescale_n = covariance_rescale_n
        self.__l_f = LikelihoodFunction (model)
        self._t = t


    def set_temperature (self, t):
        """ Defines the tempering parameter.
        
            Parameters
                t: a float with the tempering parameter.
        """
        self._t = t


    def _open_trace_file (self):
        """ Open a file to write trace. """
        trace_dir = "trace/" + self.__model.name
        Path(trace_dir).mkdir(parents=True, exist_ok=True)

        file_name = trace_dir + "/" + get_current_datetime () + "_" \
                + str (self._t) + "_" + "2nd_phase"
        self._trace_file = open (file_name, 'w')


    def __calc_jump_S (self):
        """ Calculates jump_S, an estimate of the covariance of 
            parameters. """
        sample_l_values = []
        for t in self._sample:
            log_values = [safe_log (x) for x in t.get_values ()]
            sample_l_values.append (log_values)
        self._jump_S = calc_covariance (sample_l_values)
    
    
    def get_jump_covariance (self):
        """ Returns a copy of jump_S matrix. """
        return np.array (self._jump_S)


    def _create_jump_dist (self, theta_t):
        """ Creates the jump distribution from the current point. 
        
            Parameters
                theta_t: a RandomParameterList with the current point.
            
            Returns
                A MultivariateLognormal distribution which is the jump
                distribution from the current point. This distribution
                have the normal parametrization with 
                mu = log(current_point_values) and covariance equal
                to the sample covariance of the log of the accepted
                points. 
        """
        t_vals = theta_t.get_values ()
        mu = np.array (np.log (t_vals))
        dist = MultivariateLognormal (mu, 
                self._jump_S * self._jump_scale)
        return dist


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
        if self._is_verbose:
            self._trace_file.write ("jump_S = ")
            self._trace_file.write (str (self._jump_S))

        j_gv_old = self._create_jump_dist (old_t)
        j_gv_new = self._create_jump_dist (new_t)
        try:
            log_new_gv_old = j_gv_old.log_pdf (new_t.get_values ())
            log_old_gv_new = j_gv_new.log_pdf (old_t.get_values ())
        except Exception as e:
            print (e)
            raise ValueError ("The covariance matrix is not positive " \
                    + "definite. Try using a bigger starting sample.")
        
        l_ratio = self._calc_likeli_ratio (new_l, old_l)
        prior_ratio = safe_exp_ratio (new_t.get_log_p (), 
                old_t.get_log_p ())
        jump_ratio = safe_exp_ratio (log_old_gv_new, log_new_gv_old)
        # if self._is_verbose:
            # print ("\tnew given old: " + str (log_new_gv_old))
            # print ("\told given new: " + str (log_old_gv_new))
            # print ("\tprior ratio: " + str (prior_ratio))
            # print ("\tlikelihood ratio: " + str (l_ratio))
            # print ("\tjump ratio: " + str (jump_ratio))

        if not l_ratio < float ("inf") and prior_ratio + jump_ratio > 0:
            return 1
        elif prior_ratio + jump_ratio < 1e-150:
            return 0
        else:
            mh_ratio = l_ratio * prior_ratio * jump_ratio
            return mh_ratio


    def _calc_likeli_ratio (self, log_new_l, log_old_l):
        """ Calcultes the ratio (new_l / old_l) ^ t. """
        return safe_pow_exp_ratio (log_new_l, log_old_l, self._t)


    def _calc_log_likelihood (self, theta):
        """ Calculates the log of p (experiments | theta, model). """
        return self.__l_f.get_log_likelihood (self.__experiments, theta)


    def _iteration_update (self):
        """ At the end of each sampling iteration, we should update the
            Covariance Matrix. """
        acceptance_rate = self.get_acceptance_ratio ()
        self.__calc_jump_S ()
        if acceptance_rate > .4 and self._jump_scale < 2:
            self._jump_scale *= 1.1
        if acceptance_rate < .25 and self._jump_scale > 1e-2:
            self._jump_scale *= .9
