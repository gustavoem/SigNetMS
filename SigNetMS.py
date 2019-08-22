from model.SBML import SBML
from model.SBMLtoODES import sbml_to_odes
from marginal_likelihood.MarginalLikelihood import MarginalLikelihood
from model.PriorsReader import define_sbml_params_priors
from experiment.ExperimentSet import ExperimentSet
import argparse

def perform_marginal_likelihood (sbml_file, priors_file, \
        experiment_file, burnin1_iterations, sigma_update_n, \
        burnin2_iterations, sampling_iterations, verbose=False, \
        n_process=0):
    return 0
    print  ("Performing marginal likelihood calculations of model: " + \
            sbml_file)
    sbml = SBML ()
    sbml.load_file (sbml_file)
    odes = sbml_to_odes (sbml)
    experiments = ExperimentSet (experiment_file)
    theta_priors = define_sbml_params_priors (sbml, priors_file)

    ml = MarginalLikelihood (burnin1_iterations, 
            sigma_update_n, 
            burnin2_iterations, 
            sampling_iterations, 20, 2, \
            verbose=verbose, n_process=n_process)
    log_l = ml.estimate_marginal_likelihood (experiments, odes, 
            theta_priors)
    print ("log_l = " + str (log_l))
    return log_l



if __name__ == "__main__":
    parser = argparse.ArgumentParser ()
    parser.add_argument ("model", help="SBLM file with model definition.")
    parser.add_argument ("priors", help="An XML file with the priors for" \
            + " the model parameters.")
    parser.add_argument ("experiment", help="An XML file with the" \
            + " experiments observations.")
    parser.add_argument ("first_sampling_iterations", help="How many" \
            + " iterations should be performed on the first step of the" \
            + " parameter sampling.")
    parser.add_argument ("sigma_update_n", help="Iterations between" \
            + " each time the jumping variance is updated according to" \
            + " acceptance rate.")
    parser.add_argument ("second_sampling_iterations", help="How many" \
            + " iterations should be performed on the second step of the" \
            + " parameter sampling.")
    parser.add_argument ("third_sampling_iterations", help="How many" \
            + " iterations should be performed on the third step of the" \
            + " parameter sampling.")
    parser.add_argument ('--verbose', type=bool, nargs='?', const=True, \
            help="Verbose run.")
    parser.add_argument ('--n_process', type=int, nargs='?', default=0, \
            help="Number of parallel process used on sampling step.")
    args = parser.parse_args ()
    

    # Problem input
    sbml_file = args.model
    priors_file = args.priors
    experiment_file = args.experiment

    # Algorithm parameters
    first_step_n = int (args.first_sampling_iterations)
    sigma_update_n = int (args.sigma_update_n)
    second_step_n = int (args.second_sampling_iterations)
    third_step_n = int (args.third_sampling_iterations)
    verbose = args.verbose
    n_process = args.n_process

    perform_marginal_likelihood (sbml_file, priors_file, \
            experiment_file, first_step_n, sigma_update_n, \
            second_step_n, third_step_n, verbose=verbose, \
            n_process=n_process)
