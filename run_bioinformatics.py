import SigNetMS
import sys

model = sys.argv[1]
score = SigNetMS.perform_marginal_likelihood (
        'input/bioinformatics/model' + model  + '.xml',
        'input/bioinformatics/model.priors',
        'input/bioinformatics/experiment.data',
        20000, 2000, 3000, 3000, n_process=10, verbose=False, 
        sample_output_file='bioinformatics_sample' + model + '.txt')
