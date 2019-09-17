import SigNetMS

score = SigNetMS.perform_marginal_likelihood (
        'input/simple_enzymatic/simple_enzymatic.xml',
        'input/simple_enzymatic/simple_enzymatic.priors',
        'input/simple_enzymatic/simple_enzymatic.data',
        200, 200, 30, 30, n_process=4, verbose=False, 
        sample_output_file='simple_enzymatic_sample.txt')
