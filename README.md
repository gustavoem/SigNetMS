# SigNetMS
[![Build Status](https://travis-ci.org/gustavoem/SigNetMS.svg?branch=master)](https://travis-ci.org/gustavoem/SigNetMS) 
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)


## SigNetMS stands for Signaling Network Model Selection
This is a Python program that is able to rank biochemical models given a set of experiments. These models should represent changes of concentrations of species of a signaling pathway through a set of ordinary differential equations.
The score is given is an estimate of the probability of the observed data being generated by the model, <img src="https://latex.codecogs.com/gif.latex?p(Data|Model)" />. This estimate is given in a Bayesian approach, and it is based on the work of Tian-Rui Xu et. al in "Inferring Signaling Pathway Topologies from Multiple Perturbation Measurements of Specific Biochemical Species". 

The main result used to estimate this probability is the follwing thermodynamic integral: 

<img src="https://latex.codecogs.com/gif.latex?\log(D|M)=\int_{0}^{1}E_{\theta|D,t}\[p(D|M,\theta)\]dt"/>

where <img src="https://latex.codecogs.com/gif.latex?D"/> is the data, <img src="https://latex.codecogs.com/gif.latex?M"/> is the model, <img src="https://latex.codecogs.com/gif.latex?\theta"/> are model parameters, and <img src="https://latex.codecogs.com/gif.latex?t"/> is a parameter that defines a power-posterior distribution (an intermediate distribution between prior and posterior distribution). To estimate this integral, we need samples from <img src="https://latex.codecogs.com/gif.latex?\theta|D,t"/> (power-posterior t of theta) and to do so, we use different variations of the Metropolis-Hastings algorithm.


## Sampling algorithms
To estimate the value of the power-posterior distribution we should sample from multiple power-posteriors distributions. To do so, we use the Metropolis-Hastings algorithm with multivariate truncated normal jump distributions. We use three different steps of sampling. On the first step, the proposal distribution used has a diagonal covariance matrix, i.e. the jumps of parameters are proposed independently. On the second and third steps, the covariance of the jump distribution is set as an estimate of the covariance of the current sample of the power-posterior. On the first two steps, the sampling is independent for each power-posterior, however, on the last step we
use Populational Monte Carlo Markov Chain, which allow us to mix samples of different power-posteriors.


## Running SigNetMS
To run SigNetMS, you will need to provide to the program some arguments that are related to the problem instance and some that are related to the sampling algorithms used to calculate the desired estimate. 

`SigNetMS.py [-h] [--verbose [VERBOSE]] model priors experiment first_sampling_iterations sigma_update_n second_sampling_iterations third_sampling_iterations`

The arguments related to the problem instance are:
* `model` - an SBML file with defined kinetic laws;
* `priors` - an XML file that with the prior distribution of the model parameters;
* `experiment` - an XML file with the experiments observations.
Some examples of these files are provided in the `input` folder.

The arguments related to the sampling algorithms are:
* `first_sampling_iterations ` - number of iterations on the first sampling step;
* `sigma_update_n` - number of iterations between updates of covariance matrix on the first step;
* `second_sampling_iterations` - number of iterations on the second sampling step;
* `third_sampling_iterations` - number of iterations on the third sampling step;

The program also has optinal arguments.
* `--verbose` if you'd like a verbose run.
* `--help` if you need help.

## Sample inputs
A sample command to run SigNetMS on `bioinformatics` dataset:

`python SigNetMS.py input/bioinformatics/model[1-4].xml input/bioinformatics/model.priors input/bioinformatics/experiment.data 10000 1000 5000 5000`


A sample command to run SigNetMS on `smallest` dataset:

`python SigNetMS.py input/smallest/model[1-4].xml input/smallest/model[1-4]_gamma.priors input/smallest/experiment.data 10000 1000 5000 5000`
