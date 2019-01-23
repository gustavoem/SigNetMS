# SigNetMS
[![Build Status](https://travis-ci.org/gustavoem/SigNetMS.svg?branch=master)](https://travis-ci.org/gustavoem/SigNetMS) 
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)


## SigNetMS stands for Signaling Network Model Selection
This is a Python program that is able to give a score to a model of ordinary differential equations given a set of experiments.
The score is given is an estimate of the probability of the observed data being generated by the model, <img src="https://latex.codecogs.com/gif.latex?p(Data|Model)" />. This estimate is given in a Bayesian approach, and it is based on the work of Tian-Rui Xu et. al in "Inferring Signaling Pathway Topologies from Multiple Perturbation Measurements of Specific Biochemical Species". 

The main result used to estimate this probability is the follwing thermodynamic integral: 

<img src="https://latex.codecogs.com/gif.latex?\log(D|M)=\int_{0}^{1}E_{\theta|D,t}\[p(D|M,\theta)\]dt"/>

where <img src="https://latex.codecogs.com/gif.latex?D"/> is the data, <img src="https://latex.codecogs.com/gif.latex?M"/> is the model, <img src="https://latex.codecogs.com/gif.latex?\theta"/> are model parameters, and <img src="https://latex.codecogs.com/gif.latex?t"/> is the "temperature". To estimate this integral, we need samples from <img src="https://latex.codecogs.com/gif.latex?\theta|D,t"/> and to do so, we use different variations of the Metropolis-Hastings algorithm.


## Sampling algorithms
*TODO*

## Running SigNetMS
To run SigNetMS, you will need to provide to the program some arguments that are related to the problem instance and some that are related to the sampling algorithms used to calculate the desired estimate. 

The arguments related to the problem instance are:
* `model` - an SBML file with defined kinetic laws;
* `priors` - an XML file that with the prior distribution of the model parameters;
* `experiment` - an XML file with the experiments observations.
Some examples of these files are provided in the `input` folder.

The arguments related to the sampling algorithms are:
* `first_sampling_iterations `
* `sigma_update_n`
* `second_sampling_iterations`
* `third_sampling_iterations`

The program also has optinal arguments
* `--verbose` if you'd like a verbose run
* `--help` if you need help
