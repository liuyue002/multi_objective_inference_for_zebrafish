# Parameter inference for Zebrafish ABM
Code for conducting parameter inference for an agent-based model of zebrafish patterns, associated with the paper:

* Liu & Volkening (2025). *Multi-objective Bayesian inference in an agent-based model of zebrafish patterns via topological data analysis*. In preparation.

The code in this repository uses the Zebrafish model from the following paper. The code for simulating the model is available separately upon request.

* Volkening & Sandstede (2018). *Iridophores as a source of robustness in zebrafish stripes and variability in Danio patterns*. Nature Communications 9(1):3231

The TDA (topological data analysis) code are partly adapted from 

* Cleveland, Zhu, Sandstede, & Volkening (2023). *Quantifying Different Modeling Frameworks Using Topological Data Analysis: A Case Study with Zebrafish Patterns*. SIAM Journal on Applied Dynamical Systems 22(4):3233-3266

## What does this do?
Models for biological systems tend to contain parameters, whose values are unknown *a priori*, and need to be determined by caliberating the model to data. For ODE or PDE-based models, techniques for model caliberation are well studied. But for agent-based models, model caliberation is much trickier, due to both the stochasticity, intractable likelihood function, and the cost of simulating the model. For phenomenon such as zebrafish skin patterning, an additional challenge is to extract relevent features from the pattern that encodes relevant geometric information.

The approach employed here use TDA to extract features from zebrafish patterns, which is used to establish a metric for quantifying similarity between patterns. This metric is then used in conjunction with AABC (approximate-approximate Bayesian computation) to construct approximate posterior distribution for the model parameters. These posteriors provides an estimate for the parameters, and also their identiifiability.

For details, see the Liu & Volkening paper.

## Language
Most of the code is in MATLAB, with a few components in Python.

## Key files

TODO: list the important files

## Utility files

TODO: list things like plotting functions, etc

## Troubleshooting regarding MATLAB-Python interface

TODO: talk about how to call python script from MATLAB

* Potential problem: MATLAB and Python's version are not compatible
* Potential problem: MATLAB and Python was compiled for different CPU architecture
* Potential problem: Python can't find installed packages

## Contact

If you run into problems when trying these code, contact Yue Liu at liu4194 at purdue dot edu.