# Parameter inference for Zebrafish ABM
Code for conducting parameter inference for an agent-based model of zebrafish patterns, associated with the paper:

* Liu & Volkening (2025). *Multi-objective Bayesian inference in an agent-based model of zebrafish patterns via topological data analysis*. arXiv 2605.18685. [doi:10.48550/arXiv.2605.18685](https://doi.org/10.48550/arXiv.2605.18685)

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

These Python packages need to be installed (the versions used by the authors are in parenthesis):

* numpy (2.2.4)
* persim (0.3.8) 
* ripser (0.6.12)

## Key files (look at these first)

* tda/params_to_tda.m: Master function for running the TDA pipeline from model parameters, to persistent homology computation, then compute persistence landscape/persistence images
* inference/melanophore_params_inference: demonstrates how to run the inference pipeline to infer melanophore birth parameters, reproduces Section 3.1 of the paper. This script is meant to illustrate and clarify the pipeline; it can't actually run without the zebrafish model, the code for which is available separately
* example_results/ground_truth/groundtruth_wt.mat: simulation result of the wildtype fish with ground truth parameters. This file is provided so you can try out the TDA pipeline without having to simulate the model yourself

## Utility files

* tda/compute_pimg.m: Compute persistence images from persistent homology barcodes
* tda/computePersistentHomology.py: Compute persistent homology with the Vietoris-Rips filtration. This script is meant to be called from within Matlab and can't be run directly from python
* tda/formatModelData.m: Compute pairwise distance between cells, and format simulation data to be ready for TDA computation
* tda/plot_barcode.m, plot_pland.m: plotting functions for persistent homology barcodes and persistence landscapes, respectively
* tda/tda_default_params.m: default hyper-parameter values for TDA computation


## Example

To compute TDA persistence homology, persistence landscape and persistence image for the ground truth wild-type fish: in MATLAB, run 
```
sim_params = zebrafish_default_params;
randomseed = 10001;
T = 45;
tda_params = tda_default_params;
ts=1:T;
savefolder = '/path/of/repo/zebrafish_inference_public/example_results/ground_truth/';
simname = 'groundtruth_wt';
params_to_tda(sim_params,randomseed,T,tda_params,ts,savefolder,simname,0,1);
```

## Troubleshooting regarding MATLAB-Python interface

Some scripts involves calling Python functions from within MATLAB. There are a few possible problems that can arise from this, and some possible solutions. This functionality is very finicky to use, so the actual solution might depend on the exact set up of your computer.

* Potential problem: MATLAB and Python's version are not compatible

See https://www.mathworks.com/support/requirements/python-compatibility.html . Make sure to have a Python version compatible with your MATLAB version installed.

* Potential problem: MATLAB and Python was compiled for different CPU architecture

If you are running the code on an Apple computer, it is possible that the version of MATLAB you are using is the "maci64" version. This version is compiled for Intel architecture (i.e. x86_64), while the Python interpreter is compiled for Apple's new architecture (M1, M2, etc, which is arm64). You can find which architecture your computer has with ```uname -a``` in bash. This can result in cryptic error messages. Solution is to install the "maca64" version of MATLAB, which was compiled for Apple natively. See https://www.mathworks.com/matlabcentral/answers/1977529-how-to-use-python-from-matlab-on-mac-with-apple-silicon

* Potential problem: Python can't find installed packages

When MATLAB calls Python, sometimes sys.path is not initiated correctly. Make sure the packages are actually installed, use ```pip show $PACKAGE_NAME``` to find out where it is, then manually add the path to Python's search path:
```
import sys
sys.path.append('/path/to/package/')
```

* Potential problem: Something about "GLIBCXX_3.4.29" or something similar

This is caused by MATLAB and the version of Python it calls are not using the same version of libstdc++ (GNU C++ library). Try to specity which version of Python you are calling from MATLAB by modifying the "pyenv" command, e.g. 
```
pyenv(Version="/path/to/python")
```

* Potential problem: Something about libstdc++.so.6

This is again about libstdc++. The error message is about MATLAB came with a version of libstdc++ not compatible with Python's version. Try something along the lines of
```
conda install -c conda-forge libstdcxx-ng=12
```
or the pip equivalent if you are not using Conda. Then locate the copy of libstdc++.so.6 you just installed, then do
```
export LD_PRELOAD=/path/to/libstdc++.so.6
```

## Contact

If you run into problems when trying these code, contact Yue Liu at liu4194 at purdue dot edu.