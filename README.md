# Private Density Estimation

This package contains the algorithm and experiments for the following paper: 

> [Differentially Private Learning of Structured Discrete Distributions](https://papers.nips.cc/paper/5713-differentially-private-learning-of-structured-discrete-distributions)  
> Ilias Diakonikolas, Moritz Hardt, Ludwig Schmidt  
> NIPS 2015


The code is written in Julia (v0.5). The script `src/histogram_experiment.jl` runs the experiments in the paper for synthetic distributions (mixtures of Gaussian, Beta, and Gamma distributions).

For an example of how to use the code, see the Jupyter notebook `src/histogram_approximation.ipynb`.
