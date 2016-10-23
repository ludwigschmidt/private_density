# Private Density Estimation

This package contains the algorithm and experiments for the following paper: 

> [Differentially Private Learning of Structured Discrete Distributions](https://papers.nips.cc/paper/5713-differentially-private-learning-of-structured-discrete-distributions)  
> Ilias Diakonikolas, Moritz Hardt, Ludwig Schmidt  
> NIPS 2015


The code is written in Julia (v0.5). The script [`src/histogram_experiment.jl`](https://github.com/ludwigschmidt/private_density/blob/master/src/histogram_experiment.jl) runs the experiments in the paper for synthetic distributions (mixtures of Gaussian, Beta, and Gamma distributions).

For an example of how to use the code, see the Jupyter notebook [`examples/histogram_approximation.ipynb`](https://github.com/ludwigschmidt/private_density/blob/master/examples/histogram_approximation.ipynb).

![](https://cdn.rawgit.com/ludwigschmidt/private_density/4076613c8cfd6699b2f5ba0924d3cc9f4f461c00/output/gmm.svg) ![](https://cdn.rawgit.com/ludwigschmidt/private_density/4076613c8cfd6699b2f5ba0924d3cc9f4f461c00/output/non_private.svg) ![](https://cdn.rawgit.com/ludwigschmidt/private_density/4076613c8cfd6699b2f5ba0924d3cc9f4f461c00/output/private.svg)

The three plots above are generated in the notebook [`examples/website_plots.ipynb`](https://github.com/ludwigschmidt/private_density/blob/master/examples/website_plots.ipynb).
