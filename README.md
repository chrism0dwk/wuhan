# Early analysis of Wuhan nCoV-2019 outbreak in Wuhan, China.

This R package implements an ODE-based model of the novel coronavirus
outbreak in Wuhan, China.  It presents a simulator and likelihood function
assuming Poisson-distributed increments in the number of new cases in Wuhan,
in the rest of China via the airline network, and to the rest of the world.

__Data required__: 

* `china_cases` daily case reports in all Chinese cities (see `data(package='wuhan')`)
* `world_cases` daily case reports from other countries (see `data(package='wuhan')`)
* `K` daily numbers of passengers going between cities in China via airline network, available from OAG Traffic Analyzer
* `W` daily numbers of passengers going between Chinese cities and other countries via airline network, available from OAG Traffic Analyzer
* `china_population` the population size in each Chinese city (see `data(package='wuhan')`)
    
__Parameters__:

* `beta` the human-human basic transmission rate
* `gamma` the infectious period
* `I0W` the number of initial infectives in Wuhan
* `phi` the case ascertainment rate in Wuhan
    
To use the package, assume the following workflow in R:

````r
# Load required packages
> install.packages('devtools')
> devtools::install_git('https://github.com/chrism0dwk/wuhan.git')
> library(wuhan)

# Instantiate ODE model, simulate up to day 22.
> simulator = NetworkODEModel(N=china_population, K=K, init_loc='Wuhan', alpha=1/4, max_t=22) 

# Instantiate LogLikelihood function
> llik = LogLikelihood(y=china_cases, z=world_cases, N=N, K=K, W=W, sim_fun=simulator)

# Find MLEs using optimisation
> p_hat = optim(log(par_init), llik, control=list(fnscale=-1))
````

Asymptotic assumptions for confidence intervals fail in our case, since the
parameter space is highly non-orthogonal.  Confidence intervals are therefore
calculated using parametric bootstrap.  `p_hat` is calculated on the log scale (logit scale
for the `phi` parameter), so needs to be transformed first:

````r
> p_hat[1:3] = exp(p_hat[1:3])
> p_hat[4] = exp(p_hat[4])
````

The samples can then be drawn by bootstrap, for which a computing cluster is
highly recommended (thanks Lancaster University HEC facility!).
````r
> samples = bootstrap(p_hat, K, W, alpha=1/4, max_t=22, n_samples=1000)
````

Since the airline connectivity matrices are not included in this package, samples 
from the parameters (for 4 different values of the latent period $1/\alpha$) are 
provided as in-build datasets.  See `data(package='wuhan')`.
