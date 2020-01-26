# Early analysis of Wuhan nCoV-2019 outbreak in Wuhan, China.

This R package implements an ODE-based model of the novel coronavirus
outbreak in Wuhan, China.  It presents a simulator and likelihood function
assuming Poisson-distributed increments in the number of new cases in Wuhan,
in the rest of China via the airline network, and to the rest of the world.

__Data required__: 

* `y` daily case reports in all Chinese cities
* `z` daily case reports from other countries
* `K` daily numbers of passengers going between cities in China via airline network
* `W` daily numbers of passengers going between Chinese cities and other countries via airline network
* `N` the population size in each Chinese city
    
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
> simulator = NetworkODEModel(N=N, K=K, init_loc='Wuhan', alpha=1/4, max_t=22) 

# Instantiate LogLikelihood function
> llik = LogLikelihood(y=y, z=z, N=N, K=K, W=W, sim_fun=simulator)

# Find MLEs using optimisation
> p_hat = optim(log(par_init), llik, control=list(fnscale=-1))
````
