# Early analysis of Wuhan nCoV-2019 outbreak in Wuhan, China.

This R package implements an ODE-based model of the novel coronavirus
outbreak in Wuhan, China.  It presents a simulator and likelihood function
assuming Poisson-distributed increments in the number of new cases in Wuhan,
in the rest of China via the airline network, and to the rest of the world.

__Data required__: 
    * $y$ daily case reports in all Chinese cities
    * $z$ daily case reports from other countries
    * $K$ daily numbers of passengers going between cities in China via airline network
    * $W$ daily numbers of passengers going between Chinese cities and other countries via airline network
    * $N$ the population size in each Chinese city
    
__Parameters__:
    * $\beta$ the human-human basic transmission rate
    * $\gamma$ the infectious period
    * $I0W$ the number of initial infectives in Wuhan
    * $\phi$ the case ascertainment rate in Wuhan
    
See package vignette for instructions on how to use the package.
