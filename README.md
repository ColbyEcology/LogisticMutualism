# LogisticMutualism

This repository contains code for running simulations of communities of species over time, using random interaction matrices. The purpose of these simulations is to explore the differences between two commonly used mathematical models of mutualism, the "intraspecific-benefit" model and the "equilibrium benefit" model. 

## Files
simulation_functions.R : contains the functions needed to run the simulation. 

network_simulation.Rmd : contains the code for running the simulation. This includes setting ranges of parameters, looping through many iterations of the simulation using the parameter ranges, and compiling the data produced by the simulations. The output of this code is a dataframe with data from each simulation iteration, including the model and parameters used, the number of species that persist and go extinct, and whether or not exponential growth was produced in the simulation.

simData_analysis.Rmd: contains analysis of the data produced by the simulation.

## simulation_functions contents
generate_coef_mat(n_spp, distribution, disp): generates a random coefficient matrix for a community with the given number of species. Distribution must be "rnorm" for a normal distribution of coefficients or "runif" for a uniform distribution. Coefficients range from -disp to disp.
        
generate_coef_mat2(n_spp, distribution, disp): generates a random coefficient matrix with at least one mutualism
   
get_indices(mat): takes a coefficient matrix and returns a list of the mutualism indices

get_mut_mat(mat): takes a coefficient matrix and returns a new matrix with only mutualisms

get_not_mut_mat(mat): takes a coefficient matrix and returns a new matrix with only non-mutualistic relationships

get_num_persist(out): given a deSolve ODE solver output, returns the number of species which persisted (ended with density > 0)

errorFunc(e): tryCatch error handler for exponential growth

rootfun(t, n, parms): root function to trigger event in the case of extinction or exponential growth

eventfun(t, n, parms): event function to set density to 0 in the case of extinction or throw an error in the case of exponential growth

odeSolve(y, times, func, parms, method, rootfun, events): given usual ode parameters, runs the ODE solver with tryCatch and returns a deSolve object if there is no exponential growth and an integer (1) if there is.

expGrowth(out): given the output from odeSolve(), returns T if there was exponential growth

sumExtPer(out, expGrowth, n_spp): given outputs from odeSolve() and expGrowth() and the number of species, returns the number of extinctions and persistances.

NumeratorMut(t, n, parms): deSolve function for the numerator (intraspecific benefit) mutualism model

DenominatorMut(t, n, parms): deSolve function for the denominator (equilibrium benefit) mutualism model
    
    
    
