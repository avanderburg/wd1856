# wd1856
Scripts for analyzing data collected on WD 1856+534

In Vanderburg et al. 2020, we reported the discovery of a giant planet candidate transiting a white dwarf called WD 1856+534. The code in this repository is what we used to determine the size and mass of the planet candidate. The code was not written with distribution in mind, but we release it for transparency anyway. There are two main procedures to do this analysis: fitwd1856.pro and wd1856fluxes.pro. 

fitwd1856.pro sets up and performs a Markov Chain Monte Carlo fit to light curves from the Gran Telescopio Canarias and the Spitzer Space Telescope. It uses the affinemcmcfit.pro package given here: https://github.com/avanderburg/affinemcmcfit. The transit fitting code uses software written by Jason Eastman, and routines written by John Johnson and Jason Eastman are scattered throughout. I have uploaded them to this repository for convenience.  The data files are also uploaded here - if you wish to run the code, you will need to change the paths to point to these files. The code has several "stop" statements throughout which were used for debugging. Use the .continue executive command to resume the program after each stop. 

wd1856fluxes takes the output of the MCMC fits, in particular the value of the dilution parameter d, and converts that to an upper limit on the companion's mass using the Sonora brown dwarf model grid. 
