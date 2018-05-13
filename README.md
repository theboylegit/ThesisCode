# ThesisCode
Versions of code written for my thesis

ORIGINAL *WORKING* VERSION- simsplash.c

VERSION UNDERGOING UPDATE - simsplashNewBuild.c

This repository contains the original version of the code, SIMSPLASH, written for my thesis work. The code is mostly written in C and contains some calls to python scripts also written to complete my thesis work. All together (including the databases and scripts within the appropriate directories), this code is capable of generating the results contained in my thesis. A packaged version which includes everything needed to run it not yet ready for release due to unpublished data files. 


NOTE: SIMSPLASH is currently being rewritten/updated to be more readable, portable and maintainable. The old version has some bugs, as well as long functions, confusing variable names and perhaps some other no-no's, since it was originally written to calculate the results for my thesis only. The aim is to eventually release it as an open source educational/research tool. The refactored code is here as "simsplashNewBuild.c"


SIMSPLASH has three modes of operation:
1. Orbital evolution mode: Calculates the evolution of a planet orbiting an evolving star, which depends on the stellar model given as input. The outcomes it looks for is whether the planet is engulfed(swallowed by the star it orbits during the star's evolution) on either of the giant branches (Red-GB or Asymptotic-GB) or whether it survives the evolution. The Orbital evolution mode allows the user to investigate which factors in the initial conditions of the star (by varying the input stellar models) and/or planet determine when or if the planet is engulfed. These initial conditions are specifed in the configuration file before runtime.


2. Single epoch population mode: Calculates the orbital evolution for a population of planets orbiting stars and performs some analysis to draw conclusions about the current planet-hosting star population. It reads in a file containing the population of planets and their stars and utilises Orbital Evolution mode described above for each planet-star system and logs the result of whether or not and when each planet is engulfed. Populations are either provided in the Planets_Directory or are downloaded from exoplanets.org database depending on user configuration. After the population is evolved, python scripts are called to perform analysis on the results (logistic regression). The Single Epoch population mode allows the user to investigate the effect of varying and averaging initial conditions of host stars across an entire population. The average initial conditions for the population (stellar models etc) are specified by the user in the configuration file before runtime.


3. Multi epoch population mode: Performs a population synthesis to simulate the current galactic population of stellar remnants. It utilises a star formation history and the results generated with Single Epoch mode to evolve star populations throughout the galactic history, in order to arrive at the present-day population of visible stellar remnants (known as planetary nebulae). The outcome it determines is what proportion of visible planetary nebulae should have evolved such that they engulfed a planet earlier in their lifetime.  
