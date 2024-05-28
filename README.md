# Solving Epidemic Models with Explicit Method
This repo contains a modeling implementation of SIR epidemic model with the R-Shiny interface that can be directly compiled with R-Shiny interface in the GAMS modeling language. 

SIR is a model that models the dynamic interactions between three groups of people, susceptible individuals (S), infected individuals (I) and removed individual (R), using orinary differential equation. With given starting and ending ratio of SIR proportion, we are able to find the rate of recovery by solving an optimization model by using Euler method (a discretization and numerical approximation method).

We also have a model for SIRV in this repo, which is a variation of SIR model. It is implemented with second-order Adam Bashforth method and itt aims to find the optimal control of the vaccination rate at each discretization time stamp.

## Folders and files
The application based on SIR model contains :
	1. SIR.gms : SIR model formulated in GAMS with Euler Method
	2. conf_SIR : configure file for applications
	3. data_SIR : data file for the default 
	4. renderer_sir : R-Shiny renderer files

The experiment of AB2 on SIRV model in the project is included in:
	1. SIRV_AB2 : contains SIRV epidemic model with second-order Adam Bashforth method

## More comments
More details are included in the `technical_report.pdf`
