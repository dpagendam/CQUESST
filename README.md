## Repository Description



### Relevant Files:
	- /CQUESST/CQUESST.stan
	- /CQUESST/data_cmdStan
	- /CQUESST/init_params
	- /CQUESST/read_results.R
	- run_CQUESST.slurm



### Description of individual files:

**/CQUESST/CQUESST.stan**
- This contains the Stan model for implementing CQUESST for the Millennium Tillage Trial (MTT)

**/CQUESST/data_cmdStan**
- This contains the MTT data for ingestion to CQUESST.stan

**/CQUESST/init_params**
- This contains the initial values of the parameters to use for the Markov chains
	
**/CQUESST/read_results.R**
- An R script for reading in the posterior samples (output.csv) produced by running the Stan model
	
**run_CQUESST.slurm**
- the Slurm script that was used for running a single chain on a high performance computing cluster.


### Directions for use:

- To run this model, a user needs to download cmdStan from: 
```  https://mc-stan.org/users/interfaces/cmdstan ``` 

- Place the source code in placed in a folder called cmdstan

- The CQUESST folder that sits within our repository should be put inside the cmdstan folder

- To compile the CQUESST model you will need to have gcc available to compile C++ code.  This is achieved using 

``` module load gcc``` 

within run_CQUESST.slurm, but your computing setup may be different.  Talk to your institutes computing administrator if you have problems.
  
- Compile the Stan model from within the cmdstan directory by using 

``` make CQUESST/CQUESST```  

This is also a line that is included within run_CQUESST.slurm to automate the build process.
  
- Posterior samples can be generated using the Stan model by either running the Slurm script using 

``` sbatch run_CQUESST.slurm``` 

from the cmdstan directory.  Samples will be recorded in output.csv. Usually, this file would be run multiple times with different integer inputs for the random number seed on the last line of run_CQUESST.slurm (e.g. random seed=1 or random seed=23 etc).
  
- read_results.R can be used to extract the samples from output.csv and to create plots.