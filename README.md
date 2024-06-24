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


### Directions for use (without MPI):
If you only want to run CQUESST on a single processor, then you can use these simple instructions.
Note that, CQUESST allows for parallelisation of field-plots across separate processors using MPI which provides the ability to efficiently scale the analysis to include many sites (see below).

- To run this model, a user needs to download cmdStan from: 
```  https://mc-stan.org/users/interfaces/cmdstan ``` 

- Place the source code in placed in a folder called cmdstan

- The CQUESST folder that sits within our repository should be put inside the cmdstan folder

- To compile the CQUESST model you will need to have gcc available to compile C++ code.  This is achieved using 

``` module load gcc``` 

within run_CQUESST.slurm, but your computing setup may be different.  Talk to your institutes computing administrator if you have problems.
  
- Compile the Stan model from within the cmdstan directory by using 

``` make CQUESST/CQUESST```  
  
- Posterior samples can be generated using the Stan model by running the Slurm script using 

``` sbatch run_CQUESST.slurm``` 

from the cmdstan directory.  Samples will be recorded in output.csv. Usually, this file would be run multiple times with different integer inputs for the random number seed on the last line of run_CQUESST.slurm (e.g. random seed=1 or random seed=23 etc).
  
- read_results.R can be used to extract the samples from output.csv and to create plots.


### Directions for use (with MPI):
MPI allows cmdstan to scale CQUESST across many processors across multiple nodes of a high-performance computing cluster.

- To run this model, a user needs to download cmdStan from: 
```  https://mc-stan.org/users/interfaces/cmdstan ``` 

- Place the source code in placed in a folder called cmdstan

- The CQUESST folder that sits within our repository should be put inside the cmdstan folder

- To compile the CQUESST model you will need to have gcc available to compile C++ code.  This is achieved using 

``` module load git gcc openmpi nano gcc  gtest libxml2 python```

- Next, create or edit the ./make/local file in the cmdstan directory using:

``` nano make/local ```

- Add the following lines to the make/local file

``` 
STAN_MPI=true
CXX=mpicxx
LDLIBS += -lpthread
TBB_CXX_TYPE=gcc
STAN_HAS_CXX17=true 
```

- Next, we need a couple of paths for configuring Boost.  Type the following two commands:

``` 
which mpicxx
which python
 ```

- Find and modify the file at <path_to>/cmdstan/stan/lib/stan_math/lib/boost_<version>/user-config.jam to have the following statements (inserting the paths that you got from the previous two paths as appropriate).

```
using mpi : /apps/openmpi/4.1.1-ofed51/bin/mpicxx ;
using python : 3.9 : /apps/python/3.9.4/bin/python : /apps/python/3.9.4/include/python3.9 ; 
```

(NOTE, the spaces above are important).

- Now change to the cmdstan directory and build the MPI and Stan components using

```
make clean-all
make -j8 build-mpi
make -j8 build
```

- Once that completes, you can build the CQUESST cmdstan model using

``` make CQUESST/CQUESST ```
  
- Posterior samples can be generated using the Stan model by either running the Slurm script using 

``` sbatch run_CQUESST_mpi.slurm``` 

from the cmdstan directory.  Samples will be recorded in output.csv. Usually, this file would be run multiple times with different integer inputs for the random number seed on the last line of run_CQUESST.slurm (e.g. random seed=1 or random seed=23 etc).
  
- read_results.R can be used to extract the samples from output.csv and to create plots.

