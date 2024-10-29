# Fire simulation

Developing and fiiting fire simulation models for north-western Patagonia. The models are aimed to be included in a dynamic vegetation model. Three models compose the fire module: ignition, escape and spread. Most of the code is related to the spread process, which is the most complex.  
  
As simulating fire spread is computationally expensive we wrote the spread functions in ```C++``` and use the ```Rcpp``` package to make it available in ```R```. These functions are hosted in the [```FireSpread```](https://github.com/barberaivan/FireSpread.git) package.  
  
We use large files storing the landscape data, which are gitignored. The ```data``` folder is available on request to ivanbarbera93@gmail.com, and should be placed inside the cloned repo to work (e.g.: ```/path_to_the_repo/fire_spread/data```).  
  
### Repository structure

```files/``` contains mostly large R objects, such as landscapes for spread, or expensive simulation outputs. This folder is also gitignored.  
  
In ```flammability indices/```, the vegetation and topographic flammability indices are created, and the corresponding code does not depend on any other. Similarly, in ```weather temporal scale/```, the lengthscale for FWI temporal correlation is fitted to fire size data (in ```FWI temporal scale.R```). Output from the code in these folders is used in the spread model, and the flammability indices are also used in the ignition and escape models.  
  
Within ```spread/```, objects used for fire spread ("landscapes") are created in ```landscapes_preparation.R ```. Fire-wise posterior distributions are sampled in ```sampling_fire_wise_posteriors_(stage1).R```, and the hierarchical model is fitted in ```hierarchical model fitting.R```. The remaining code is mostly MCMC functions and sampling trials.  
  
The code assumes ```.../fire_spread``` as the working directory, not the lower-level directories.
