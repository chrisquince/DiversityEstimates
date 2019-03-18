**Installation**


The software currently only runs on Linux computers with a recent version of the Gnu Science Library installed. To install the software:

```
git clone https://github.com/chrisquince/DiversityEstimates.git
cd DiversityEstimates
make
make install
```

This creates binaries in DiversityEstimates/bin

You may need to change the compiler flags in the makefile. In particular remove the flag "-m64" on 32 bit machines.

Each executable (MetroLogNormal, MetroIG, MetroLogStudent, MetroSichel) fits a different abundance distribution (Log-normal, inverse Gaussian, log-Student's t, Sichel) respectively.


**Usage**

Typing the executable name, e.g. "MetroIG" gives a list of command line arguments.

The directory "data" contains the sample abundance distributions fitted in the [article](https://www.ncbi.nlm.nih.gov/pubmed/18650928).

First run the executable on a sample without MCMC sampling. This will perform a maximum likelihood fit using the Simplex algorithm. For instance for the GOS data type "./MetroIG -in GOS.sample -out GOS".

Then run a short trial MCMC run of 1000 iterations with guessed std. dev.s for the proposal distributions say about 10% of the parameter values e.g. "./MetroIG -in GOS.sample -out GOS -s 1000 -sigmaA 0.05 -sigmaB 0.3 -sigmaS 200".

Adjust the std. dev.s untill the acceptance ratios are about 0.5. Then perform a longer run of say 250,000 iterations: "./MetroIG -in GOS.sample -out GOS -s 250000 -sigmaA 0.01 -sigmaB 0.1 -sigmaS 100"

Three data files with posterior samples for three different sets of parameter values will be generated, in the above case "GOS_0.sample", "GOS_1.sample", "GOS_2.sample".
