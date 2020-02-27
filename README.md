# High speed implementation of the Malliavin-Mancino estimator

## Authors:
- Patrick Chang
- Etienne Pienaar
- Tim Gebbie

## Link to the paper:


## Steps for Replication:
- Change Directories for the NUFFTcorr functions found under `/Correlation Estimators/Dirichlet` and `/Correlation Estimators/Fejer`. Currently the directories are set as: `cd("/Users/patrickchang1/PCEPTG-MM-NUFFT")`. Change this to where you have stored the file `PCEPTG-MM-NUFFT`. 
- Change the Directories for the files which reproduce the required results.
	- Run [/Timing/Dirichlet Timing](https://github.com/CHNPAT005/PCEPTG-MM-NUFFT/blob/master/Timing/Dirichlet\%20Timing) and [/Timing/Fejer Timing](https://github.com/CHNPAT005/PCEPTG-MM-NUFFT/blob/master/Timing/Fejer\%20Timing) to reproduce Fig. 4.
	- Run [/Timing/Error Timing](https://github.com/CHNPAT005/PCEPTG-MM-NUFFT/blob/master/Timing/Error\%20Timing) to reproduce Fig. 5.
	- Run [/Accuracy/AccSynDS](https://github.com/CHNPAT005/PCEPTG-MM-NUFFT/blob/master/Accuracy/AccSynDS) and [/Accuracy/AccRE](https://github.com/CHNPAT005/PCEPTG-MM-NUFFT/blob/master/Accuracy/AccRE) to reproduce Figs. 6 and 7.
	- Run [/Time Scales/MMZandMM](https://github.com/CHNPAT005/PCEPTG-MM-NUFFT/blob/master/Time\%20Scales/MMZandMM) to reproduce Fig. 8.
- To reproduce the Empirical analysis - download the processed dataset from DOI: 10.25375/uct.11903442 and put the datasets into the folder `/Time Scales`.
	- Run [/Time Scales/Empirical](https://github.com/CHNPAT005/PCEPTG-MM-NUFFT/blob/master/Time\%20Scales/Empirical) to reproduce Figs. 9, C.10 and C.11.
- We have included the plots under `/Plots` and Computed results under `/Computed Data` if one does not wish to re-run everything.


## Using the functions for other purposes:
### NUFFT

We have included 1 Dimensional Type 1 non-uniform fast Fourier transforms, implemented using three different kernels: Gaussian, Kaiser-Bessel and exponential of semi-circle under [/NUFFT](https://github.com/CHNPAT005/PCEPTG-MM-NUFFT/tree/master/NUFFT).

#### Example
```julia

include("../../NUFFT/NUFFT-FGG")

# Simulate some non-uniform data
nj = 10
x = (collect(0:nj-1) + 0.5 .* rand(nj))
xj = (x .- minimum(x)) .* (2*pi / (maximum(x) - minimum(x))) 	# Re-scale s.t. xj \in [0, 2\pi]
cj = rand(nj) + 0im*rand(nj)

# Parameter settings
M = 11 # Output size
tol = 10^-12 # Tolerance

# Output 
fk = NUFFTFGG(cj, xj, M, tol)

```

### Malliavin-Mancino estimators using non-uniform FFTs



### Site Navigation:
- All the script files for the estimators are under `/Correlation Estimators` which includes:
  - Various Dirichlet implementations under `/Correlation Estimators/Dirichlet`
  - Various Fejer implementations under `/Correlation Estimators/Fejer`
  - The trig implementation from MM2002 under `/Correlation Estimators/Trig`
  - We also have the Hayashi-Yoshida implementation under `/Correlation Estimators/HY`
  
- The script files for the simulation of SDEs are under `/Monte Carlo Simulation Algorithms` 
- The script files for our implementation of the Type 1 Non-Uniform Fast Fourier Transform (NUFFT) are under `/NUFFT`. These include: Naive Gaussian (NG) [GL2004] and the Fast Gaussian Gridding (FGG) [GL2004].
- Test files are under `/Test`. These test the consistency across various implementations of the Dirichlet and Fejer estiamtors.

### TODO:
- Timing
- Code up the Kaiser-Besel kernel implementation of the NUFFT and integrate into the estimators
