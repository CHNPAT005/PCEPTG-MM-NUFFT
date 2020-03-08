# Malliavin-Mancino estimators implemented with non-uniform fast Fourier transforms

## Authors:
- Patrick Chang
- Etienne Pienaar
- Tim Gebbie

## Link to resources:

Link to paper:

Link to the Dataset: [Link]{https://zivahub.uct.ac.za/articles/Malliavin-Mancino_estimators_implemented_with_the_non-uniform_fast_Fourier_transform_Dataset/11903442}

## Steps for Replication:
- Change directories for all the files under [/Scripts](https://github.com/CHNPAT005/PCEPTG-MM-NUFFT/tree/master/Scripts). Currently the directories are set as: `cd("/Users/patrickchang1/PCEPTG-MM-NUFFT")`. Change this to where you have stored the file `PCEPTG-MM-NUFFT`. 

	- Run [/Scripts/Timing/Dirichlet Timing](https://github.com/CHNPAT005/PCEPTG-MM-NUFFT/blob/master/Scripts/Timing/Dirichlet\%20Timing) and [/Scripts/Timing/Fejer Timing](https://github.com/CHNPAT005/PCEPTG-MM-NUFFT/blob/master/Scripts/Timing/Fejer\%20Timing) to reproduce Fig. 4.
	 Run [/Scripts/Timing/Error Timing](https://github.com/CHNPAT005/PCEPTG-MM-NUFFT/blob/master/Scripts/Timing/Error\%20Timing) to reproduce Fig. 5.
	 - Run [/Scripts/Accuracy/AccSynDS](https://github.com/CHNPAT005/PCEPTG-MM-NUFFT/blob/master/Scripts/Accuracy/AccSynDS) and [/Scripts/Accuracy/AccRE](https://github.com/CHNPAT005/PCEPTG-MM-NUFFT/blob/master/Scripts/Accuracy/AccRE) to reproduce Figs. 6 and 7.
	 - Run [/Scripts/Time Scales/MMZandMM](https://github.com/CHNPAT005/PCEPTG-MM-NUFFT/blob/master/Scripts/Time\%20Scales/MMZandMM) to reproduce Fig. 8.
	 - Run [/Scripts/Time Scales/3Asset](https://github.com/CHNPAT005/PCEPTG-MM-NUFFT/blob/master/Scripts/Time\%20Scales/3Asset) to reproduce Fig. C.12.
 
 - To reproduce the Empirical analysis - download the processed dataset from ZivaHub and put the two csv files into the folder `/Real Data`.
	 - Run [/Scripts/Time Scales/Empirical](https://github.com/CHNPAT005/PCEPTG-MM-NUFFT/blob/master/Scripts/Time\%20Scales/Empirical) to reproduce Figs. 9, B.10 and B.11.
	 
- We have included the plots under `/Plots` and Computed results under `/Computed Data` if one does not wish to re-run everything.

## Using the functions for other purposes:
### NUFFT

We have included 1 Dimensional Type 1 non-uniform fast Fourier transforms, implemented using three different kernels: Gaussian, Kaiser-Bessel and exponential of semi-circle under [/Functions/NUFFT](https://github.com/CHNPAT005/PCEPTG-MM-NUFFT/tree/master/Functions/NUFFT).

The functions require 4 input variables:
- cj: vector of source strength
- xj: vector of source points
- M: the number of Fourier coefficients you want returned (integer)
- tol: the precision requested

#### Example

```julia

include("Functions/NUFFT/NUFFT-FGG")

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

We implement Malliavin-Mancino estimators which include the [/Functions/Correlation Estimators/Dirichlet](https://github.com/CHNPAT005/PCEPTG-MM-NUFFT/tree/master/Functions/Correlation\%20Estimators/Dirichlet) basis kernel and the [/Functions/Correlation Estimators/Fejer](https://github.com/CHNPAT005/PCEPTG-MM-NUFFT/tree/master/Functions/Correlation\%20Estimators/Fejer) basis kernel. 

Implementation methods include the ``for-loop'' implementation (MS), vectorised implementation (CFT), fast Fourier transform implementation (FFT), zero-padded fast Fourier transform implementation (FFTZP) and the non-uniform fast Fourier transform (NUFFT). 

The functions require 2 input variables:
- p: (nxm) matrix of prices, with non-trade times represented as NaNs.
- t: (nxm) corresponding matrix of trade times, with non-trade times represented as NaNs.
and two optional input variables.
- N: (optional input) for the number of Fourier coefficients used in the convolution of the Malliavin-Mancino estimator (integer) - defaults to the Nyquist frequency.
- tol: tolerance requested - applies only to NUFFT implementations - defaults to 10^-12.

Note that the zero-padded FFT implementation has no optional argument for N, and the FFT implementation takes only 1 input variable, P.

#### Example

```julia

include("Functions/Correlation Estimators/Dirichlet/NUFFTcorrDK-FGG")
include("Functions/Monte Carlo Simulation Algorithms/GBM")

# Create some data
mu = [0.01/86400, 0.01/86400]
sigma = [0.1/86400 sqrt(0.1/86400)*0.35*sqrt(0.2/86400);
        sqrt(0.1/86400)*0.35*sqrt(0.2/86400) 0.2/86400]

P = GBM(10000, mu, sigma, seed = 10)
t = reshape([collect(1:1:10000.0); collect(1:1:10000.0)], 10000, 2)

# Obtain results
output1 = NUFFTcorrDKFGG(P, t) # Nyquist, tol = 1e-12
output2 = NUFFTcorrDKFGG(P, t, N = 500, tol = 10^-8) 

# Extract results
cor1 = output1[1]
cov1 = output1[2]

cor2 = output2[1]
cov2 = output2[2]

```

