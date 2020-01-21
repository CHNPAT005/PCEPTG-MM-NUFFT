# High speed implementation of the Malliavin-Mancino estimator

### Authors:
- Patrick Chang
- Etienne Pienaar
- Tim Gebbie

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
