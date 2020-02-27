## Author: Patrick Chang
# Script file for the MM Fast Fourier Transform
# Supporting Algorithms are at the start of the script
#  Include:
#           - Scale function to re-scale time to [0, 2 \pi]
# Number of Fourier Coefficients used is length of data

#---------------------------------------------------------------------------

### Data Format:
## p = [n x m] matrix of prices, log returns are computed in the function
# non-trading times are indicated by NaNs

#---------------------------------------------------------------------------

using ArgCheck; using LinearAlgebra;
using FFTW

#---------------------------------------------------------------------------

# Fast Fourier Transform implementaion of the Fejer Kernel
# Can ONLY be used for fully synchronous data
# No option for optional cutoff (N)

function FFTcorrFK(p)
    ## Pre-allocate arrays and check Data
    np = size(p)[1]
    mp = size(p)[2]
    #nt = size(t)[1]

    #@argcheck size(p) == size(t) DimensionMismatch

    # # Re-scale trading times
    # tau = scale(t)
    # # Computing minimum time change
    # # minumum step size to avoid smoothing
    # dtau = diff(filter(!isnan, tau))
    # taumin = minimum(filter((x) -> x>0, dtau))
    # taumax = 2*pi
    # # Sampling Freq.
    # N0 = taumax/taumin

    if (iseven(np))
        k = collect(0:1:(np-1)/2)
    else
        k = collect(0:1:(np+1)/2)
    end

    Den = length(k)

    #------------------------------------------------------
    c_pos = zeros(ComplexF64, mp, Den + (Den-1))
    c_neg = zeros(ComplexF64, mp, Den + (Den-1))

    for i in 1:mp
        DiffP = diff(log.(p[:,i]))

        C = rfft(DiffP)

        C_pos = C
        C_neg = conj(C)

        c_pos[i,:] = [C_pos; C_neg[2:end]]
        c_neg[i,:] = [C_neg; C_pos[2:end]]
    end

    k = [k; k[2:end]]
    N = Den-1

    Sigma = zeros(ComplexF64, mp, mp)
    for i in 1:mp-1
        for j in i+1:mp
            Sigma[i,i] = sum( (1 .- abs.(k)./N) .* c_pos[i,:] .* c_neg[i,:] ) / (N+1)
            Sigma[j,j] = sum( (1 .- abs.(k)./N) .* c_pos[j,:] .* c_neg[j,:] ) / (N+1)
            Sigma[i,j] = Sigma[j,i] = sum( (1 .- abs.(k)./N) .* c_pos[i,:] .* c_neg[j,:] ) / (N+1)
        end
    end

    Sigma = real(Sigma)
    var = diag(Sigma)
    sigma = sqrt.(var)
    rho = Sigma ./ (sigma * sigma')

    return rho, Sigma, var
end
