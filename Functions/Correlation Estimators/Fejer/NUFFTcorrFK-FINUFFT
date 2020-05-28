## Author: Patrick Chang
# Script file for the MM NUFFT
# Supporting Algorithms are at the start of the script
#  Include:
#           - Scale function to re-scale time to [0, 2 \pi]
# Number of Fourier Coefficients used is length of data

## Implementation uses FINUFFT package

#---------------------------------------------------------------------------

### Data Format:
## p = [n x m] matrix of prices, log returns are computed in the function
# non-trading times are indicated by NaNs
## t = [n x m] matrix of trading times, non-trading times are indicated by NaNs
# dimensions of p and t must match.
## N = Optional input for cutoff frequency
## tol = error tolerance for NUFFT - determines how much spreading

#---------------------------------------------------------------------------

using ArgCheck; using LinearAlgebra; using FINUFFT

#---------------------------------------------------------------------------
### Supporting functions

function scale(t)
    maxt = maximum(filter(!isnan, t))
    mint = minimum(filter(!isnan, t))

    tau = (2*pi) .* (t .- mint) ./ (maxt - mint)
    return tau
end

#---------------------------------------------------------------------------

# Non-uniform Fast Fourier Transform implementaion of the Fejer Kernel

function NUFFTcorrFKFINUFFT(p, t; kwargs...)
    ## Pre-allocate arrays and check Data
    np = size(p)[1]
    mp = size(p)[2]
    nt = size(t)[1]

    @argcheck size(p) == size(t) DimensionMismatch

    # Re-scale trading times
    tau = scale(t)
    # Computing minimum time change
    dtau = zeros(mp,1)
    for i in 1:mp
        dtau[i] = minimum(diff(filter(!isnan, tau[:,i])))
    end
    # maximum of minumum step size to avoid aliasing
    taumin = maximum(dtau)
    taumax = 2*pi
    # Sampling Freq.
    N0 = taumax/taumin

    # Optional Cutoff - if not specified we use Nyquist Cutoff
    kwargs = Dict(kwargs)

    if haskey(kwargs, :N)
        k = collect(-kwargs[:N]:1:kwargs[:N])
    else
        k = collect(-floor(N0/2):1:floor(N0/2))
    end

    if haskey(kwargs, :tol)
        tol = kwargs[:tol]
    else
        tol = 10^-12
    end

    Den = length(k)

    #------------------------------------------------------
    e_pos = zeros(ComplexF64, mp, Den)
    e_neg = zeros(ComplexF64, mp, Den)

    for i in 1:mp
        psii = findall(!isnan, p[:,i])
        P = p[psii, i]
        Time = tau[psii, i]
        DiffP = complex(diff(log.(P)))
        Time = Time[1:(end-1)]

        C = nufft1d1(Time, DiffP, -1, tol, Den)

        e_pos[i,:] = C
        e_neg[i,:] = conj(C)
    end

    N = (Den-1)/2

    Sigma = zeros(ComplexF64, mp, mp)
    for i in 1:mp-1
        for j in i+1:mp
            Sigma[i,i] = sum( (1 .- abs.(k)./N) .* e_pos[i,:] .* e_neg[i,:] ) / (N+1)
            Sigma[j,j] = sum( (1 .- abs.(k)./N) .* e_pos[j,:] .* e_neg[j,:] ) / (N+1)
            Sigma[i,j] = Sigma[j,i] = sum( (1 .- abs.(k)./N) .* e_pos[i,:] .* e_neg[j,:] ) / (N+1)
        end
    end

    Sigma = real(Sigma)
    var = diag(Sigma)
    sigma = sqrt.(var)
    rho = Sigma ./ (sigma * sigma')

    return rho, Sigma, var
end
