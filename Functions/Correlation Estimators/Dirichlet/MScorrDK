## Author: Patrick Chang
# Script file for the MM Complex Fourier Transform (Dirichlet)
# Supporting Algorithms are at the start of the script
#  Include:
#           - Scale function to re-scale time to [0, 2 \pi]
# Number of Fourier Coefficients automatically chosen so that events
# are not aliased

#---------------------------------------------------------------------------

### Data Format:
## p = [n x m] matrix of prices, log returns are computed in the function
# non-trading times are indicated by NaNs
## t = [n x m] matrix of trading times, non-trading times are indicated by NaNs
# dimensions of p and t must match.
## N = Optional input for cutoff frequency

#---------------------------------------------------------------------------

using ArgCheck; using LinearAlgebra

#---------------------------------------------------------------------------
### Supporting functions

function scale(t)
    maxt = maximum(filter(!isnan, t))
    mint = minimum(filter(!isnan, t))

    tau = (2*pi) .* (t .- mint) ./ (maxt - mint)
    return tau
end

#---------------------------------------------------------------------------

# Mancino Sanfelici implementation using for loops
# with the Dirichlet Kernel

function MScorrDK(p, t; kwargs...)
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
        N = kwargs[:N]
    else
        N = trunc(Int, floor(N0/2))
    end

    #------------------------------------------------------

    c_pos = zeros(ComplexF64, mp, 2*N + 1)
    c_neg = zeros(ComplexF64, mp, 2*N + 1)

    for i in 1:mp
        psii = findall(!isnan, p[:,i])
        P = p[psii, i]
        Time = tau[psii, i]
        DiffP = diff(log.(P))

        for k=1:(2*N+1)
            s=k-N-1;
            c_neg[i, k]=sum(exp.((-1im*s).*Time[1:(end-1)]).*DiffP);
            c_pos[i, k]=sum(exp.((1im*s).*Time[1:(end-1)]).*DiffP);
        end
    end

    Sigma = 0.5 / (2*N + 1) .* (c_pos*c_pos' + c_neg*c_neg')

    Sigma = real(Sigma)
    var = diag(Sigma)
    sigma = sqrt.(var)
    rho = Sigma ./ (sigma * sigma')

    return rho, Sigma, var
end
