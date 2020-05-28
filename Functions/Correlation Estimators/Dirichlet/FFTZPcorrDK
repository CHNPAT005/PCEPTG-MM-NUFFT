## Author: Patrick Chang
# Script file for the MM Fast Fourier Transform with Zero Padding
# Supporting Algorithms are at the start of the script
#  Include:
#           - Scale function to re-scale time to [0, 2 \pi]
#           - Zero-padding function
# Number of Fourier Coefficients used is length of data

#---------------------------------------------------------------------------

### Data Format:
## p = [n x m] matrix of prices, log returns are computed in the function
# non-trading times are indicated by NaNs
## t = [n x m] matrix of trading times, non-trading times are indicated by NaNs
# dimensions of p and t must match.

#---------------------------------------------------------------------------

using ArgCheck; using LinearAlgebra;
using FFTW

#---------------------------------------------------------------------------
### Supporting functions

function scale(t)
    maxt = maximum(filter(!isnan, t))
    mint = minimum(filter(!isnan, t))

    tau = (2*pi) .* (t .- mint) ./ (maxt - mint)
    return tau
end

function zero_pad(cj, xj, N0)
    res = zeros(N0, 1)
    nj = length(xj)
    for j in 1:nj
        pos = Int(round(xj[j] * (N0/(2*pi))))
        res[pos+1] = cj[j]
    end
    return res
end
#ZPP = zero_pad(DiffP, tau[psii, 1], Int(ceil(N0)))
#---------------------------------------------------------------------------

# Fast Fourier Transform (with Zero Padding) implementaion of the Dirichlet Kernel
# Can be used for asynchronous data
# No option for optional cutoff (N)

function FFTZPcorrDK(p, t)
    ## Pre-allocate arrays and check Data
    np = size(p)[1]
    mp = size(p)[2]
    nt = size(t)[1]

    #@argcheck size(p) == size(t) DimensionMismatch

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

    k = collect(-floor(N0/2):1:floor(N0/2))

    Den = length(k)

    #------------------------------------------------------
    c_pos = zeros(ComplexF64, mp, Den)
    c_neg = zeros(ComplexF64, mp, Den)

    for i in 1:mp
        psii = findall(!isnan, p[:,i])
        P = p[psii, i]
        DiffP = diff(log.(P))
        psii = psii[1:(end-1)]
        Time = tau[psii, i]
        ZP = zeros(Den, 1)

        ZP[1:Int(floor(N0))] = zero_pad(DiffP, Time, Int(floor(N0)))

        C = fft(ZP)

        C_pos = C
        C_neg = conj(C)

        c_pos[i,:] = C_pos
        c_neg[i,:] = C_neg
    end

    Sigma = zeros(ComplexF64, mp, mp)

    Sigma = 0.5 / Den .* (c_pos*c_pos' + c_neg*c_neg')

    Sigma = real(Sigma)
    var = diag(Sigma)
    sigma = sqrt.(var)
    rho = Sigma ./ (sigma * sigma')

    return rho, Sigma, var
end
