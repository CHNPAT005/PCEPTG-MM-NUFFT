## Author: Patrick Chang
# Script file for the MM Complex Fourier Transform (Dirichlet)
# Supporting Algorithms are at the start of the script
#  Include:
#           - Scale function to re-scale time to [0, 2 \pi]
# Number of Fourier Coefficients automatically chosen so that events
# are not aliased

#---------------------------------------------------------------------------

### Data Format:
## p1;p2 = vector of stock prices - can include or not include NaNs
## t1;t2 = vector of trading times - can include or not include NaNs
## N = Optional input for cutoff frequency

#---------------------------------------------------------------------------

using ArgCheck; using LinearAlgebra

#---------------------------------------------------------------------------
### Supporting functions

function scaleV2(t1, t2)
    maxt = maximum(filter(!isnan, [t1;t2]))
    mint = minimum(filter(!isnan, [t1;t2]))

    tau1 = (2*pi) .* (t1 .- mint) ./ (maxt - mint)
    tau2 = (2*pi) .* (t2 .- mint) ./ (maxt - mint)
    return tau1, tau2
end

function KanataniDK(t1, t2, N)
    t1 = filter(!isnan, t1)
    t2 = filter(!isnan, t2)
    L1 = length(t1)
    L2 = length(t2)

    W = zeros(L1-1, L2-1)

    for i in 1:(L1-1)
        for j in 1:(L2-1)
            if (t1[i] == t2[j])
                W[i,j] = 1
            else
                W[i,j] = (1/(2*N+1)) * sin((N + 0.5) * (t1[i] - t2[j]))/sin((t1[i] - t2[j])/2)
            end
        end
    end
    return W
end

#---------------------------------------------------------------------------

# Kanatani Weight matrix implementaion using the Dirichlet Kernel

function KANcorrDK(p1, p2, t1, t2; kwargs...)
    t1 = filter(!isnan, t1)
    t2 = filter(!isnan, t2)

    p1 = filter(!isnan, p1)
    p2 = filter(!isnan, p2)

    # Re-scale trading times
    tau = scaleV2(t1, t2)
    tau1 = tau[1]
    tau2 = tau[2]
    # Computing minimum time change
    # minumum step size to avoid smoothing
    dtau1 = diff(tau1)
    dtau2 = diff(tau2)
    taumin = maximum([dtau1; dtau2])
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

    DiffP1 = diff(log.(p1))
    DiffP2 = diff(log.(p2))

    W1 = KanataniDK(tau1, tau1, N)
    W2 = KanataniDK(tau2, tau2, N)
    W12 = KanataniDK(tau1, tau2, N)

    Sigma = zeros(Float64, 2, 2)
    Sigma[1, 1] = DiffP1' * W1 * DiffP1
    Sigma[2, 2] = DiffP2' * W2 * DiffP2
    Sigma[1, 2] = DiffP1' * W12 * DiffP2
    Sigma[2, 1] = Sigma[1, 2]

    var = diag(Sigma)
    sigma = sqrt.(var)
    rho = Sigma ./ (sigma * sigma')

    return rho, Sigma, var
end
