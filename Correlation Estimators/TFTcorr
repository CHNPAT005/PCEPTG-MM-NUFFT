## Author: Patrick Chang
# Script file for the MM Complex Fourier Transform
# Supporting Algorithms are at the start of the script
#  Include:
#           - Scale function to re-scale time to [0, 2 \pi]
# Number of Fourier Coefficients automatically chosen so that events
# are not aliased

#---------------------------------------------------------------------------

### Data Format:
## p = [n x 2] matrix of prices, log returns are computed in the function
# non-trading times are indicated by NaNs
## t = [n x 2] matrix of trading times, non-trading times are indicated by NaNs
# dimensions of p and t must match.

#---------------------------------------------------------------------------

using ArgCheck
using LinearAlgebra

#---------------------------------------------------------------------------
### Supporting functions

function scale(t)
    maxt = maximum(filter(!isnan, t))
    mint = minimum(filter(!isnan, t))

    tau = (2*pi) .* (t .- mint) ./ (maxt - mint)
    return tau
end

function scaleV2(t1, t2)
    maxt = maximum(filter(!isnan, [t1;t2]))
    mint = minimum(filter(!isnan, [t1;t2]))

    tau1 = (2*pi) .* (t1 .- mint) ./ (maxt - mint)
    tau2 = (2*pi) .* (t2 .- mint) ./ (maxt - mint)
    return tau1, tau2
end

#---------------------------------------------------------------------------

function TFTcorr(p, t; kwargs...)
    ## Pre-allocate arrays and check Data
    np = size(p)[1]
    mp = size(p)[2]
    nt = size(t)[1]

    @argcheck size(p) == size(t) DimensionMismatch

    # Re-scale trading times
    tau = scale(t)
    # Computing minimum time change
    # minumum step size to avoid smoothing
    dtau = diff(filter(!isnan, tau))
    taumin = minimum(filter((x) -> x>0, dtau))
    taumax = 2*pi
    # Sampling Freq.
    N0 = taumax/taumin

    # Optional Cutoff - if not specified we use Nyquist Cutoff
    kwargs = Dict(kwargs)

    if haskey(kwargs, :N)
        k = collect(1:1:kwargs[:N])
    else
        #k = collect(1:1:round(2*N0))
        k = collect(1:1:floor(N0/2))
    end

    Den = length(k)

    #------------------------------------------------------
    fca = zeros(Float64, mp, Den)
    fcb = zeros(Float64, mp, Den)

    for i in 1:mp
        psii = findall(!isnan, p[:,i])
        P = p[psii, i]
        Time = tau[psii, i]
        P = log.(P)
        D = (P[end] - P[1]) / pi    # drift
        #DiffP = diff(log.(P))

        CSA = cos.(Time[1:(end-1),:] * k') - cos.(Time[2:end,:] * k')
        CSB = sin.(Time[1:(end-1),:] * k') - sin.(Time[2:end,:] * k')

        fca[i,:] = D .+ (1/pi) .* P[1:(end-1),:]' * CSA
        fcb[i,:] = (1/pi) .* P[1:(end-1),:]' * CSB
    end

    Sigma = (pi^2 / Den) .* (fca * fca' + fcb * fcb')

    Sigma = real(Sigma)
    var = diag(Sigma)
    sigma = sqrt.(var)
    rho = Sigma ./ (sigma * sigma')
    Sigma_rescaled = 1 / N0 .* Sigma

    #dict = Dict("correlation" => rho, "variance" => var)
    return rho, Sigma, var
end

d = TFTcorr(P, t)

#---------------------------------------------------------------------------

function TFTcorrV2(p1, p2, t1, t2; kwargs...)
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
    taumin = minimum([dtau1; dtau2])
    taumax = 2*pi
    # Sampling Freq.
    N0 = taumax/taumin

    # Optional Cutoff - if not specified we use Nyquist Cutoff
    kwargs = Dict(kwargs)

    if haskey(kwargs, :N)
        k = collect(-kwargs[:N]:1:kwargs[:N])
    else
        k = collect(1:1:round(2*N0))
    end

    Den = length(k)

    #------------------------------------------------------
    fca = zeros(Float64, 2, Den)
    fcb = zeros(Float64, 2, Den)

    DP = [(log.(p1)), (log.(p2))]
    for i in 1:2
        Time = tau[i]
        P = DP[i]
        D = (P[end] - P[1]) / pi    # drift
        #DiffP = diff(log.(P))

        CSA = cos.(Time[1:(end-1),:] * k') - cos.(Time[2:end,:] * k')
        CSB = sin.(Time[1:(end-1),:] * k') - sin.(Time[2:end,:] * k')

        fca[i,:] = D .+ (1/pi) .* P[1:(end-1),:]' * CSA
        fcb[i,:] = (1/pi) .* P[1:(end-1),:]' * CSB
    end

    Sigma = (pi^2 / Den) .* (fca * fca' + fcb * fcb')

    Sigma = real(Sigma)
    var = diag(Sigma)
    sigma = sqrt.(var)
    rho = Sigma ./ (sigma * sigma')
    Sigma_rescaled = 1 / N0 .* Sigma

    #dict = Dict("correlation" => rho, "variance" => var)
    return rho, Sigma, var
end

#---------------------------------------------------------------------------
# Complex Exp implementaion with cos and sin - not the true trig implementaion
function TFTcorrV3(p, t; kwargs...)
    ## Pre-allocate arrays and check Data
    np = size(p)[1]
    mp = size(p)[2]
    nt = size(t)[1]

    @argcheck size(p) == size(t) DimensionMismatch

    # Re-scale trading times
    tau = scale(t)
    # Computing minimum time change
    # minumum step size to avoid smoothing
    dtau = diff(filter(!isnan, tau))
    taumin = minimum(filter((x) -> x>0, dtau))
    taumax = 2*pi
    # Sampling Freq.
    N0 = taumax/taumin

    # Optional Cutoff - if not specified we use Nyquist Cutoff
    kwargs = Dict(kwargs)

    if haskey(kwargs, :N)
        k = collect(-kwargs[:N]:1:kwargs[:N])
    else
        #k = collect(1:1:round(2*N0))
        #k = collect(-round(2*N0):1:round(2*N0))
        k = [collect(-round(2*N0):1:-1); collect(1:1:round(2*N0))]
        ## first and third give same answer
        #k = collect(0:1:round(2*N0))
    end

    Den = length(k)

    #------------------------------------------------------
    c_pos = zeros(ComplexF64, mp, Den)
    c_neg = zeros(ComplexF64, mp, Den)

    for i in 1:mp
        psii = findall(!isnan, p[:,i])
        P = p[psii, i]
        Time = tau[psii, i]
        DiffP = diff(log.(P))

        #c_pos[i,:] = DiffP' * exp.(1im * Time[1:(end-1),:] * k')
        #c_neg[i,:] = DiffP' * exp.(-1im * Time[1:(end-1),:] * k')
        c_pos[i,:] = DiffP' * (cos.(Time[1:(end-1),:] * k')
                                + 1im .* sin.(Time[1:(end-1),:] * k'))
        c_neg[i,:] = DiffP' * (cos.(-Time[1:(end-1),:] * k')
                                + 1im .* sin.(-Time[1:(end-1),:] * k'))
    end

    Sigma = zeros(ComplexF64, mp, mp)

    Sigma[1, 1] = 1 / (Den) * (sum(c_pos[1,:]
                                .* c_neg[1,:]))
    Sigma[1, 2] = 1 / (Den) * (sum(c_pos[1,:]
                                .* c_neg[2,:]))
    Sigma[2, 1] = Sigma[1, 2]
    Sigma[2, 2] = 1 / (Den) * (sum(c_pos[2,:]
                                .* c_neg[2,:]))

    Sigma = real(Sigma)
    var = diag(Sigma)
    sigma = sqrt.(var)
    rho = Sigma ./ (sigma * sigma')
    Sigma_rescaled = 1 / N0 .* Sigma

    #dict = Dict("correlation" => rho, "variance" => var)
    return rho, Sigma, var
end
