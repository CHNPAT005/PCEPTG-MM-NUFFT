## Author: Patrick Chang
# Script file to simulate a multivariate Geometric Brownian Motion
# includes an example to demonstrate the code

#------------------------------------------------------------------------------

## Code for Geometric Brownian Motion using the method from Paul Glasserman
#  in his book - Monte Carlo Methods in Financial Engineering

# timesteps (dt) can be controlled by scaling mu and sigma accordingly

using Random; using LinearAlgebra; using Plots
# https://docs.juliaplots.org/latest/

function GBM(n, mu, sigma; kwargs...)
    # n - simlulation length
    # mu - vector input of the drift component
    # sigma - covariance matrix of the stocks
    # startprice - starting price for the assets

    # all inputs must have appropriate dimensions

    k = size(sigma)[1]

    kwargs = Dict(kwargs)

    if haskey(kwargs, :startprice)
        startprice = kwargs[:startprice]
    else
        startprice = fill(100.0, (k,1))
    end

    if haskey(kwargs, :seed)
        seed = kwargs[:seed]
    else
        seed = 1
    end

    mu = reshape(mu, k, 1)
    sigma = reshape(sigma, k, k)
    sigma2 = reshape(diag(sigma), k, 1)

    P = zeros(n, k)
    P[1,:] = startprice

    A = cholesky(sigma).L
    b = mu - sigma2./2

    Random.seed!(seed)
    Z = randn(k, n-1)

    for i in 2:n
        z = Z[:,i-1]
        X = b + A * z
        P[i,:] = P[i-1,:] .* exp.(X)
    end
    return P
end

function GBM2(n, mu, sigma; kwargs...)
    # n - simlulation length
    # mu - vector input of the drift component
    # sigma - covariance matrix of the stocks
    # startprice - starting price for the assets

    # all inputs must have appropriate dimensions

    k = size(sigma)[1]

    kwargs = Dict(kwargs)

    if haskey(kwargs, :startprice)
        startprice = kwargs[:startprice]
    else
        startprice = fill(100.0, (k,1))
    end

    if haskey(kwargs, :seed)
        seed = kwargs[:seed]
    else
        seed = 1
    end

    if haskey(kwargs, :dt)
        dt = kwargs[:dt]
    else
        dt = 1
    end

    mu = reshape(mu, k, 1)
    sigma = reshape(sigma, k, k)
    sigma2 = reshape(diag(sigma), k, 1)

    P = zeros(n, k)
    P[1,:] = startprice

    A = cholesky(sigma).L
    b = mu - sigma2./2

    Random.seed!(seed)
    Z = randn(k, n-1)

    for i in 2:n
        z = Z[:,i-1]
        X = b.*dt + (A * z).*sqrt(dt)
        P[i,:] = P[i-1,:] .* exp.(X)
    end
    return P
end

#------------------------------------------------------------------------------

## Uncomment for plots

# n = 500
# mu = [0.01/86400, 0.01/86400]
# sigma = [0.1/86400 sqrt(0.1/86400)*0.35*sqrt(0.2/86400);
#         sqrt(0.1/86400)*0.35*sqrt(0.2/86400) 0.2/86400]
#
#
# P = GBM(n, mu, sigma)
# P2 = GBM(n, mu, sigma, seed = 2)
#
# # important thing to note about plotting is that [mxn] matrix is m series with n data points
# p1 = plot(1:n, P, label = ["Line 1" "Line 2"])
# title!(p1, "Correlated GBM")
# xlabel!(p1, "Seconds")
# ylabel!(p1, "Price")
#
# p2 = plot(1:n, P2, label = ["Line 1" "Line 2"])
# title!(p2, "Correlated GBM")
# xlabel!(p2, "Seconds")
# ylabel!(p2, "Price")
#
# plot(p1, p2, layout = (1,2), legend = false, lw = 2)
#
# #plt = plot(1, P[1,:], xlim = (0,500), ylim = (94, 102), title = "Correlated GBM", legend = false)
# @gif for i in 1:500
#     plot(1:i, P[1:i,:], xlim = (0,500), ylim = (94, 102), title = "Correlated GBM", legend = false)
# end every 10
