## Author: Patrick Chang
# Script file to generate random correlation/covariance matrices

function gencormatrix(d)
    A = rand(d,d);

    A = A*A';

    A = A + d.*Matrix(I, d, d)

    var = diag(A)
    sigma = sqrt.(var)
    rho = A ./ (sigma * sigma')
    return rho
end

function gencovmatrix(d, σ)
    length(σ) != d && throw(DimensionMismatch("length(σ) doesn't match d"))
    r = gencormatrix(d)
    Σ = zeros(d, d)
    @inbounds for i in 1:d, j in 1:d
        Σ[i,j] = r[i,j] * σ[i] * σ[j]
    end
    return Σ
end
