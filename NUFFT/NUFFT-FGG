## Author: Patrick Chang
# Script file for the Nonuniform Fast Fourier Transform

# Fast Gaussian Gridding - code is adapted from [GL2004] Fortran code

using ArgCheck; using FFTW; using LinearAlgebra

#---------------------------------------------------------------------------

# cj = the source strength
# xj = the corresponding time for cj. xj ∈ [0, 2π]
# M = number of Fourier modes to return:
#   M even returns [-M/2, (M-1)/2]; M odd returns [-M/2, M/2]
# tol = tolerance level for error - controls the spreading

## Optional kwargs:
# R = oversampling ratio - default set to R=2

function NUFFTFGG(cj, xj, M, tol; kwargs...)

    @argcheck size(cj) == size(xj) DimensionMismatch

    nj = size(cj)[1]

    if haskey(kwargs, :R)
        R = kwargs[:R]
    else
        R = 2
    end

    Mr = R*M

    nspread = Int(floor(-log(tol)/(pi*(R-1)/(R-0.5)) + 0.5))

    r2lamb = R^2 * nspread / (R*(R-0.5))
    hx = 2*pi/Mr

    t1 = pi/r2lamb
    fw = zeros(Float64, Mr, 1)

    for k in 1:nspread
        fw[k] = exp(-t1*k^2)
    end

    ftau = zeros(ComplexF64, Mr, 1)

    for j in 1:nj
        ccj = cj[j]

        jb1 = floor(xj[j]/hx)
        diff1 = xj[j]/hx - jb1
        jb1 = Int(jb1)

        xc = zeros(Float64, 2*nspread, 1)
        xc[nspread] = exp(-t1*diff1^2)
        cross = xc[nspread]
        cross1 = exp(2*t1 * diff1)
        for k in 1:nspread
            cro = cross * cross1^k
            xc[k+nspread] = fw[k] * cro
        end
        for k in 1:(nspread-1)
            cro = cross * cross1^(-k)
            xc[-k+nspread] = fw[k] * cro
        end

        jb1d = min(nspread-1, jb1)
        jb1u = min(nspread, Mr-jb1-1)
        for k in (-nspread+1):(-jb1d-1)
            istart = jb1 + k + Mr + 1
            zz = xc[k+nspread] * ccj
            ftau[istart] = ftau[istart] + zz
        end
        for k in -jb1d:jb1u
            istart = jb1 + k + 1
            zz = xc[k+nspread] * ccj
            ftau[istart] = ftau[istart] + zz
        end
        for k in jb1u+1:nspread
            istart = jb1 + k - Mr + 1
            zz = xc[k+nspread] * ccj
            ftau[istart] = ftau[istart] + zz
        end
    end

    freq = fftfreq(Mr, 1) .* Mr
    pos = findall(abs.(freq) .<= M/2)

    Ftau = fft(ftau)[pos] ./ Mr

    k = freq[pos]

    tau = pi * r2lamb / Mr^2
    Fk = sqrt(pi/tau) .* exp.(k.^2 .* tau) .* Ftau

    return Fk
end

# # Simple test case
#
# nj = 10000
# x = (collect(0:nj-1) + 0.5 .* rand(nj))
# xj = (x .- minimum(x)) .* (2*pi / (maximum(x) - minimum(x)))
# cj = rand(nj) + 0im*rand(nj)
#
# M = 5001
# Mr = 2*M
# freq = fftfreq(Mr, 1) .* Mr
# pos = findall(abs.(freq) .<= M/2)
# k = freq[pos]
#
# test = (cj' * exp.(1im .* xj * k'))'
# test2 = NUFFTFGG(cj, xj, M, 10^-12)
#
# # Check accruacy with the L2 norm
# norm(test-test2) / norm(test)
