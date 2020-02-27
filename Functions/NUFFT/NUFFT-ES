## Author: Patrick Chang
# Script file for the Nonuniform Fast Fourier Transform

# Kaiser-Bessel - spreading to grid points near the source points

using ArgCheck; using FFTW; using QuadGK

# Supporting functions
#---------------------------------------------------------------------------

function ES(z, ω)
    B = 2.3*ω
    if abs(z) <= 1
        return exp(B*(sqrt(1-z^2) - 1))
    else
        return 0
    end
end

function ES_kernel(x, ω, Mr)
    α = pi * ω / Mr
    return ES(x/α, ω)
end

function ES_kernel_hat(k, ω, Mr)
    α = pi * ω / Mr
    return quadgk(x -> α * exp(1im*α*k*x) * ES(x, ω), -1, 1)[1]
end

#---------------------------------------------------------------------------

# cj = the source strength
# xj = the corresponding time for cj. xj ∈ [-π, π]
# M = number of Fourier modes to return:
#   M even returns [-M/2, (M-1)/2]; M odd returns [-M/2, M/2]
# tol = tolerance level for error - controls the spreading

function NUFFTES(cj, xj, M, tol; kwargs...)

    @argcheck size(cj) == size(xj) DimensionMismatch

    nj = size(cj)[1]

    R = 2
    Mr = R*M
    nspread = Int(floor((ceil(log(10, 1/tol)) + 2)/2) + 2)
    # nspread = Int(floor((ceil(log(10, 1/tol)) + 2)/2))
    ω = 2*nspread + 1
    hx = 2*pi/Mr

    ftau = zeros(ComplexF64, Mr, 1)

    for j in 1:nj
        ccj = cj[j]

        jb1 = floor(xj[j]/hx)
        jb1 = Int(jb1)

        diff = xj[j] - jb1*hx

        jb1d = min(nspread, jb1)
        jb1u = min(nspread, Mr-jb1-1)
        for k in -nspread:(-jb1d-1)
            istart = jb1 + k + Mr + 1
            zz = ccj * ES_kernel(diff - k*hx, ω, Mr)
            ftau[istart] = ftau[istart] + zz
        end
        for k in -jb1d:jb1u
            istart = jb1 + k + 1
            zz = ccj * ES_kernel(diff - k*hx, ω, Mr)
            ftau[istart] = ftau[istart] + zz
        end
        for k in jb1u+1:nspread
            istart = jb1 + k - Mr + 1
            zz = ccj * ES_kernel(diff - k*hx, ω, Mr)
            ftau[istart] = ftau[istart] + zz
        end
    end

    freq = fftfreq(Mr, 1) .* Mr
    pos = findall(abs.(freq) .<= M/2)

    Ftau = fft(ftau)[pos]

    k = freq[pos]

    Fk = Ftau .* (hx ./ ES_kernel_hat.(k, ω, Mr))

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
# test2 = NUFFTES(cj, xj, M, 10^-12)
#
# # Check accruacy with the L2 norm
# norm(test-test2, 2) / norm(test, 2)
