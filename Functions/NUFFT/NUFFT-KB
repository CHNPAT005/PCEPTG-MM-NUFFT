## Author: Patrick Chang
# Script file for the Nonuniform Fast Fourier Transform

# Kaiser-Bessel - spreading to grid points near the source points

using ArgCheck; using FFTW; using SpecialFunctions

# Supporting functions
#---------------------------------------------------------------------------

function kaiser_bessel(x,n,m,sigma)
    b = pi*(2-1/sigma)
    arg = m^2-n^2*x^2
    if abs(x) < m/n
        y = sinh(b*sqrt(arg))/sqrt(arg)/pi
    elseif abs(x) > m/n
        y = zero(x)
    else
        y = b/pi
    end
    return y
end

function kaiser_bessel_hat(k,n,m,sigma)
    b = pi*(2-1/sigma)
    return besseli(0,m*sqrt(b^2-(2*pi*k/n)^2))
end

function kaiser_bessel_conv(x,n,m,sigma)
    store = 0
    # for r in -m:m
    for r in -1:1   # this actually works, only need it to get the end points working!
        store = store + kaiser_bessel(x+r,n,m,sigma)
    end
    return store
end

#---------------------------------------------------------------------------

# cj = the source strength
# xj = the corresponding time for cj. xj ∈ [-0.5, 0.5]
# M = number of Fourier modes to return:
#   M even returns [-M/2, (M-1)/2]; M odd returns [-M/2, M/2]
# tol = tolerance level for error - controls the spreading

function NUFFTKB(cj, xj, M, tol; kwargs...)

    @argcheck size(cj) == size(xj) DimensionMismatch

    nj = size(cj)[1]

    R = 2
    Mr = R*M
    nspread = Int(floor((ceil(log(10, 1/tol)) + 2)/2))
    #nspread = Int(floor((ceil(log(10, 1/tol)) + 1)/2))
    # hx = 2*pi/Mr
    # EpsGrid = collect(0:1:Mr-1) ./ Mr .- 0.5

    ftau = zeros(ComplexF64, Mr, 1)

    for j in 1:nj
        ccj = cj[j]

        jb1 = floor((xj[j])*Mr)
        jb1 = Int(jb1)

        diff = xj[j] - jb1/Mr

        jb1d = min(nspread, jb1)
        jb1u = min(nspread, Mr-jb1-1)
        for k in -nspread:(-jb1d-1)
            istart = jb1 + k + Mr + 1
            # zz = ccj * kaiser_bessel_conv(xj[j]-EpsGrid[istart], Mr, nspread, R)
            zz = ccj * kaiser_bessel(diff - k/Mr, Mr, nspread, R)
            ftau[istart] = ftau[istart] + zz
        end
        for k in -jb1d:jb1u
            istart = jb1 + k + 1
            # zz = ccj * kaiser_bessel_conv(xj[j]-EpsGrid[istart], Mr, nspread, R)
            zz = ccj * kaiser_bessel(diff - k/Mr, Mr, nspread, R)
            ftau[istart] = ftau[istart] + zz
        end
        for k in jb1u+1:nspread
            istart = jb1 + k - Mr + 1
            # zz = ccj * kaiser_bessel_conv(xj[j]-EpsGrid[istart], Mr, nspread, R)
            zz = ccj * kaiser_bessel(diff - k/Mr, Mr, nspread, R)
            ftau[istart] = ftau[istart] + zz
        end
    end

    freq = fftfreq(Mr, 1) .* Mr
    pos = findall(abs.(freq) .<= M/2)

    Ftau = fft(ftau)[pos]

    k = freq[pos]

    Fk = Ftau ./ kaiser_bessel_hat.(k, Mr, nspread, R)

    return Fk
end

# # Simple test case
# 
# nj = 10
# x = (collect(0:nj-1) + 0.5 .* rand(nj))
# xj = (x .- minimum(x)) .* (2*pi / (maximum(x) - minimum(x)))
# cj = rand(nj) + 0im*rand(nj)
# xjj = xj ./ (2*pi)
#
# M = 11
# Mr = 2*M
# freq = fftfreq(Mr, 1) .* Mr
# pos = findall(abs.(freq) .<= M/2)
# k = freq[pos]
#
# test = (cj' * exp.(1im .* xj * k'))'
# test2 = NUFFTKB(cj, xjj, M, 10^-12)
#
# test - test2
