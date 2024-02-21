using Primes
using Test
using DataStructures: SortedDict

import Primes: isprime, primes, primesmask, factor, ismersenneprime, isrieselprime, Factorization

@test primes(10000) == primes(2, 10000) == [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
    73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151,
    157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233,
    239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317,
    331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419,
    421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503,
    509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607,
    613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701,
    709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811,
    821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911,
    919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013,
    1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, 1087, 1091,
    1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153, 1163, 1171, 1181,
    1187, 1193, 1201, 1213, 1217, 1223, 1229, 1231, 1237, 1249, 1259, 1277,
    1279, 1283, 1289, 1291, 1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361,
    1367, 1373, 1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451,
    1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, 1523, 1531,
    1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, 1597, 1601, 1607, 1609,
    1613, 1619, 1621, 1627, 1637, 1657, 1663, 1667, 1669, 1693, 1697, 1699,
    1709, 1721, 1723, 1733, 1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789,
    1801, 1811, 1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889,
    1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987, 1993, 1997,
    1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053, 2063, 2069, 2081, 2083,
    2087, 2089, 2099, 2111, 2113, 2129, 2131, 2137, 2141, 2143, 2153, 2161,
    2179, 2203, 2207, 2213, 2221, 2237, 2239, 2243, 2251, 2267, 2269, 2273,
    2281, 2287, 2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357,
    2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423, 2437, 2441,
    2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531, 2539, 2543, 2549, 2551,
    2557, 2579, 2591, 2593, 2609, 2617, 2621, 2633, 2647, 2657, 2659, 2663,
    2671, 2677, 2683, 2687, 2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729,
    2731, 2741, 2749, 2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819,
    2833, 2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903, 2909, 2917,
    2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999, 3001, 3011, 3019, 3023,
    3037, 3041, 3049, 3061, 3067, 3079, 3083, 3089, 3109, 3119, 3121, 3137,
    3163, 3167, 3169, 3181, 3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251,
    3253, 3257, 3259, 3271, 3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331,
    3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413, 3433, 3449,
    3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511, 3517, 3527, 3529, 3533,
    3539, 3541, 3547, 3557, 3559, 3571, 3581, 3583, 3593, 3607, 3613, 3617,
    3623, 3631, 3637, 3643, 3659, 3671, 3673, 3677, 3691, 3697, 3701, 3709,
    3719, 3727, 3733, 3739, 3761, 3767, 3769, 3779, 3793, 3797, 3803, 3821,
    3823, 3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907, 3911, 3917,
    3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989, 4001, 4003, 4007, 4013,
    4019, 4021, 4027, 4049, 4051, 4057, 4073, 4079, 4091, 4093, 4099, 4111,
    4127, 4129, 4133, 4139, 4153, 4157, 4159, 4177, 4201, 4211, 4217, 4219,
    4229, 4231, 4241, 4243, 4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297,
    4327, 4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409, 4421, 4423,
    4441, 4447, 4451, 4457, 4463, 4481, 4483, 4493, 4507, 4513, 4517, 4519,
    4523, 4547, 4549, 4561, 4567, 4583, 4591, 4597, 4603, 4621, 4637, 4639,
    4643, 4649, 4651, 4657, 4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729,
    4733, 4751, 4759, 4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831,
    4861, 4871, 4877, 4889, 4903, 4909, 4919, 4931, 4933, 4937, 4943, 4951,
    4957, 4967, 4969, 4973, 4987, 4993, 4999, 5003, 5009, 5011, 5021, 5023,
    5039, 5051, 5059, 5077, 5081, 5087, 5099, 5101, 5107, 5113, 5119, 5147,
    5153, 5167, 5171, 5179, 5189, 5197, 5209, 5227, 5231, 5233, 5237, 5261,
    5273, 5279, 5281, 5297, 5303, 5309, 5323, 5333, 5347, 5351, 5381, 5387,
    5393, 5399, 5407, 5413, 5417, 5419, 5431, 5437, 5441, 5443, 5449, 5471,
    5477, 5479, 5483, 5501, 5503, 5507, 5519, 5521, 5527, 5531, 5557, 5563,
    5569, 5573, 5581, 5591, 5623, 5639, 5641, 5647, 5651, 5653, 5657, 5659,
    5669, 5683, 5689, 5693, 5701, 5711, 5717, 5737, 5741, 5743, 5749, 5779,
    5783, 5791, 5801, 5807, 5813, 5821, 5827, 5839, 5843, 5849, 5851, 5857,
    5861, 5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927, 5939, 5953, 5981,
    5987, 6007, 6011, 6029, 6037, 6043, 6047, 6053, 6067, 6073, 6079, 6089,
    6091, 6101, 6113, 6121, 6131, 6133, 6143, 6151, 6163, 6173, 6197, 6199,
    6203, 6211, 6217, 6221, 6229, 6247, 6257, 6263, 6269, 6271, 6277, 6287,
    6299, 6301, 6311, 6317, 6323, 6329, 6337, 6343, 6353, 6359, 6361, 6367,
    6373, 6379, 6389, 6397, 6421, 6427, 6449, 6451, 6469, 6473, 6481, 6491,
    6521, 6529, 6547, 6551, 6553, 6563, 6569, 6571, 6577, 6581, 6599, 6607,
    6619, 6637, 6653, 6659, 6661, 6673, 6679, 6689, 6691, 6701, 6703, 6709,
    6719, 6733, 6737, 6761, 6763, 6779, 6781, 6791, 6793, 6803, 6823, 6827,
    6829, 6833, 6841, 6857, 6863, 6869, 6871, 6883, 6899, 6907, 6911, 6917,
    6947, 6949, 6959, 6961, 6967, 6971, 6977, 6983, 6991, 6997, 7001, 7013,
    7019, 7027, 7039, 7043, 7057, 7069, 7079, 7103, 7109, 7121, 7127, 7129,
    7151, 7159, 7177, 7187, 7193, 7207, 7211, 7213, 7219, 7229, 7237, 7243,
    7247, 7253, 7283, 7297, 7307, 7309, 7321, 7331, 7333, 7349, 7351, 7369,
    7393, 7411, 7417, 7433, 7451, 7457, 7459, 7477, 7481, 7487, 7489, 7499,
    7507, 7517, 7523, 7529, 7537, 7541, 7547, 7549, 7559, 7561, 7573, 7577,
    7583, 7589, 7591, 7603, 7607, 7621, 7639, 7643, 7649, 7669, 7673, 7681,
    7687, 7691, 7699, 7703, 7717, 7723, 7727, 7741, 7753, 7757, 7759, 7789,
    7793, 7817, 7823, 7829, 7841, 7853, 7867, 7873, 7877, 7879, 7883, 7901,
    7907, 7919, 7927, 7933, 7937, 7949, 7951, 7963, 7993, 8009, 8011, 8017,
    8039, 8053, 8059, 8069, 8081, 8087, 8089, 8093, 8101, 8111, 8117, 8123,
    8147, 8161, 8167, 8171, 8179, 8191, 8209, 8219, 8221, 8231, 8233, 8237,
    8243, 8263, 8269, 8273, 8287, 8291, 8293, 8297, 8311, 8317, 8329, 8353,
    8363, 8369, 8377, 8387, 8389, 8419, 8423, 8429, 8431, 8443, 8447, 8461,
    8467, 8501, 8513, 8521, 8527, 8537, 8539, 8543, 8563, 8573, 8581, 8597,
    8599, 8609, 8623, 8627, 8629, 8641, 8647, 8663, 8669, 8677, 8681, 8689,
    8693, 8699, 8707, 8713, 8719, 8731, 8737, 8741, 8747, 8753, 8761, 8779,
    8783, 8803, 8807, 8819, 8821, 8831, 8837, 8839, 8849, 8861, 8863, 8867,
    8887, 8893, 8923, 8929, 8933, 8941, 8951, 8963, 8969, 8971, 8999, 9001,
    9007, 9011, 9013, 9029, 9041, 9043, 9049, 9059, 9067, 9091, 9103, 9109,
    9127, 9133, 9137, 9151, 9157, 9161, 9173, 9181, 9187, 9199, 9203, 9209,
    9221, 9227, 9239, 9241, 9257, 9277, 9281, 9283, 9293, 9311, 9319, 9323,
    9337, 9341, 9343, 9349, 9371, 9377, 9391, 9397, 9403, 9413, 9419, 9421,
    9431, 9433, 9437, 9439, 9461, 9463, 9467, 9473, 9479, 9491, 9497, 9511,
    9521, 9533, 9539, 9547, 9551, 9587, 9601, 9613, 9619, 9623, 9629, 9631,
    9643, 9649, 9661, 9677, 9679, 9689, 9697, 9719, 9721, 9733, 9739, 9743,
    9749, 9767, 9769, 9781, 9787, 9791, 9803, 9811, 9817, 9829, 9833, 9839,
    9851, 9857, 9859, 9871, 9883, 9887, 9901, 9907, 9923, 9929, 9931, 9941,
    9949, 9967, 9973 ]

for n = 100:100:1000
    @test primes(n, 10n) == primes(10n)[(length(primes(n)) + 1):end]
    @test primesmask(n, 10n) == primesmask(10n)[n:end]
end

for T in [Int, BigInt], n = [1:1000;1000000]
    n = convert(T, n)
    f = factor(n)
    @test n == prod(T[p^k for (p, k) = f])
    prime = n != 1 && length(f) == 1 && get(f, n, 0) == 1
    @test isprime(n) == prime || n

    s = primesmask(n)
    for k = 1:n
        @test s[k] == isprime(k)
        @test s[k] == primesmask(k, k)[1]
    end
end

# Issue 25
@test primes(-2, 7) == primes(0, 7) == primes(2, 7)

@test !isprime(1000000003)
@test !isprime(1000000005)
@test  isprime(1000000007)
@test  isprime(1000000009)
@test !isprime(1000000011)
@test !isprime(1000000013)

@test !isprime(10000000015)
@test !isprime(10000000017)
@test  isprime(10000000019)
@test !isprime(10000000021)
@test !isprime(10000000023)

@test !isprime(9223372036854775779)
@test !isprime(9223372036854775781)
@test  isprime(9223372036854775783)
@test !isprime(9223372036854775785)
@test !isprime(9223372036854775787)

@test !isprime(0xffffffffffffffc1)
@test !isprime(0xffffffffffffffc3)
@test  isprime(0xffffffffffffffc5)
@test !isprime(0xffffffffffffffc7)
@test !isprime(0xffffffffffffffc9)

for T in [Int8, UInt8, Int16, UInt16, Int128, UInt128]
    @test isprime(T(2))
    @test !isprime(T(4))
end


# Base issue #5210
@test prod([k^v for (k, v) in factor(typemax(UInt32))]) == typemax(UInt32)
@test prod([k^v for (k, v) in factor(typemax(Int8))]) == typemax(Int8)

# factorization of factors > 2^16
@test factor((big(2)^31 - 1)^2) == Dict(big(2^31 - 1) => 2)
@test factor((big(2)^31 - 1) * (big(2)^17 - 1)) == Dict(big(2^31 - 1) => 1, big(2^17 - 1) => 1)

# fast factorization of Int128 (Base issue #11477)
@test factor((Int128(2)^39 - 7)^2) == Dict(Int128(2)^39 - 7 => 2)

# Base issue #9611
@test factor(Int128(2)^101 + 1) == Dict(3 => 1, 845100400152152934331135470251 => 1)

# test second branch, after all small primes in list have been searched
@test factor(10009 * Int128(1000000000000037)) == Dict(10009 => 1, 1000000000000037 => 1)

for n = 1:100
    m = 1
    for (p,k) in factor(n)
        m *= p^k
    end
    @test n == m
end

let i = rand(1:2^(3 * min(Sys.WORD_SIZE,64) ÷ 4))
    @test primes(i, i + 300) == filter(isprime, i:(i + 300)) == filter(isprime, big(i):big(i + 300))
end


@test isprime(BigInt(1000000007))
@test isprime(BigInt(1000000007), 1)
@test isprime(BigInt(10000000019))
@test isprime(parse(BigInt,"359334085968622831041960188598043661065388726959079837"))
@test !isprime(BigInt(1))
@test !isprime(BigInt(10000000020))

# project euler 3: 6857
@test maximum(keys(factor(600851475143))) == 6857

# project euler 7: 104743
euler7(n) = primes(floor(Int, n * log(n * log(n))))[n]
@test euler7(10001) == 104743

# project euler 10: 142913828922
@test sum(map(Int64,primes(2000000))) == 142913828922

# factor(Vector, n)
for V in (Vector, Vector{Int}, Vector{Int128})
    @test factor(V, 1) == Int[]
    @test factor(V, 3) == [3]
    @test factor(V, 4) == [2, 2]
end

# factor with non-default associative containers
@test factor(SortedDict, 100) == factor(Dict, 100) == factor(100)

# factor sets
@test factor(Set, 100) == Set([2, 5])
@test factor(BitSet, 100) == BitSet([2, 5])

# factor other things and fail
@test_throws MethodError factor(Int, 10)
@test_throws MethodError factor(Any, 10)
@test_throws MethodError factor(Tuple, 10)

# factor non-positive numbers:
@test factor(0)  == Dict( 0 => 1)
@test factor(-1) == Dict(-1 => 1)
@test factor(-9) == Dict(-1 => 1, 3 => 2)

@test factor(typemin(Int32))  == Dict(-1 => 1, 2 => 31)
@test factor(typemin(Int64))  == Dict(-1 => 1, 2 => 63)
@test factor(typemin(Int128)) == Dict(-1 => 1, 2 => 127)

@test factor(1) == Dict{Int,Int}()

# correctly return sign of factored numbers
@test sign(factor(100)) == 1
@test sign(factor(1)) == 1
@test sign(factor(0)) == 0
@test sign(factor(-1)) == -1
@test sign(factor(-100)) == -1

# factor returns a sorted dict
@test all([issorted(collect(factor(rand(Int)))) for x in 1:100])

# test eachfactor iteration
for T in (Int32, Int64, BigInt)
    @test iterate(eachfactor(T(36))) == ((T(2), 2), T.((9, 3)))
    @test iterate(eachfactor(T(7^2*5^3))) == ((T(5), 3), T.((49, 5)))
    @test iterate(eachfactor(T(257))) == ((T(257), 1), T.((1, 257)))
    @test iterate(eachfactor(T(nextprime(2^16)))) == ((T(65537), 1), T.((1, 65537)))
    for (p,e) in eachfactor(T(901800900))
        @test (p,e) isa Tuple{T, Int}
    end
end

# Lucas-Lehmer
@test !ismersenneprime(2047)
@test ismersenneprime(8191)
@test_throws ArgumentError ismersenneprime(9)
@test_throws ArgumentError ismersenneprime(9, check=true)
# test the following does not throw
ismersenneprime(9, check=false)

#  Lucas-Lehmer-Riesel
@test_throws ArgumentError isrieselprime(1000, 511)
@test_throws ArgumentError isrieselprime(0, 1)
@test isrieselprime(1, 8191) == ismersenneprime(8191)  # Case 1
@test isrieselprime(3, BigInt(2)^607 - 1)              # Case 2
@test_throws ErrorException isrieselprime(20, 31)      # Case `else`

# @testset "Factorization{$T} as an AbstractDict" for T = (Int, UInt, BigInt)
for T = (Int, UInt, BigInt)
    d = Dict(map(Pair, rand(T(1):T(100), 30), 1:30))
    f = Factorization{T}(d)
    @test f == d == Dict(f) == Factorization(d) == convert(Factorization, d)
    @test collect(f) == sort!(collect(d)) # test start/next/done
    @test length(f) == length(d)
    @test get(f, T(101), nothing) === nothing
    @test f[101] == 0
    @test f[0] == 0
    f[0] = 1
    @test get(f, T(0), 0) == 1
    @test f[0] == 1
    @test f[101] == 0
    f = Factorization{T}()
    for i in 100:-1:1
        f[T(i)] = i
    end
    @test length(f) == 100
    @test issorted(f.pe)
end

d = Dict(:a=>1,:b=>2)
f = Factorization(d)
@test f == d == Dict(f) == convert(Factorization, d)
@test collect(f) == sort!(collect(d)) # test start/next/done
@test length(f) == length(d)
@test get(f, :c, nothing) === nothing
@test f[:c] == 0
@test f[:a] == 1
f[:c] = 1
@test get(f, :c, 0) == 1

# dumb implementation of Euler totient for correctness tests
ϕ(n) = count(m -> gcd(n, m) == 1, 1:n)

@testset "totient(::$T)" for T in [Int16, Int32, Int64, BigInt]
    for n in 1:1000
        @test ϕ(T(n)) == totient(T(n)) == totient(-T(n))
    end
    @test_throws ArgumentError @inferred(totient(T(0)))
end

# check some big values
@testset "totient() correctness" begin
    @test totient(2^4 * 3^4 * 5^4) == 216000
    @test totient(big"2"^1000) == big"2"^999

    some_coprime_numbers = BigInt[
        450000000, 1099427429702334733, 200252151854851, 1416976291499, 7504637909,
        1368701327204614490999, 662333585807659, 340557446329, 1009091
    ]

    for i in some_coprime_numbers
        for j in some_coprime_numbers
            if i ≠ j
                @test totient(i*j) == totient(i) * totient(j)
            end
        end
        # can use directly with Factorization
        @test totient(i) == totient(factor(i))
    end
end

# brute-force way to get divisors. Same elements as divisors(n), but order may differ.
divisors_brute_force(n) = [d for d in one(n):n if iszero(n % d)]

@testset "divisors(::$T)" for T in [Int16, Int32, Int64, BigInt]
    # 1 and 0 are handled specially
    @test divisors(one(T)) == divisors(-one(T)) == T[one(T)]
    @test divisors(zero(T)) == T[]

    for n in 2:1000
        ds = divisors(T(n))
        @test ds == divisors(-T(n))
        @test sort!(ds) == divisors_brute_force(T(n))
    end
end

@testset "divisors(::Factorization)" begin
    # divisors(n) calls divisors(factor(abs(n))), so the previous testset covers most cases.
    # We just need to verify that the following cases are also handled correctly:
    @test divisors(factor(1)) == divisors(factor(-1)) == [1] # factorizations of 1 and -1
    @test divisors(factor(-56)) == divisors(factor(56)) == [1, 2, 4, 8, 7, 14, 28, 56] # factorizations of negative numbers
    @test isempty(divisors(factor(0)))
end

# check copy property for big primes relied upon in nextprime/prevprime
for n = rand(big(-10):big(10), 10)
    @test n+0 !== n
end

@testset "nextprime(::$T)" for T in (Int64, Int32, BigInt)
    for N in zip(T[rand(-1000:1), 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
                   2^20, 2^30],
                 T[2, 2, 3, 5, 5, 7, 7, 11, 11, 11, 11, 13,
                   1048583, 1073741827],
                 T[3, 3, 5, 7, 7, 11, 11, 13, 13, 13, 13, 17,
                   1048589, 1073741831])
        @test nextprime(N[1]) == N[2]
        @test nextprime(N[1], 1) == N[2]
        @test nextprime(N[1], 2) == N[3]
    end
    @test nextprime(Int64(2)^60) == 1152921504606847009
    @test nextprime(Int64(2)^60, 2) == 1152921504606847067
    @test_throws DomainError nextprime(rand(-100:100), 0)
    for i = rand(-100:-1, 100),
        n = rand(600:2^20)
        @test nextprime(n, i) == prevprime(n, -i)
    end

    # interval
    @test nextprime(0, interval=2) == 2
    @test nextprime(0, interval=3) == 3
    @test_throws DomainError nextprime(0, 2, interval=2)
    @test_throws DomainError nextprime(0, interval=4)
    @test nextprime(-20, interval=2) == 2
    @test nextprime(4, interval=5) == 19
    @test nextprime(4, 2, interval=5) == 29
    @test nextprime(4, 3, interval=5) == 59
    @test gcd(nextprime(2^17-1024+1; interval=1024) - 1, 1024) == 1024
    @test gcd(nextprime(1024*rand(1:2^6)+1; interval=1024) - 1, 1024) == 1024
end

@testset "prevprime(::$T)" for T in (Int64, Int32, BigInt)
    for N in zip(T[4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,
                   2^20, 2^30],
                 T[3, 5, 5, 7, 7, 7, 7, 11, 11, 13, 13,
                   1048573, 1073741789],
                 T[2, 3, 3, 5, 5, 5, 5, 7, 7, 11, 11,
                   1048571, 1073741783])
        @test prevprime(N[1]) == N[2]
        @test prevprime(N[1], 1) == N[2]
        @test prevprime(N[1], 2) == N[3]
    end
    @test prevprime(Int64(2)^60) == 1152921504606846883
    @test prevprime(Int64(2)^60, 2) == 1152921504606846869
    @test prevprime(3) == 3
    @test prevprime(3, 2) == 2
    @test prevprime(2) == 2
    @test prevprime(typemax(UInt8)) == prevprime(big(typemax(UInt8)))
    @test typeof(prevprime(typemax(UInt8))) === UInt8
    @test_throws ArgumentError prevprime(2, 2)

    @test_throws DomainError prevprime(rand(-100:100), 0)
    for n = rand(-100:1, 10)
        @test_throws ArgumentError prevprime(n)
        @test_throws ArgumentError prevprime(n, 1)
        @test_throws ArgumentError prevprime(n, 2)
    end

    for i = rand(-100:-1, 100),
        n = rand(600:2^20)
        @test prevprime(n, i) == nextprime(n, -i)
    end

    # interval
    @test gcd(prevprime(2^17-1024+1; interval=1024) - 1, 1024) == 1024
    @test gcd(prevprime(1024*rand(12:2^6)+1; interval=1024) - 1, 1024) == 1024
    @test prevprime(10, interval=4) == 2
    @test [prevprime(11, i, interval=2) for i=1:4] == [11, 7, 5, 3]
    @test [prevprime(11, i, interval=3) for i=1:3] == [11, 5, 2]
end

@testset "prime(::$T)" for T = (Int64, Int32, BigInt)
    for (n, p) = zip([1, 2, 3, 4, 5, 6, 7, 100, 1000, 10000],
                     T[2, 3, 5, 7, 11, 13, 17, 541, 7919, 104729])
        @test prime(n) == p
        @test prime(T, n) == p
    end
end

for T in (Int, UInt, BigInt)
    for n in [T(1); rand(T(2):T(100000), 10)]
        # for n=T(1), must not error out (#51)
        for C = (Factorization, Vector, Dict)
            @test prodfactors(factor(C, n)) == n
        end
        if Primes.radical(n) == n
            for C = (Set, BitSet)
                @test prodfactors(factor(C, n)) == n
            end
        end
    end
    @test prodfactors(factor(Set, T(123456))) == 3858
    @test prod(factor(T(123456))) == 123456
end

@testset "nextprimes(::$T)" for T = (Int32, Int64, BigInt)
    for (i, p) in enumerate(nextprimes(T))
        @test nextprime(0, i) == p
        i > 20 && break
    end
    @test nextprimes() == nextprimes(Int)
    for (i, p) in enumerate(nextprimes(T(5)))
        @test nextprime(T(5), i) == p
        i > 20 && break
    end
    @test nextprimes(T(5), 10) == [nextprime(T(5), i) for i=1:10]
    @test nextprimes(1, 1)[1] == nextprimes(2, 1)[1] == 2
    @test nextprimes(3, 1)[1] == 3
    @test nextprimes(4, 1)[1] == nextprimes(5, 1)[1] == 5
    @test eltype(nextprimes(10)) == Int
    @test eltype(nextprimes(big(10))) == BigInt
    @test Base.IteratorEltype(nextprimes(10)) == Base.HasEltype()
    @test Base.IteratorSize(nextprimes(10)) == Base.IsInfinite()

end


@testset "prevprimes(::$T)" for T = (Int32, Int64, BigInt)
    for (i, p) in enumerate(prevprimes(T(500)))
        @test prevprime(T(500), i) == p
        i > 20 && break
    end
    @test prevprimes(T(500), 10) == [prevprime(T(500), i) for i=1:10]
    @test prevprimes(6, 1)[1] == prevprimes(5, 1)[1] == 5
    @test prevprimes(4, 1)[1] == prevprimes(3, 1)[1] == 3
    @test prevprimes(2, 1)[1] == 2
    @test isempty(prevprimes(1, 1))
    let p8 = collect(prevprimes(typemax(Int8)))
        @test length(p8) == 31
        @test p8[end] == 2
        @test p8[1] == 127
        @test eltype(p8) == Int8
    end
    @test eltype(prevprimes(10)) == Int
    @test eltype(prevprimes(big(10))) == BigInt
    @test Base.IteratorEltype(prevprimes(10)) == Base.HasEltype()
    @test Base.IteratorSize(prevprimes(10)) == Base.SizeUnknown()
end

@testset "primes with huge arguments" begin
    if Base.Sys.WORD_SIZE == 64
        @test primes(2^63-200, 2^63-1) == [9223372036854775643, 9223372036854775783]
    end
    @test primes(2^31-20, 2^31-1) == [2147483629, 2147483647]
end
