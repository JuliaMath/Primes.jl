# Sieve tests, split out for a fast standalone dev loop:
#   julia --project -e 'using TestEnv; TestEnv.activate(); include("test/sieve_tests.jl")'
# Also included from runtests.jl so CI runs them too. Self-contained imports make
# both entry points work (re-importing when run from runtests.jl is harmless).
using Primes
using Test
import Primes: isprime, primes, primesmask

@testset "segmented sieve" begin
    # brute-force reference independent of the sieve
    bruteprimes(lo, hi) = [n for n in max(2, lo):hi if isprime(n)]
    brutemask(lo, hi)   = Bool[isprime(n) for n in lo:hi]

    # small ranges, lo==7 boundary, sub-49 ranges
    for (lo, hi) in [(1, 100), (7, 100), (2, 7), (7, 7), (10, 48), (7, 49), (1, 1000)]
        @test primes(lo, hi) == bruteprimes(lo, hi)
        @test primesmask(lo, hi) == brutemask(lo, hi)
    end

    # ranges that span multiple segment windows (a window covers ~10^6 integers)
    for (lo, hi) in [(1, 3_000_000), (2_000_000, 3_000_000), (999_983, 2_500_000)]
        @test primes(lo, hi) == bruteprimes(lo, hi)
    end
    @test primesmask(2_000_000, 2_100_000) == brutemask(2_000_000, 2_100_000)

    # consistency: primes(n) and primesmask(n) agree
    let pm = primesmask(2_000_000)
        @test primes(2_000_000) == [n for n in 2:2_000_000 if pm[n]]
    end

    # direct core tests: the SegmentedSieve iterator (yields itself per filled window) and
    # the each_lane_prime consumer. The sieve masks values < max(7, lo) and > hi.
    function collect_sieve(lo, hi)
        out = Int[]
        for s in Primes.SegmentedSieve(max(7, lo), hi)
            Primes.each_lane_prime(s) do p
                push!(out, p)
            end
        end
        out
    end
    @test collect_sieve(7, 1000) == [n for n in 7:1000 if isprime(n)]
    # a range crossing several segment windows
    @test collect_sieve(1_000_000, 1_500_000) == [n for n in 1_000_000:1_500_000 if isprime(n)]
    # a window starting mid-chunk (lo not 64-block aligned)
    @test collect_sieve(999_983, 1_100_000) == [n for n in 999_983:1_100_000 if isprime(n)]

    # eachprime streams the same primes as primes(), and is lazily takeable
    @test eltype(eachprime(100)) == Int
    @test collect(eachprime(100)) == primes(100)
    for (lo, hi) in [(1, 100), (7, 100), (2, 7), (10, 48), (8, 10), (1, 3_000_000), (999_983, 1_100_000)]
        @test collect(eachprime(lo, hi)) == primes(lo, hi)
    end
    @test collect(eachprime(6)) == [2, 3, 5]
    @test isempty(collect(eachprime(8, 10)))
    @test collect(Iterators.take(eachprime(10^6), 6)) == [2, 3, 5, 7, 11, 13]
    @test_throws ArgumentError eachprime(10, 5)
    # NOTE: the widemul / `q < 0` overflow guards only trigger for hi near
    # typemax(Int), which needs base primes up to ~isqrt(typemax) ≈ 3e9 — far too
    # expensive to sieve in a test. Those guards are carried over verbatim from the
    # existing implementation and are not unit-tested here.
end
