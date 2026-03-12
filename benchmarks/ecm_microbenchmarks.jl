# ECM Micro-Benchmarks
#
# Benchmarks individual ECM building blocks:
#   - _mulmod / _mulmod! (modular multiplication)
#   - _ecm_add / _ecm_add! (differential point addition)
#   - _ecm_double / _ecm_double! (point doubling)
#   - _ecm_scalar_mul / _ecm_scalar_mul! (Montgomery ladder)
#
# Each operation is tested in four modes:
#   1. BigInt GMP in-place (specialised path)
#   2. BigInt generic (allocating, functional)
#   3. UInt128 generic (note: widen(UInt128) == BigInt, so still allocates)
#   4. UInt256 generic via BitIntegers.jl (widen stays in LLVM fixed-width)
#
# Usage:
#   julia --project=benchmarks benchmarks/ecm_microbenchmarks.jl

using Primes
using BenchmarkTools
using BitIntegers

# ---------------------------------------------------------------------------
# Setup: representative values for a ~40-digit modulus
# ---------------------------------------------------------------------------
const N_BIG = big"824633720831" * big"1000000007"
const A24_BIG = mod(big"123456789012345", N_BIG)
const PX_BIG  = mod(big"314159265358979323846", N_BIG)
const PZ_BIG  = mod(big"271828182845904523536", N_BIG)
const QX_BIG  = mod(big"141421356237309504880", N_BIG)
const QZ_BIG  = mod(big"173205080756887729352", N_BIG)
const DX_BIG  = mod(big"223606797749978969640", N_BIG)
const DZ_BIG  = mod(big"264575131106459059050", N_BIG)

# UInt128 versions (note: widen(UInt128) == BigInt in base Julia)
const N_U128   = UInt128(824633720831) * UInt128(1000000007)
const A24_U128 = UInt128(mod(big"123456789012345", big(N_U128)))
const PX_U128  = UInt128(mod(big"314159265358979323846", big(N_U128)))
const PZ_U128  = UInt128(mod(big"271828182845904523536", big(N_U128)))
const QX_U128  = UInt128(mod(big"141421356237309504880", big(N_U128)))
const QZ_U128  = UInt128(mod(big"173205080756887729352", big(N_U128)))
const DX_U128  = UInt128(mod(big"223606797749978969640", big(N_U128)))
const DZ_U128  = UInt128(mod(big"264575131106459059050", big(N_U128)))

# UInt256 versions (widen(UInt256) == UInt512, stays in LLVM fixed-width)
const N_U256   = UInt256(824633720831) * UInt256(1000000007)
const A24_U256 = UInt256(mod(big"123456789012345", big(N_U256)))
const PX_U256  = UInt256(mod(big"314159265358979323846", big(N_U256)))
const PZ_U256  = UInt256(mod(big"271828182845904523536", big(N_U256)))
const QX_U256  = UInt256(mod(big"141421356237309504880", big(N_U256)))
const QZ_U256  = UInt256(mod(big"173205080756887729352", big(N_U256)))
const DX_U256  = UInt256(mod(big"223606797749978969640", big(N_U256)))
const DZ_U256  = UInt256(mod(big"264575131106459059050", big(N_U256)))

# GMP buffers
const BUF = Primes.ECMBuffers()
const TMP_MUL = BigInt()

println("=" ^ 70)
println("  ECM Micro-Benchmarks")
println("=" ^ 70)
println()
println("  Note: widen(UInt128) == BigInt   → still allocates (no LLVM benefit)")
println("        widen(UInt256) == UInt512  → pure LLVM fixed-width arithmetic")

# ---------------------------------------------------------------------------
# 1. Modular multiplication
# ---------------------------------------------------------------------------
println("\n--- _mulmod ---")

print("  BigInt GMP in-place (_mulmod!):   ")
dst = BigInt()
b1 = @benchmark Primes._mulmod!($dst, $PX_BIG, $PZ_BIG, $N_BIG, $TMP_MUL)
println("$(round(median(b1).time; digits=1)) ns,  $(median(b1).memory) bytes")

print("  BigInt generic    (_mulmod):      ")
b2 = @benchmark Primes._mulmod($PX_BIG, $PZ_BIG, $N_BIG)
println("$(round(median(b2).time; digits=1)) ns,  $(median(b2).memory) bytes")

print("  UInt128 generic   (_mulmod):      ")
b3 = @benchmark Primes._mulmod($PX_U128, $PZ_U128, $N_U128)
println("$(round(median(b3).time; digits=1)) ns,  $(median(b3).memory) bytes")

print("  UInt256 generic   (_mulmod):      ")
b3b = @benchmark Primes._mulmod($PX_U256, $PZ_U256, $N_U256)
println("$(round(median(b3b).time; digits=1)) ns,  $(median(b3b).memory) bytes")

# ---------------------------------------------------------------------------
# 2. Differential point addition
# ---------------------------------------------------------------------------
println("\n--- _ecm_add ---")

res_X, res_Z = BigInt(), BigInt()

print("  BigInt GMP in-place (_ecm_add!):  ")
b4 = @benchmark Primes._ecm_add!($res_X, $res_Z,
    $PX_BIG, $PZ_BIG, $QX_BIG, $QZ_BIG, $DX_BIG, $DZ_BIG, $N_BIG, $BUF)
println("$(round(median(b4).time; digits=1)) ns,  $(median(b4).memory) bytes")

print("  BigInt generic    (_ecm_add):     ")
b5 = @benchmark Primes._ecm_add($PX_BIG, $PZ_BIG, $QX_BIG, $QZ_BIG, $DX_BIG, $DZ_BIG, $N_BIG)
println("$(round(median(b5).time; digits=1)) ns,  $(median(b5).memory) bytes")

print("  UInt128 generic   (_ecm_add):     ")
b6 = @benchmark Primes._ecm_add($PX_U128, $PZ_U128, $QX_U128, $QZ_U128, $DX_U128, $DZ_U128, $N_U128)
println("$(round(median(b6).time; digits=1)) ns,  $(median(b6).memory) bytes")

print("  UInt256 generic   (_ecm_add):     ")
b6b = @benchmark Primes._ecm_add($PX_U256, $PZ_U256, $QX_U256, $QZ_U256, $DX_U256, $DZ_U256, $N_U256)
println("$(round(median(b6b).time; digits=1)) ns,  $(median(b6b).memory) bytes")

# ---------------------------------------------------------------------------
# 3. Point doubling
# ---------------------------------------------------------------------------
println("\n--- _ecm_double ---")

print("  BigInt GMP in-place (_ecm_double!): ")
b7 = @benchmark Primes._ecm_double!($res_X, $res_Z, $PX_BIG, $PZ_BIG, $N_BIG, $A24_BIG, $BUF)
println("$(round(median(b7).time; digits=1)) ns,  $(median(b7).memory) bytes")

print("  BigInt generic    (_ecm_double):    ")
b8 = @benchmark Primes._ecm_double($PX_BIG, $PZ_BIG, $N_BIG, $A24_BIG)
println("$(round(median(b8).time; digits=1)) ns,  $(median(b8).memory) bytes")

print("  UInt128 generic   (_ecm_double):    ")
b9 = @benchmark Primes._ecm_double($PX_U128, $PZ_U128, $N_U128, $A24_U128)
println("$(round(median(b9).time; digits=1)) ns,  $(median(b9).memory) bytes")

print("  UInt256 generic   (_ecm_double):    ")
b9b = @benchmark Primes._ecm_double($PX_U256, $PZ_U256, $N_U256, $A24_U256)
println("$(round(median(b9b).time; digits=1)) ns,  $(median(b9b).memory) bytes")

# ---------------------------------------------------------------------------
# 4. Montgomery ladder scalar multiplication
# ---------------------------------------------------------------------------
println("\n--- _ecm_scalar_mul (k=997, prime) ---")

const K = 997

print("  BigInt GMP in-place (_ecm_scalar_mul!): ")
b10 = @benchmark Primes._ecm_scalar_mul!($res_X, $res_Z, $K, $PX_BIG, $PZ_BIG, $N_BIG, $A24_BIG, $BUF)
println("$(round(median(b10).time / 1e3; digits=2)) μs,  $(median(b10).memory) bytes")

print("  BigInt generic    (_ecm_scalar_mul):    ")
b11 = @benchmark Primes._ecm_scalar_mul($K, $PX_BIG, $PZ_BIG, $N_BIG, $A24_BIG)
println("$(round(median(b11).time / 1e3; digits=2)) μs,  $(median(b11).memory) bytes")

print("  UInt128 generic   (_ecm_scalar_mul):    ")
b12 = @benchmark Primes._ecm_scalar_mul($K, $PX_U128, $PZ_U128, $N_U128, $A24_U128)
println("$(round(median(b12).time / 1e3; digits=2)) μs,  $(median(b12).memory) bytes")

print("  UInt256 generic   (_ecm_scalar_mul):    ")
b12b = @benchmark Primes._ecm_scalar_mul($K, $PX_U256, $PZ_U256, $N_U256, $A24_U256)
println("$(round(median(b12b).time / 1e3; digits=2)) μs,  $(median(b12b).memory) bytes")

# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------
println("\n" * "=" ^ 70)
println("  Summary: median time")
println("  (GMP in-place → Generic BigInt → UInt128 → UInt256)")
println("=" ^ 70)

function fmt_speedup(base, other)
    s = base / other
    "$(round(s; digits=2))×"
end

println("  _mulmod:         $(round(median(b1).time; digits=0)) ns → $(round(median(b2).time; digits=0)) ns ($(fmt_speedup(median(b1).time, median(b2).time))) → $(round(median(b3).time; digits=0)) ns ($(fmt_speedup(median(b1).time, median(b3).time))) → $(round(median(b3b).time; digits=0)) ns ($(fmt_speedup(median(b1).time, median(b3b).time)))")
println("  _ecm_add:        $(round(median(b4).time; digits=0)) ns → $(round(median(b5).time; digits=0)) ns ($(fmt_speedup(median(b4).time, median(b5).time))) → $(round(median(b6).time; digits=0)) ns ($(fmt_speedup(median(b4).time, median(b6).time))) → $(round(median(b6b).time; digits=0)) ns ($(fmt_speedup(median(b4).time, median(b6b).time)))")
println("  _ecm_double:     $(round(median(b7).time; digits=0)) ns → $(round(median(b8).time; digits=0)) ns ($(fmt_speedup(median(b7).time, median(b8).time))) → $(round(median(b9).time; digits=0)) ns ($(fmt_speedup(median(b7).time, median(b9).time))) → $(round(median(b9b).time; digits=0)) ns ($(fmt_speedup(median(b7).time, median(b9b).time)))")
println("  _ecm_scalar_mul: $(round(median(b10).time/1e3; digits=1)) μs → $(round(median(b11).time/1e3; digits=1)) μs ($(fmt_speedup(median(b10).time, median(b11).time))) → $(round(median(b12).time/1e3; digits=1)) μs ($(fmt_speedup(median(b10).time, median(b12).time))) → $(round(median(b12b).time/1e3; digits=1)) μs ($(fmt_speedup(median(b10).time, median(b12b).time)))")
println("=" ^ 70)
