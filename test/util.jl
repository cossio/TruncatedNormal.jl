using TruncatedNormal, Test, SpecialFunctions
using TruncatedNormal: xerfcx_asym, xerfcx_asym_pi

@test isnan(erf(NaN, NaN))
@test isnan(erf(0, NaN))
@test isnan(erf(NaN, 0))

@test erf(0, 0) == 0
@test erf(+1, +1) == 0
@test erf(-1, -1) == 0

@test erf(0, +Inf) == +1
@test erf(-Inf, 0) == +1
@test erf(+Inf, 0) == -1
@test erf(0, -Inf) == -1

@test erf(-1, 1) ≈ 1.68540158589942973868244127017
@test erf(-2, 2) ≈ 1.99064453003790546832413851273
@test erf(-5, 5) ≈ 1.99999999999692508041114393030

@test xerfcx_asym(500, 3) ≈ 500erfcx(500)
@test xerfcx_asym(1e6, 3) ≈ 1e6erfcx(1e6)
@test xerfcx_asym_pi(500, 3) ≈ 1.12837239688821180302500663804e-6
@test xerfcx_asym_pi(9e5, 3) ≈ 3.48265175028834253964811258438e-13
@test xerfcx_asym_pi(1e6, 3) ≈ 2.82094791773455001286379966421e-13
