# Scan-angle unit test (PLAN.md, Section 3). The converted DR4 files and the installed
# Octofitter must agree on the unit of scan_pos_angle: Octofitter v9 takes degrees, the
# archive's unit (v8 took radians).
#
# The test builds noise-free along-scan abscissae for a single star with a known
# position offset, proper motion, and parallax, computed here independently of
# Octofitter with psi in degrees (eta = dRA* sin psi + dDec cos psi + plx f_AL), on the
# real epochs, scan angles, and parallax factors of the pre-release files. Octofitter's
# log-likelihood at the true parameters must then equal its normalization alone
# (chi^2 = 0). Feeding the same scan angles in radians, or swapping sin and cos, must
# fail, so the test can tell the conventions apart.
#
# Run from the project root:
#     julia --project=. tests/test_scan_angle.jl

using Octofitter, CSV, Test
using Octofitter.TypedTables: Table

const DATA_DIR = joinpath(dirname(@__DIR__), "data_conversion")
const FILES = ["GaiaDR4_data_Gaia4.csv", "GaiaDR4_data_HD114762.csv"]
const TRUTH = (ra0=0.7, dec0=-1.3, pmra=12.0, pmdec=-8.0, plx=20.0, ref_epoch=57388.5)

"""Along-scan abscissae [mas] of the true single-star model, with `psi` in degrees;
`swap=true` uses the wrong projection (cos for RA, sin for Dec)."""
function abscissae(epoch, psi_deg, f_al; swap=false)
    dt = (epoch .- TRUTH.ref_epoch) ./ 365.25
    dra = TRUTH.ra0 .+ TRUTH.pmra .* dt
    ddec = TRUTH.dec0 .+ TRUTH.pmdec .* dt
    s, c = swap ? (cosd.(psi_deg), sind.(psi_deg)) : (sind.(psi_deg), cosd.(psi_deg))
    return dra .* s .+ ddec .* c .+ TRUTH.plx .* f_al
end

"""chi^2 of Octofitter's likelihood at the true parameters: twice the gap between the
normalization and the log-likelihood."""
function chi2(epoch, psi_given, f_al, al, err)
    tab = Table(epoch=epoch, scan_pos_angle=psi_given, parallax_factor_al=f_al,
                centroid_pos_al=al, centroid_pos_error_al=err,
                outlier_flag=zeros(Int, length(epoch)))
    obs = GaiaDR4AstromObs(tab; target=Photocentre, ref=Barycentre, name="GaiaDR4",
        variables=@variables begin
            ra_offset_mas = $(TRUTH.ra0)
            dec_offset_mas = $(TRUTH.dec0)
            pmra = $(TRUTH.pmra)
            pmdec = $(TRUTH.pmdec)
            ref_epoch = $(TRUTH.ref_epoch)
        end)
    ## The host only: the companion is massless and dark, so the photocentre is the
    ## barycentre
    A = Octofitter.Body(name="A", variables=@variables begin
        mass = 1.0
        flux = 1.0
    end)
    b = Octofitter.Body(name="b", about=A, variables=@variables begin
        mass = 1e-12
        flux = 0.0
        a = 1.0; e = 0.0; ω = 0.0; i = 0.5; Ω = 0.0; θ = 0.0
        epoch = 57388.5
    end)
    sys = Octofitter.System(name="scan_angle_test", bodies=[A, b], observations=[obs],
        variables=(@variables begin
            plx = $(TRUTH.plx)
        end), verbosity=0)
    ll = Octofitter.make_ln_like(sys)(sys, Octofitter.make_arr2nt(sys)(Float64[]))
    return 2*(obs.lognorm0 - ll)
end

@testset "scan-angle unit" begin
    @test pkgversion(Octofitter) >= v"9"
    for file in FILES
        d = CSV.File(joinpath(DATA_DIR, file))
        epoch = collect(Float64, d.epoch)
        psi = collect(Float64, d.scan_pos_angle)
        f_al = collect(Float64, d.parallax_factor_al)
        err = collect(Float64, d.centroid_pos_error_al)
        al = abscissae(epoch, psi, f_al)
        @testset "$file" begin
            ## The converter wrote degrees (a radian column never exceeds 2 pi)
            @test maximum(abs, psi) > 2pi
            ## Octofitter reproduces the degree-based abscissae exactly...
            @test chi2(epoch, psi, f_al, al, err) < 1e-6
            ## ...and rejects radians and a swapped projection
            @test chi2(epoch, deg2rad.(psi), f_al, al, err) > 1e3
            @test chi2(epoch, psi, f_al, abscissae(epoch, psi, f_al; swap=true), err) > 1e3
        end
    end
end
