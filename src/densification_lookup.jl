"""
Calibrated densification coefficients for the Ligtenberg model.

Each entry maps a region/calibration identifier to a matrix of coefficients:
  [M0_550_offset M0_550_slope M0_830_offset M0_830_slope]
or (for regions with bare-ice calibration):
  [M0_550_offset M0_550_slope M0_830_offset M0_830_slope;
   M1_550_offset M1_550_slope M1_830_offset M1_830_slope]
"""
const DENSIFICATION_COEFFS_M01 = Dict{String,Matrix{Float64}}(
    # --- Antarctic ---
    "Ant_ERA5_GS_SW0"    => [1.5131 0.1317 0.1317 0.2158;
                             1.8422 0.1688 2.4979 0.3225],
    "Ant_ERA5v4_Paolo23" => [2.84 0.32 3.10 0.37],
    "Ant_ERA5_BF_SW1"    => [2.2191 0.2301 2.2917 0.2710],
    "Ant_RACMO_GS_SW0"   => [1.6383 0.1691 1.9991 0.2414],
    "Ant_Ligtenberg"     => [1.435 0.151 2.366 0.293],
    # --- Greenland ---
    "Gre_ERA5_GS_SW0"    => [1.3566 0.1350 1.8705 0.2290;
                             1.4318 0.1055 2.0453 0.2137],
    "Gre_RACMO_GS_SW0"   => [1.2691 0.1184 1.9983 0.2511],
    "Gre_RACMO_GB_SW1"   => [1.7834 0.1409 1.9260 0.1527],
    "Gre_KuipersMunneke"  => [1.042 0.0916 1.734 0.2039],
)

"""
    densification_lookup_M01(densification_coeffs_M01::Symbol)

Return calibrated densification coefficients for the Ligtenberg model.
Matches MATLAB's `densification_lookup_M01.m`.
"""
function densification_lookup_M01(densification_coeffs_M01::Symbol)
    key = String(densification_coeffs_M01)
    haskey(DENSIFICATION_COEFFS_M01, key) ||
        error("Unrecognized densification coefficients: $densification_coeffs_M01")
    return DENSIFICATION_COEFFS_M01[key]
end
