# Synthetic forcing for thin snow over bare ice at low sun, used by the third synthetic regression
# pin. Every series is a closed-form function of time (no random numbers), so the forcing is
# identical on every platform and Julia version.
#
# Ninety days of hourly forcing from 1 April at 72° N and 90 kPa, where the sun stays more than
# 50° from the zenith. Air temperature `temperature_mean` (278 K) with a ±4 K diurnal cycle and
# ±2 K of weather from two incommensurate periods (1.7 and 4.1 d). Snow falls only in the cold,
# night-time phase of the weather, at `snow_rate` (1.5 kg m-2 h-1), and coarsens and melts within
# days, so bare ice is mostly under thin, coarse snow. The cloud optical thickness is never below
# `tau_min`; at 3 or more it sets the albedo's effective zenith cosine to 0.64 (Gardner & Sharp,
# 2010, eq. 10).

using Dates: DateTime, Hour, Day, value

function synthetic_thin_snow_forcing(; days::Int=90, tau_min::Real=3.0,
    temperature_mean::Real=278.0, snow_rate::Real=1.5)
    t0 = DateTime(2001, 4, 1)
    time = collect(t0:Hour(1):(t0 + Day(days) - Hour(1)))
    day = [value(t - t0) / 86_400_000 for t in time]          # days since t0
    hour = [24 * (d - floor(d)) for d in day]
    doy = [90.0 + d for d in day]                              # day of year

    diurnal = [cos(2π * (h - 15) / 24) for h in hour]          # warmest at 15 h
    weather = [0.7 * sin(2π * d / 1.7) + 0.3 * sin(2π * d / 4.1 + 1.0) for d in day]
    T = [temperature_mean + 4.0 * c + 2.0 * w for (c, w) in zip(diurnal, weather)]

    # Snow at `snow_rate` while the weather is in its cold phase at night.
    P = [(w < -0.5 && c < -0.3) ? snow_rate : 0.0 for (w, c) in zip(weather, diurnal)]
    cloud = [0.3 + 0.4 * max(-w, 0.0) for w in weather]       # 0.3 to 0.7

    # Shortwave at 72° N: top-of-atmosphere irradiance from the declination and hour angle,
    # transmissivity 0.75 falling with cloud. The solar zenith angle is that of the same cosine
    # (90° or more at night). The cloud optical thickness rises from `tau_min` with cloud, and
    # the diffuse share of the shortwave with it.
    lat = deg2rad(72.0)
    μ = map(doy, hour) do n, h
        δ = deg2rad(23.44) * sin(2π * (n - 80) / 365.25)
        sin(lat) * sin(δ) + cos(lat) * cos(δ) * cos(deg2rad(15 * (h - 12)))
    end
    sw = [1361.0 * max(m, 0.0) * 0.75 * (1 - 0.7 * c) for (m, c) in zip(μ, cloud)]
    sza = [acosd(clamp(m, -1.0, 1.0)) for m in μ]
    tau = [tau_min + 10.0 * (c - 0.3) for c in cloud]
    sw_diffuse = sw .* [0.2 + c for c in cloud]
    σ = 5.67e-8
    lw = [σ * t^4 * (0.72 + 0.26 * c) for (t, c) in zip(T, cloud)]
    # Vapor pressure from relative humidity 0.80 to 0.92 with respect to water, Bolton (1980)
    # eq. 10.
    vapor = [(0.80 + 0.3 * c) * 611.2 * exp(17.67 * (t - 273.15) / (t - 29.65))
             for (t, c) in zip(T, cloud)]
    wind = [3.0 + 1.5 * abs(w) for w in weather]
    pressure = fill(90000.0, length(time))

    return initialize_forcing(time, T, pressure, P, wind, sw, lw, vapor;
        temperature_air_mean=sum(T) / length(T), wind_speed_mean=sum(wind) / length(wind),
        precipitation_mean=sum(P) * 365.25 / days, temperature_observation_height=2.0,
        wind_observation_height=10.0, solar_zenith_angle=sza, cloud_optical_thickness=tau,
        shortwave_downward_diffuse=sw_diffuse)
end
