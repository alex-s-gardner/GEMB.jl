# Synthetic forcing for a maritime temperate glacier site, used by the second synthetic regression
# pin. Every series is a closed-form function of time (no random numbers), so the forcing is
# identical on every platform and Julia version.
#
# Four years of hourly forcing at 61° N and 90 kPa: mean air temperature 270 K with a ±6 K
# seasonal cycle, ±2 K diurnal cycle and ±3 K of weather from three incommensurate periods
# (4.3, 9.7 and 23 d); 3000 kg m-2 yr-1 of precipitation falling only in the warm phase of the
# weather, 43% of it as rain. Seasonal snow over ice, with rain on snow, refreezing, pore water,
# runoff, and vapor exchange at both melting and sub-freezing surfaces. The sun's zenith angle and
# the cloud optical thickness vary, so the albedo's zenith and cloud terms act.

using Dates: DateTime, Hour, Year, value

function synthetic_maritime_forcing(; years::Int=4)
    t0 = DateTime(2001, 1, 1)
    time = collect(t0:Hour(1):(t0 + Year(years) - Hour(1)))
    day = [value(t - t0) / 86_400_000 for t in time]          # days since t0
    hour = [24 * (d - floor(d)) for d in day]

    seasonal = [-cos(2π * (d - 15) / 365.25) for d in day]     # coldest mid-January
    diurnal = [cos(2π * (h - 15) / 24) for h in hour]          # warmest at 15 h
    # Weather: three incommensurate periods (4.3, 9.7 and 23 d), weights summing to 1.
    weather = [0.6 * sin(2π * d / 4.3) + 0.3 * sin(2π * d / 9.7 + 1.0) +
               0.1 * sin(2π * d / 23.0 + 2.0) for d in day]
    storm = [max(w, 0.0) for w in weather]                     # 0 between storms, ≤ 1

    T = [270.0 + 6.0 * s + 2.0 * c + 3.0 * w for (s, c, w) in zip(seasonal, diurnal, weather)]
    # 3000 kg m-2 yr-1, falling only in storms, as the square of storm strength.
    storm_sq_mean = sum(s -> s^2, storm) / length(storm)
    P = [3000.0 / 8766.0 * s^2 / storm_sq_mean for s in storm]
    cloud = [min(1.0, 0.2 + 1.6 * s) for s in storm]

    # Shortwave at 61° N: top-of-atmosphere irradiance from the declination and hour angle,
    # transmissivity 0.75 under clear sky falling to 0.225 under full cloud. The solar zenith
    # angle is that of the same cosine (90° or more at night); the cloud optical thickness and
    # the diffuse share of the shortwave rise linearly from 0 at the clearest sky (cloud 0.2)
    # to 20 and 1 under full cloud.
    lat = deg2rad(61.0)
    μ = map(day, hour) do d, h
        δ = deg2rad(23.44) * sin(2π * (d - 80) / 365.25)
        sin(lat) * sin(δ) + cos(lat) * cos(δ) * cos(deg2rad(15 * (h - 12)))
    end
    sw = [1361.0 * max(m, 0.0) * 0.75 * (1 - 0.7 * c) for (m, c) in zip(μ, cloud)]
    sza = [acosd(clamp(m, -1.0, 1.0)) for m in μ]
    overcast = [(c - 0.2) / 0.8 for c in cloud]
    tau = 20.0 .* overcast
    sw_diffuse = sw .* overcast
    σ = 5.67e-8
    lw = [σ * t^4 * (0.72 + 0.26 * c) for (t, c) in zip(T, cloud)]
    # Vapor pressure from relative humidity 0.80 (clear) to 0.95 (overcast) with respect to
    # water, Bolton (1980) eq. 10.
    vapor = [(0.80 + 0.15 * c) * 611.2 * exp(17.67 * (t - 273.15) / (t - 29.65))
             for (t, c) in zip(T, cloud)]
    wind = [4.0 + 2.5 * s + sin(2π * d / 3.1) for (s, d) in zip(storm, day)]
    pressure = fill(90000.0, length(time))

    return initialize_forcing(time, T, pressure, P, wind, sw, lw, vapor;
        temperature_air_mean=sum(T) / length(T), wind_speed_mean=sum(wind) / length(wind),
        precipitation_mean=sum(P) / years, temperature_observation_height=2.0,
        wind_observation_height=10.0, solar_zenith_angle=sza, cloud_optical_thickness=tau,
        shortwave_downward_diffuse=sw_diffuse)
end
