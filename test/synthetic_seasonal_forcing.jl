# Synthetic forcing for seasonal snow on soil, used by the `lower_boundary = :ground` tests. Every
# series is a closed-form function of time (no random numbers).
#
# Hourly forcing from 1 October at 40° N and 70 kPa (about 3000 m): air temperature
# `temperature_mean` (273.5 K) with a ±`amplitude` (10 K) seasonal cycle coldest in mid-January,
# ±4 K diurnal and ±3 K of weather from two incommensurate periods (5.3 and 11.9 d). 700 kg m-2
# yr-1 of precipitation falling only in storms, as snow when the air is at or below 0 °C. Snow
# builds from November, peaks in spring and melts out by early summer, leaving the ground bare
# into autumn.

using Dates: DateTime, Hour, Year, value

function synthetic_seasonal_forcing(; years::Int=2, temperature_mean::Real=273.5,
    amplitude::Real=10.0, precipitation_annual::Real=700.0)
    t0 = DateTime(2001, 10, 1)
    time = collect(t0:Hour(1):(t0 + Year(years) - Hour(1)))
    day = [value(t - t0) / 86_400_000 for t in time]          # days since t0
    hour = [24 * (d - floor(d)) for d in day]
    doy = [mod(d + 273.0, 365.25) for d in day]               # day of year

    seasonal = [-cos(2π * (n - 15) / 365.25) for n in doy]     # coldest mid-January
    diurnal = [cos(2π * (h - 15) / 24) for h in hour]
    weather = [0.7 * sin(2π * d / 5.3) + 0.3 * sin(2π * d / 11.9 + 1.0) for d in day]
    storm = [max(w, 0.0) for w in weather]
    T = [temperature_mean + amplitude * s + 4.0 * c - 3.0 * w
         for (s, c, w) in zip(seasonal, diurnal, weather)]
    storm_sq_mean = sum(s -> s^2, storm) / length(storm)
    P = [precipitation_annual / 8766.0 * s^2 / storm_sq_mean for s in storm]
    cloud = [min(1.0, 0.2 + 1.6 * s) for s in storm]

    lat = deg2rad(40.0)
    μ = map(doy, hour) do n, h
        δ = deg2rad(23.44) * sin(2π * (n - 80) / 365.25)
        sin(lat) * sin(δ) + cos(lat) * cos(δ) * cos(deg2rad(15 * (h - 12)))
    end
    sw = [1361.0 * max(m, 0.0) * 0.75 * (1 - 0.7 * c) for (m, c) in zip(μ, cloud)]
    sza = [acosd(clamp(m, -1.0, 1.0)) for m in μ]
    overcast = [(c - 0.2) / 0.8 for c in cloud]
    σ = 5.67e-8
    lw = [σ * t^4 * (0.70 + 0.26 * c) for (t, c) in zip(T, cloud)]
    vapor = [(0.55 + 0.4 * c) * 611.2 * exp(17.67 * (t - 273.15) / (t - 29.65))
             for (t, c) in zip(T, cloud)]
    wind = [2.5 + 2.0 * s + 0.5 * sin(2π * d / 3.1) for (s, d) in zip(storm, day)]
    pressure = fill(70000.0, length(time))

    return initialize_forcing(time, T, pressure, P, wind, sw, lw, vapor;
        temperature_air_mean=sum(T) / length(T), wind_speed_mean=sum(wind) / length(wind),
        precipitation_mean=sum(P) / years, temperature_observation_height=2.0,
        wind_observation_height=10.0, solar_zenith_angle=sza,
        cloud_optical_thickness=20.0 .* overcast, shortwave_downward_diffuse=sw .* overcast)
end
