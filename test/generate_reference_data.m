% generate_reference_data.m
% Generates reference data from MATLAB GEMB for Julia validation testing.
% Run from the GEMB.jl/test directory.
%
% This script calls individual GEMB functions with controlled inputs
% and saves outputs to .mat files for comparison in Julia tests.

addpath('/Users/gardnera/Documents/GitHub/GEMB/src')

%% Test thermal_conductivity
fprintf('Generating thermal_conductivity reference data...\n')

% Test case 1: Mixed snow and ice
temperature = [260; 255; 250; 245; 270; 268; 265; 260; 255; 250];
density = [300; 400; 500; 600; 700; 800; 910; 910; 910; 910];
ModelParam.density_ice = 910;

ModelParam.thermal_conductivity_method = "Sturm";
K_sturm = thermal_conductivity(temperature, density, ModelParam);

ModelParam.thermal_conductivity_method = "Calonne";
K_calonne = thermal_conductivity(temperature, density, ModelParam);

save('reference_data/thermal_conductivity.mat', 'temperature', 'density', 'K_sturm', 'K_calonne');

%% Test turbulent_heat_flux
fprintf('Generating turbulent_heat_flux reference data...\n')

T_surface = 265.0;
density_air = 1.225;
z0 = 0.00012;
zT = z0 * 0.10;
zQ = z0 * 0.10;

CFS_thf.wind_speed = 5.0;
CFS_thf.pressure_air = 80000.0;
CFS_thf.temperature_air = 268.0;
CFS_thf.vapor_pressure = 300.0;
CFS_thf.temperature_observation_height = 2.0;
CFS_thf.wind_observation_height = 10.0;
CFS_thf.dt = 10800;

[shf, lhf, lh] = turbulent_heat_flux(T_surface, density_air, z0, zT, zQ, CFS_thf);

% Test unstable case
T_surface_unstable = 275.0;
CFS_thf.temperature_air = 260.0;
[shf_unstable, lhf_unstable, lh_unstable] = turbulent_heat_flux(T_surface_unstable, density_air, z0, zT, zQ, CFS_thf);

save('reference_data/turbulent_heat_flux.mat', ...
    'T_surface', 'density_air', 'z0', 'zT', 'zQ', ...
    'shf', 'lhf', 'lh', ...
    'T_surface_unstable', 'shf_unstable', 'lhf_unstable', 'lh_unstable');

%% Test initialize_profile
fprintf('Generating initialize_profile reference data...\n')

MP = model_initialize_parameters();
CF_init = simulate_climate_forcing("test_1", 3);
Profile = model_initialize_profile(MP, CF_init);

dz_init = Profile.dz;
z_center_init = Profile.z_center;
n_layers_init = height(Profile);

save('reference_data/initialize_profile.mat', 'dz_init', 'z_center_init', 'n_layers_init');

%% Test calculate_shortwave_radiation
fprintf('Generating calculate_shortwave_radiation reference data...\n')

n = 10;
dz_sw = ones(n,1) * 0.1;
density_sw = ones(n,1) * 400;
grain_radius_sw = ones(n,1) * 0.5;
albedo_surface_sw = 0.85;
albedo_diffuse_surface_sw = 0.85;

CFS_sw.shortwave_downward = 500;
CFS_sw.shortwave_downward_diffuse = 100;

MP_sw.shortwave_subsurface_absorption = false;
MP_sw.albedo_method = "GardnerSharp";
MP_sw.density_ice = 910;

swf_surface = calculate_shortwave_radiation(dz_sw, density_sw, grain_radius_sw, albedo_surface_sw, albedo_diffuse_surface_sw, CFS_sw, MP_sw);

MP_sw.shortwave_subsurface_absorption = true;
swf_subsurface = calculate_shortwave_radiation(dz_sw, density_sw, grain_radius_sw, albedo_surface_sw, albedo_diffuse_surface_sw, CFS_sw, MP_sw);

save('reference_data/calculate_shortwave_radiation.mat', ...
    'dz_sw', 'density_sw', 'grain_radius_sw', ...
    'albedo_surface_sw', 'albedo_diffuse_surface_sw', ...
    'swf_surface', 'swf_subsurface');

%% Test gemb_core (single timestep integration)
fprintf('Generating gemb_core reference data...\n')

MP_core = model_initialize_parameters();
CF_core = simulate_climate_forcing("test_1", 3);
Profile_core = model_initialize_profile(MP_core, CF_core);

% Extract profile
temperature_core = Profile_core.temperature;
dz_core = Profile_core.dz;
density_core = Profile_core.density;
water_core = Profile_core.water;
grain_radius_core = Profile_core.grain_radius;
grain_dendricity_core = Profile_core.grain_dendricity;
grain_sphericity_core = Profile_core.grain_sphericity;
albedo_core = Profile_core.albedo;
albedo_diffuse_core = Profile_core.albedo_diffuse;

% Build CFS for first timestep
dt_core = (datenum(CF_core.time(2)) - datenum(CF_core.time(1))) * 86400;
dt_core = round(dt_core);
MP_core.dt_divisors = fast_divisors(dt_core * 10000)/10000;

CFS_core.dt = dt_core;
CFS_core.temperature_air = CF_core.temperature_air(1);
CFS_core.pressure_air = CF_core.pressure_air(1);
CFS_core.precipitation = CF_core.precipitation(1);
CFS_core.wind_speed = CF_core.wind_speed(1);
CFS_core.shortwave_downward = CF_core.shortwave_downward(1);
CFS_core.longwave_downward = CF_core.longwave_downward(1);
CFS_core.vapor_pressure = CF_core.vapor_pressure(1);
CFS_core.wind_observation_height = CF_core.Properties.CustomProperties.wind_observation_height;
CFS_core.temperature_observation_height = CF_core.Properties.CustomProperties.temperature_observation_height;
CFS_core.temperature_air_mean = CF_core.Properties.CustomProperties.temperature_air_mean;
CFS_core.wind_speed_mean = CF_core.Properties.CustomProperties.wind_speed_mean;
CFS_core.precipitation_mean = CF_core.Properties.CustomProperties.precipitation_mean;
CFS_core.black_carbon_snow = MP_core.black_carbon_snow;
CFS_core.black_carbon_ice = MP_core.black_carbon_ice;
CFS_core.cloud_optical_thickness = MP_core.cloud_optical_thickness;
CFS_core.solar_zenith_angle = MP_core.solar_zenith_angle;
CFS_core.shortwave_downward_diffuse = MP_core.shortwave_downward_diffuse;
CFS_core.cloud_fraction = MP_core.cloud_fraction;

evaporation_condensation_core = 0;
melt_surface_core = 0;

[temperature_out, dz_out, density_out, water_out, grain_radius_out, grain_dendricity_out, grain_sphericity_out, ...
    albedo_out, albedo_diffuse_out, ec_out, ms_out, sw_net_out, shf_out, lhf_out, lw_up_out, ...
    rain_out, melt_out, runoff_out, refreeze_out, mass_added_out, E_added_out, dc_out, dm_out] = ...
    gemb_core(temperature_core, dz_core, density_core, water_core, grain_radius_core, ...
    grain_dendricity_core, grain_sphericity_core, albedo_core, albedo_diffuse_core, ...
    evaporation_condensation_core, melt_surface_core, CFS_core, MP_core, false);

save('reference_data/gemb_core.mat', ...
    'temperature_core', 'dz_core', 'density_core', 'water_core', ...
    'grain_radius_core', 'grain_dendricity_core', 'grain_sphericity_core', ...
    'albedo_core', 'albedo_diffuse_core', ...
    'temperature_out', 'dz_out', 'density_out', 'water_out', ...
    'grain_radius_out', 'grain_dendricity_out', 'grain_sphericity_out', ...
    'albedo_out', 'albedo_diffuse_out', ...
    'ec_out', 'ms_out', 'sw_net_out', 'shf_out', 'lhf_out', 'lw_up_out', ...
    'rain_out', 'melt_out', 'runoff_out', 'refreeze_out', ...
    'mass_added_out', 'E_added_out', 'dc_out', 'dm_out', ...
    'dt_core');

%% Done
fprintf('All reference data generated successfully.\n')
