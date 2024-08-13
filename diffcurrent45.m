% Load the data
data = readtable('Pier_Data_BE.csv', 'VariableNamingRule', 'preserve');

% Remove rows where Wdir.range is 'Offshore' or Wdir.deg is 9999
valid_idx = ~strcmp(data.('Wdir.range'), 'Offshore') & data.('Wdir.deg') ~= 9999;

% Apply data quality control
valid_idx = valid_idx & data.Cdz <= 0.01;
valid_idx = valid_idx & ~isnan(data.Cdz) & ~isnan(data.('Wdir.deg')) & ~isnan(data.('Sdir.measured')) & ~isnan(data.Uz);
valid_idx = valid_idx & data.r2_uw >= 0.9;

% Filtered data
filtered_Wdir_measured = data.('Wdir.deg')(valid_idx);
filtered_Sdir_measured = data.('Sdir.measured')(valid_idx);

% Adjust wave and current direction data
valid_wave_direction_idx = valid_idx & ~strcmp(data.('waveMeanDirection_awc4.5'), 'NA');
valid_current_direction_idx = valid_idx & ~strcmp(data.('currentDirection_awc4.5'), 'NA');

filtered_data_wave = data(valid_wave_direction_idx, :);
filtered_data_current = data(valid_current_direction_idx, :);

filtered_waveMeanDirection_awc45 = str2double(filtered_data_wave.('waveMeanDirection_awc4.5'));
filtered_currentDirection_awc45 = str2double(filtered_data_current.('currentDirection_awc4.5'));

% Interpolation of wave and current data
time_idx = find(valid_idx);
wave_time_idx = find(valid_wave_direction_idx);
current_time_idx = find(valid_current_direction_idx);

% Interpolate wave directions
adjusted_wave_direction = interp1(wave_time_idx, filtered_waveMeanDirection_awc45, time_idx, 'linear', 'extrap');

% Interpolate current directions
adjusted_current_direction = interp1(current_time_idx, filtered_currentDirection_awc45, time_idx, 'linear', 'extrap');


% Convert directions from "coming from" into "going to" by adding 180 degrees
adjusted_current_direction = adjusted_current_direction + 180;
% Normalize directions to the range [-180, 180)
adjusted_current_direction(adjusted_current_direction > 180) = adjusted_current_direction(adjusted_current_direction > 180) - 360;

% Adjust directions to align with the coastline (shoreward)
coastline_direction_deg = 72;
adjusted_wave_direction = adjusted_wave_direction - coastline_direction_deg;
adjusted_current_direction = adjusted_current_direction - coastline_direction_deg;
% Ensure directions are within the range [-180, 180)
adjusted_wave_direction(adjusted_wave_direction < -180) = adjusted_wave_direction(adjusted_wave_direction < -180) + 360;
adjusted_wave_direction(adjusted_wave_direction > 180) = adjusted_wave_direction(adjusted_wave_direction > 180) - 360;
adjusted_current_direction(adjusted_current_direction < -180) = adjusted_current_direction(adjusted_current_direction < -180) + 360;
adjusted_current_direction(adjusted_current_direction > 180) = adjusted_current_direction(adjusted_current_direction > 180) - 360;

% Filter data where adjusted_current_direction > 0
number = 0;
idx_current_gt_0 = adjusted_current_direction > number;
Wdir_gt_0 = filtered_Wdir_measured(idx_current_gt_0);
Sdir_gt_0 = filtered_Sdir_measured(idx_current_gt_0);

% Plot wind direction versus stress direction for current direction > 0
createFilledPlot(Wdir_gt_0, Sdir_gt_0, 'Wind Direction vs. Stress Direction (Current Direction = South)');

% Filter data where adjusted_current_direction < 0
idx_current_lt_0 = adjusted_current_direction < number;
Wdir_lt_0 = filtered_Wdir_measured(idx_current_lt_0);
Sdir_lt_0 = filtered_Sdir_measured(idx_current_lt_0);

% Plot wind direction versus stress direction for current direction < 0
createFilledPlot(Wdir_lt_0, Sdir_lt_0, 'Wind Direction vs. Stress Direction (Current Direction = North)');