% Load the data
data = readtable('Pier_Data_BE.csv', 'VariableNamingRule', 'preserve');

% Remove rows where Wdir.range is 'Offshore' or Wdir.deg is 9999
valid_idx = ~strcmp(data.('Wdir.range'), 'Offshore') & data.('Wdir.deg') ~= 9999;

% Apply data quality control
valid_idx = valid_idx & data.Cdz <= 0.01;
valid_idx = valid_idx & ~isnan(data.Cdz) & ~isnan(data.('Wdir.deg')) & ~isnan(data.('Sdir.measured')) & ~isnan(data.Uz);

% Add the additional quality control criteria
valid_idx = valid_idx & data.r2_uw >= 0.9;

% Remove data from April 2022 to July 2022
date_nums = datenum(data.year, data.month, data.day, data.hour, data.minute, zeros(size(data.hour)));
remove_dates = date_nums >= datenum(2022, 4, 1) & date_nums <= datenum(2022, 7, 31);
valid_idx = valid_idx & ~remove_dates;

% Filtered data
filtered_Wdir_measured = data.('Wdir.deg')(valid_idx);
filtered_Sdir_measured = data.('Sdir.measured')(valid_idx);

% Adjust wave and current direction data
valid_wave_direction_idx = valid_idx & ~strcmp(data.waveMeanDirection_awc6, 'NA');
valid_current_direction_idx = valid_idx & ~strcmp(data.currentDirection_awc6, 'NA');
filtered_data_wave = data(valid_wave_direction_idx, :);
filtered_data_current = data(valid_current_direction_idx, :);
filtered_waveMeanDirection_awc6 = str2double(filtered_data_wave.waveMeanDirection_awc6);
filtered_currentDirection_awc6 = str2double(filtered_data_current.currentDirection_awc6);

% Interpolation of wave and current data
time_idx = find(valid_idx);
wave_time_idx = find(valid_wave_direction_idx);
current_time_idx = find(valid_current_direction_idx);

% Interpolate wave directions
adjusted_wave_direction = interp1(wave_time_idx, filtered_waveMeanDirection_awc6, time_idx, 'linear', 'extrap');

% Interpolate current directions
adjusted_current_direction = interp1(current_time_idx, filtered_currentDirection_awc6, time_idx, 'linear', 'extrap');

% Interpolate wave heights
wave_height = str2double(data.waveHs_awc6);
valid_wave_height_idx = ~isnan(wave_height);
adjusted_wave_height = interp1(find(valid_wave_height_idx), wave_height(valid_wave_height_idx), time_idx, 'linear', 'extrap');

% Interpolate current speeds
current_speed = str2double(data.currentSpeed_awc6);
valid_current_speed_idx = ~isnan(current_speed);
adjusted_current_speed = interp1(find(valid_current_speed_idx), current_speed(valid_current_speed_idx), time_idx, 'linear', 'extrap');

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

% Find times when wave direction and wind direction are approximately 1 degree apart
threshold = 1; % Define threshold for approximation
similar_wave_wind_idx = abs(adjusted_wave_direction - filtered_Wdir_measured) <= threshold;

% Extract corresponding current direction and stress direction
current_direction_at_similar = adjusted_current_direction(similar_wave_wind_idx);
stress_direction_at_similar = filtered_Sdir_measured(similar_wave_wind_idx);

% Analyze the impact of current direction on stress direction
figure;
scatter(current_direction_at_similar, stress_direction_at_similar, 'filled');
xlabel('Current Direction (degrees)');
ylabel('Stress Direction (degrees)');
title('Impact of Current Direction on Stress Direction when Wind and Wave Directions are Approximately 1 Degree Apart');
grid on;

% Display statistical analysis
disp('Statistical Analysis:');
corr_coeff = corr(current_direction_at_similar, stress_direction_at_similar, 'Rows','complete');
fprintf('Correlation coefficient between current direction and stress direction: %.2f\n', corr_coeff);