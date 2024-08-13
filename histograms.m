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

% Compute the sum of stress direction and wind direction
sum_stress_wind_direction = filtered_Sdir_measured + filtered_Wdir_measured;
% Ensure the sum is within the range [-180, 180)
sum_stress_wind_direction(sum_stress_wind_direction < -180) = sum_stress_wind_direction(sum_stress_wind_direction < -180) + 360;
sum_stress_wind_direction(sum_stress_wind_direction > 180) = sum_stress_wind_direction(sum_stress_wind_direction > 180) - 360;

% Plot histograms on the same figure with different colors and transparency
figure;
hold on;

% Define histogram properties
binEdges = -120:5:120; % Narrower bins
alphaValue = 0.6;

% Plot wind direction histogram
hist_wind = histogram(filtered_Wdir_measured, 'BinEdges', binEdges, 'FaceColor', 'r', 'FaceAlpha', alphaValue);
legend_entries{1} = 'Wind Direction';

% Plot stress direction histogram
hist_stress = histogram(filtered_Sdir_measured, 'BinEdges', binEdges, 'FaceColor', 'g', 'FaceAlpha', alphaValue);
legend_entries{2} = 'Stress Direction';

% Plot wave direction histogram
hist_wave = histogram(adjusted_wave_direction, 'BinEdges', binEdges, 'FaceColor', 'b', 'FaceAlpha', alphaValue);
legend_entries{3} = 'Wave Direction';

% Plot current direction histogram
hist_current = histogram(adjusted_current_direction, 'BinEdges', binEdges, 'FaceColor', 'm', 'FaceAlpha', alphaValue);
legend_entries{4} = 'Current Direction';

% Plot stress + wind direction histogram
hist_sum = histogram(sum_stress_wind_direction, 'BinEdges', binEdges, 'FaceColor', 'y', 'FaceAlpha', alphaValue);
legend_entries{5} = 'Stress + Wind Direction';

% Add legend and titles
legend(legend_entries, 'Location', 'northwest');
title('Histograms of Wind, Stress, Wave, Current, and Stress+Wind Directions');
xlabel('Direction (degrees)');
ylabel('Frequency');
xlim([-120 120]);
grid on;
hold off;