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
filtered_Wdir_measured = -(data.('Wdir.deg')(valid_idx));
filtered_Sdir_measured = -(data.('Sdir.measured')(valid_idx));

% Adjust wave and current direction data to be centered around 0 and align with the shoreline
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

% Convert directions from "going to" into "coming from"
adjusted_wave_direction(adjusted_wave_direction > 180) = 180 - (360 - adjusted_wave_direction(adjusted_wave_direction > 180));
adjusted_wave_direction(adjusted_wave_direction <= 180) = 180 + adjusted_wave_direction(adjusted_wave_direction <= 180);
adjusted_current_direction(adjusted_current_direction > 180) = 180 - (360 - adjusted_current_direction(adjusted_current_direction > 180));
adjusted_current_direction(adjusted_current_direction <= 180) = 180 + adjusted_current_direction(adjusted_current_direction <= 180);

% Adjust directions to align with the coastline (shoreward)
coastline_direction_deg = 72;
for i = 1:length(adjusted_wave_direction)
    if adjusted_wave_direction(i) >= 0 && adjusted_wave_direction(i) <= 72
        adjusted_wave_direction(i) = adjusted_wave_direction(i) - coastline_direction_deg;
    elseif adjusted_wave_direction(i) > 342 && adjusted_wave_direction(i) <= 360
        adjusted_wave_direction(i) = adjusted_wave_direction(i) + (360 - coastline_direction_deg);
    else
        adjusted_wave_direction(i) = adjusted_wave_direction(i) - coastline_direction_deg;
    end
    % Ensure direction is within 0-360 degrees
    if adjusted_wave_direction(i) < 0
        adjusted_wave_direction(i) = adjusted_wave_direction(i) + 360;
    elseif adjusted_wave_direction(i) >= 360
        adjusted_wave_direction(i) = adjusted_wave_direction(i) - 360;
    end
end

for i = 1:length(adjusted_current_direction)
    if adjusted_current_direction(i) >= 0 && adjusted_current_direction(i) <= 72
        adjusted_current_direction(i) = adjusted_current_direction(i) - coastline_direction_deg;
    elseif adjusted_current_direction(i) > 342 && adjusted_current_direction(i) <= 360
        adjusted_current_direction(i) = adjusted_current_direction(i) + (360 - coastline_direction_deg);
    else
        adjusted_current_direction(i) = adjusted_current_direction(i) - coastline_direction_deg;
    end
    % Ensure direction is within 0-360 degrees
    if adjusted_current_direction(i) < 0
        adjusted_current_direction(i) = adjusted_current_direction(i) + 360;
    elseif adjusted_current_direction(i) >= 360
        adjusted_current_direction(i) = adjusted_current_direction(i) - 360;
    end
end

% Center directions around 0 degrees
adjusted_wave_direction = -(adjusted_wave_direction - 180);
adjusted_current_direction = -(adjusted_current_direction - 180);

% Select data for September 16-18, 2021
start_date = datenum(2021, 9, 16);
end_date = datenum(2021, 9, 18);
selected_dates_idx = date_nums >= start_date & date_nums <= end_date;

% Extract relevant data for the selected dates
time_selected = datetime(data.year(selected_dates_idx), data.month(selected_dates_idx), data.day(selected_dates_idx), data.hour(selected_dates_idx), data.minute(selected_dates_idx), 0);
wind_direction_selected = filtered_Wdir_measured(selected_dates_idx);
wave_direction_selected = adjusted_wave_direction(selected_dates_idx);
current_direction_selected = adjusted_current_direction(selected_dates_idx);
stress_direction_selected = filtered_Sdir_measured(selected_dates_idx);

% Calculate the mean directions
mean_wind_direction = mean(wind_direction_selected);
mean_wave_direction = mean(wave_direction_selected);
mean_current_direction = mean(current_direction_selected);
mean_stress_direction = mean(stress_direction_selected);

% Calculate the mean magnitudes (assuming the magnitude is 1 for each direction)
mean_magnitude = 1;

% Circular plot
figure;
compass(deg2rad(mean_wind_direction), mean_magnitude, 'r');
hold on;
compass(deg2rad(mean_wave_direction), mean_magnitude, 'b');
compass(deg2rad(mean_current_direction), mean_magnitude, 'g');
compass(deg2rad(mean_stress_direction), mean_magnitude, 'm');
legend('Wind Direction', 'Wave Direction', 'Current Direction', 'Stress Direction');
title('Mean Directions (September 16-18, 2021)');
hold off;

