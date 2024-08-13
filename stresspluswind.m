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

% Flip data to have south be positive and north negative
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

% Identify periods where stress direction is continuously above or below threshold
threshold = 10;
above_threshold = filtered_Sdir_measured > threshold;
below_threshold = filtered_Sdir_measured < -threshold;

% Combine above and below threshold to find periods
combined_threshold = above_threshold | below_threshold;

% Find contiguous periods of at least 10 data points
min_duration = 4;
contiguous_periods = {};
current_period = [];
for i = 1:length(combined_threshold)
    if combined_threshold(i)
        current_period = [current_period, i];
    else
        if length(current_period) >= min_duration
            contiguous_periods{end+1} = current_period;
        end
        current_period = [];
    end
end
if length(current_period) >= min_duration
    contiguous_periods{end+1} = current_period;
end

% Check if no periods were found and exit if so
if isempty(contiguous_periods)
    disp('No periods found with the specified criteria.');
    return;
end

% Create plots for each identified period
for i = 1:length(contiguous_periods)
    period_indices = contiguous_periods{i};
    
    % Determine the start and end dates for the 3-day window around the period
    start_index = max(1, period_indices(1) - 48);
    end_index = min(length(date_nums), period_indices(end) + 48); 
    
    period_data = data(start_index:end_index, :);
    period_time = datetime(period_data.year, period_data.month, period_data.day, period_data.hour, period_data.minute, 0);
    
    % Extract relevant columns and convert to numeric where needed
    stress_direction = filtered_Sdir_measured(start_index:end_index);
    wave_direction = adjusted_wave_direction(start_index:end_index);
    current_direction = adjusted_current_direction(start_index:end_index);
    wind_direction = filtered_Wdir_measured(start_index:end_index);
    combined_stress_wind_direction = stress_direction + wind_direction;
    
    % Ensure all vectors are the same length
    min_length = min([length(period_time), length(combined_stress_wind_direction), length(wave_direction), length(current_direction), length(wind_direction)]);
    period_time = period_time(1:min_length);
    combined_stress_wind_direction = combined_stress_wind_direction(1:min_length);
    wave_direction = wave_direction(1:min_length);
    current_direction = current_direction(1:min_length);
    wind_direction = wind_direction(1:min_length);
    
    % Create plot
    figure;
    hold on;
    plot(period_time, combined_stress_wind_direction, 'r-o', 'DisplayName', 'Wind+Stress Direction');
    plot(period_time, wave_direction, 'b-s', 'DisplayName', 'Wave Direction');
    plot(period_time, current_direction, 'g-^', 'DisplayName', 'Current Direction');
    plot(period_time, wind_direction, 'm-*', 'DisplayName', 'Wind Direction');
    plot(period_time, stress_direction,'k','DisplayName','Stress Direction')
    xlabel('Time');
    ylabel('Direction (degrees)');
    title(['Directions over Time for Period ', num2str(i)]);
    legend('show');
    hold off;
end