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

% Adjust current direction data to be centered around 0 and align with the shoreline
valid_current_direction_idx = valid_idx & ~strcmp(data.currentDirection_awc6, 'NA');
filtered_data_current = data(valid_current_direction_idx, :);
filtered_currentDirection_awc6 = str2double(filtered_data_current.currentDirection_awc6);

% Interpolation of current data
time_idx = find(valid_idx);
current_time_idx = find(valid_current_direction_idx);

% Interpolate current directions
adjusted_current_direction = interp1(current_time_idx, filtered_currentDirection_awc6, time_idx, 'linear', 'extrap');

% Convert directions from "going to" into "coming from"
adjusted_current_direction(adjusted_current_direction > 180) = 180 - (360 - adjusted_current_direction(adjusted_current_direction > 180));
adjusted_current_direction(adjusted_current_direction <= 180) = 180 + adjusted_current_direction(adjusted_current_direction <= 180);

% Adjust directions to align with the coastline (shoreward)
coastline_direction_deg = 72;
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
% Interpolate current speeds
current_speed = str2double(data.currentSpeed_awc6);
valid_current_speed_idx = ~isnan(current_speed);
adjusted_current_speed = interp1(find(valid_current_speed_idx), current_speed(valid_current_speed_idx), time_idx, 'linear', 'extrap');

% Center directions around 0 degrees
adjusted_current_direction = -(adjusted_current_direction - 180);

% Create the quiver plot for current direction over time
figure;
quiver(date_nums(valid_idx), zeros(size(date_nums(valid_idx))), u, v);
datetick('x', 'keepticks');
xlabel('Date');
ylabel('Current Direction');
title('Current Direction and Strength Over Time');
grid on;