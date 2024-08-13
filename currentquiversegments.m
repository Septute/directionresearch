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

% Ensure the sizes match
adjusted_current_direction = adjusted_current_direction(:);
adjusted_current_speed = adjusted_current_speed(:);

% Convert polar to Cartesian coordinates for quiver plot
u = adjusted_current_speed .* cosd(adjusted_current_direction); % x-component (current speed)
v = adjusted_current_speed .* sind(adjusted_current_direction); % y-component (directionality)

% Check dimensions of variables
disp(['Size of date_nums(valid_idx): ', num2str(length(date_nums(valid_idx)))]);
disp(['Size of u: ', num2str(length(u))]);
disp(['Size of v: ', num2str(length(v))]);

% Ensure u, v, and date_nums(valid_idx) have the same size
if length(date_nums(valid_idx)) == length(u) && length(date_nums(valid_idx)) == length(v)
    disp('Sizes match, proceeding with plot.');
else
    error('Size mismatch detected, please check the data.');
end

% Plot only a small segment to check
start_idx = 1;
end_idx = 1000; % Adjust as needed

% Correct the segment selection to ensure consistent indexing
segment_date_nums = date_nums(valid_idx);
segment_u = u;
segment_v = v;

% Ensure the selected segment sizes match
disp(['Size of segment_date_nums: ', num2str(length(segment_date_nums(start_idx:end_idx)))]);
disp(['Size of segment_u: ', num2str(length(segment_u(start_idx:end_idx)))]);
disp(['Size of segment_v: ', num2str(length(segment_v(start_idx:end_idx)))]);

% Create the quiver plot for the specified segment
figure;
quiver(segment_date_nums(start_idx:end_idx), zeros(size(segment_date_nums(start_idx:end_idx))), segment_u(start_idx:end_idx), segment_v(start_idx:end_idx), 0.5);
datetick('x', 'keepticks');
xlabel('Date');
ylabel('Current Speed and Direction');
title('Current Speed and Direction Over Time (Segment)');
grid on;