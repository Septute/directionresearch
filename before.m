% Load the data
data = readtable('Pier_Data_BE.csv', 'VariableNamingRule', 'preserve');
load('All_Pier_Data_30_Min_Bins.mat');
% Remove rows where Wdir.range is 'Offshore' or Wdir.deg is 9999
valid_idx = ~strcmp(data.('Wdir.range'), 'Offshore') & data.('Wdir.deg') ~= 9999;
windStressTimestamps = datetime(Date_p, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');

% Apply data quality control
valid_idx = valid_idx & data.Cdz <= 0.01; % Filter out data where Cdz is lower than or equal to 0.01
valid_idx = valid_idx & ~isnan(data.Cdz) & ~isnan(data.('Wdir.deg')) & ~isnan(data.('Sdir.measured')) & ~isnan(data.Uz);

% Add the additional quality control criteria
valid_idx = valid_idx & data.r2_uw >= 0.9;

% Filter the data for wind coming 90 degrees from onshore
onshore90_idx = onshore90_p == 1 & valid_idx;
filtered_Wdir_measured_p = Wdir_measured_p(onshore90_idx);
filtered_Sdir_measured = Sdir_measured_p(onshore90_idx);
filtered_Uz_p = Uz_p(onshore90_idx);


% Shift wind direction data
shifted_Wdir_measured_p = filtered_Wdir_measured_p;
shifted_Wdir_measured_p(shifted_Wdir_measured_p < 60) = shifted_Wdir_measured_p(shifted_Wdir_measured_p < 60) + 360;

% Find the mean of the shifted wind direction data
mean_shifted_Wdir = mean(shifted_Wdir_measured_p);

% Determine the shift needed to center data between -90 and 90 degrees
shift_needed = 90 - max(shifted_Wdir_measured_p);
shifted_Wdir_measured_p = shifted_Wdir_measured_p + shift_needed;
% Adjust wave and current direction data to be centered around 0 and align with the shoreline
valid_wave_direction_idx = valid_idx & ~strcmp(data.waveMeanDirection_awc6, 'NA');

filtered_data_wave = data(valid_wave_direction_idx, :);
filtered_data_current = data(valid_wave_direction_idx, :);
filtered_waveMeanDirection_awc6 = str2double(filtered_data_wave.waveMeanDirection_awc6);
filtered_currentDirection_awc6 = str2double(filtered_data_current.currentDirection_awc6);

adjusted_wave_direction = mod(filtered_waveMeanDirection_awc6 - 72 + 360, 360);
adjusted_current_direction = mod(filtered_currentDirection_awc6 - 72 + 360, 360);

% Center the wave direction around 0 degrees
adjusted_wave_direction(adjusted_wave_direction > 180) = adjusted_wave_direction(adjusted_wave_direction > 180) - 360;

% Center the current direction around 0 degrees
adjusted_current_direction(adjusted_current_direction > 180) = adjusted_current_direction(adjusted_current_direction > 180) - 360;


% Separate wind and stress data into before April 2022
beforeApril2022_wind = (windStressTimestamps < datetime(2022, 4, 1));

windStressTimestamps_before = windStressTimestamps(beforeApril2022_wind);
Sdir_measured_p_before = Sdir_measured_p(beforeApril2022_wind);
Wdir_measured_p_before = Wdir_measured_p(beforeApril2022_wind);

wavetimestamps = windStressTimestamps_before(valid_wave_direction_idx);
waveMeanDirection_shoreline_before = adjusted_wave_direction;
currentDirection_shoreline_before = adjusted_current_direction;


% Define thresholds for high stress direction
high_threshold = 20;  
low_threshold = -20;  

% Find indices where stress direction is notably high or low before April 2022
high_stress_indices_before = find(Sdir_measured_p_before > high_threshold | Sdir_measured_p_before < low_threshold);

% Extract periods of at least 3 days for high stress before April 2022
significant_periods_before = {};
date_nums_before = datenum(windStressTimestamps_before(high_stress_indices_before));
date_diffs_before = diff(date_nums_before);
split_indices_before = find(date_diffs_before > 1);
split_indices_before = [0; split_indices_before'; length(high_stress_indices_before)];

for i = 1:length(split_indices_before)-1
    start_idx = split_indices_before(i) + 1;
    end_idx = split_indices_before(i+1);
    period_indices = high_stress_indices_before(start_idx:end_idx);
    unique_days = unique(floor(date_nums_before(start_idx:end_idx)));
    if length(unique_days) >= 3
        significant_periods_before{end+1} = period_indices;
    end
end

% Create graphs for each significant period before April 2022
for i = 1:length(significant_periods_before)
    period_indices = significant_periods_before{i};
    period_timestamps = windStressTimestamps_before(period_indices);
    period_Wdir = Wdir_measured_p_before(period_indices);
    period_Sdir = Sdir_measured_p_before(period_indices);
    period_wave = waveMeanDirection_shoreline_before(period_indices);
    period_current = currentDirection_shoreline_before(period_indices);
    
    figure;
    hold on;
    % Plot wind direction
    plot(period_timestamps, period_Wdir, 'DisplayName', 'Wind Direction');
    % Plot stress direction
    plot(period_timestamps, period_Sdir, 'DisplayName', 'Stress Direction');
    % Plot wave direction
    plot(period_timestamps, period_wave, 'DisplayName', 'Wave Direction');
    % Plot current direction
    plot(period_timestamps, period_current, 'DisplayName', 'Current Direction');
    xlabel('Time');
    ylabel('Direction (degrees, shoreline reference frame)');
    title(['High Stress Period Before April 2022: ', num2str(i)]);
    legend;
    grid on;
    hold off;
end
