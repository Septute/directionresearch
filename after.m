% Load the .mat file containing wind and stress data
load('All_Pier_Data_30_Min_Bins.mat');

% Load the wave and current data from CSV files
waveData = readtable('awac6m_wave_202108to202301_BE.csv');
currentData = readtable('awac6m_current_202108to202301_BE.csv');

% Define the constant correction angle based on shoreline orientation
Correction_Angle = 12; % degrees

% Define the bin size for the moving mean
bin_size = 1; % You can adjust this value to change the smoothing window size

% Extract relevant variables from the wave and current data
waveMeanDirection = waveData.waveMeanDirection;
currentDirection = currentData.currentDirection;
windStressTimestamps = datetime(Date_p, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');

% Convert wave and current timestamps
waveTimestamps = datetime(waveData.year, waveData.month, waveData.day, waveData.hour, waveData.minute, 0);
currentTimestamps = datetime(currentData.year, currentData.month, currentData.day, currentData.hour, currentData.minute, 0);

% Separate wind and stress data into after July 2022
afterJuly2022_wind = (windStressTimestamps > datetime(2022, 7, 31));

windStressTimestamps_after = windStressTimestamps(afterJuly2022_wind);
Sdir_measured_p_after = Sdir_measured_p(afterJuly2022_wind);
Wdir_measured_p_after = Wdir_measured_p(afterJuly2022_wind);

% Separate wave and current data into after July 2022
waveTimestamps_after = waveTimestamps(waveTimestamps > datetime(2022, 7, 31));
waveMeanDirection_after = waveMeanDirection(waveTimestamps > datetime(2022, 7, 31));

currentTimestamps_after = currentTimestamps(currentTimestamps > datetime(2022, 7, 31));
currentDirection_after = currentDirection(currentTimestamps > datetime(2022, 7, 31));

% Convert data from 0-360 to -180 to 180
convertTo180 = @(x) mod(x + 180, 360) - 180;
waveMeanDirection_after = convertTo180(waveMeanDirection_after);
currentDirection_after = convertTo180(currentDirection_after);

% Ensure unique timestamps for interpolation
[uniqueWaveTimestamps_after, uniqueWaveIndices_after] = unique(datenum(waveTimestamps_after));
waveMeanDirection_after = waveMeanDirection_after(uniqueWaveIndices_after);
[uniqueCurrentTimestamps_after, uniqueCurrentIndices_after] = unique(datenum(currentTimestamps_after));
currentDirection_after = currentDirection_after(uniqueCurrentIndices_after);

% Convert wind stress timestamps to datenum
windStressDateNum_after = datenum(windStressTimestamps_after);

% Interpolate wave and current directions to match wind and stress data timestamps
waveMeanDirection_interp_after = interp1(uniqueWaveTimestamps_after, waveMeanDirection_after, windStressDateNum_after, 'linear', 'extrap');
currentDirection_interp_after = interp1(uniqueCurrentTimestamps_after, currentDirection_after, windStressDateNum_after, 'linear', 'extrap');

% Convert wave and current directions from compass frame to shoreline reference frame
waveMeanDirection_shoreline_after = waveMeanDirection_interp_after - Correction_Angle;
currentDirection_shoreline_after = currentDirection_interp_after - Correction_Angle;

% Apply moving mean to smooth the data
smoothData = @(data, bin_size) movmean(data, bin_size);

waveMeanDirection_shoreline_after_smooth = smoothData(waveMeanDirection_shoreline_after, bin_size);
currentDirection_shoreline_after_smooth = smoothData(currentDirection_shoreline_after, bin_size);
Wdir_measured_p_after_smooth = smoothData(Wdir_measured_p_after, bin_size);
Sdir_measured_p_after_smooth = smoothData(Sdir_measured_p_after, bin_size);

% Define thresholds for high stress direction
high_threshold = 50;  
low_threshold = -50;  

% Find indices where stress direction is notably high or low after July 2022
high_stress_indices_after = find(Sdir_measured_p_after > high_threshold | Sdir_measured_p_after < low_threshold);

% Extract periods of at least 3 days for high stress after July 2022
significant_periods_after = {};
date_nums_after = datenum(windStressTimestamps_after(high_stress_indices_after));
date_diffs_after = diff(date_nums_after);
split_indices_after = find(date_diffs_after > 1);
split_indices_after = [0; split_indices_after'; length(high_stress_indices_after)];

for i = 1:length(split_indices_after)-1
    start_idx = split_indices_after(i) + 1;
    end_idx = split_indices_after(i+1);
    period_indices = high_stress_indices_after(start_idx:end_idx);
    unique_days = unique(floor(date_nums_after(start_idx:end_idx)));
    if length(unique_days) >= 3
        significant_periods_after{end+1} = period_indices;
    end
end

% Create graphs for each significant period after July 2022
for i = 1:length(significant_periods_after)
    period_indices = significant_periods_after{i};
    period_timestamps = windStressTimestamps_after(period_indices);
    period_Wdir = Wdir_measured_p_after_smooth(period_indices);
    period_Sdir = Sdir_measured_p_after_smooth(period_indices);
    period_wave = waveMeanDirection_shoreline_after_smooth(period_indices);
    period_current = currentDirection_shoreline_after_smooth(period_indices);
    
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
    title(['High Stress Period After July 2022: ', num2str(i)]);
    legend;
    grid on;
    hold off;
end