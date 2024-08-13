% Load the .mat file containing wind and stress data
load('All_Pier_Data_30_Min_Bins.mat');

% Load the wave and current data from CSV files
waveData = readtable('awac6m_wave_202108to202301_BE.csv');
currentData = readtable('awac6m_current_202108to202301_BE.csv');

% Calculate the correction angle
Correction_Angle = Wdir_deg_N_p - Wdir_measured_p;

% Extract relevant variables from the wave and current data
waveMeanDirection = waveData.waveMeanDirection;
currentDirection = currentData.currentDirection;
windStressTimestamps = datetime(Date_p, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');

% Convert wave and current timestamps
waveTimestamps = datetime(waveData.year, waveData.month, waveData.day, waveData.hour, waveData.minute, 0);
currentTimestamps = datetime(currentData.year, currentData.month, currentData.day, currentData.hour, currentData.minute, 0);

% Ensure data is finite and valid for interpolation
validWaveIndices = ~isnan(waveMeanDirection) & isfinite(waveMeanDirection);
validCurrentIndices = ~isnan(currentDirection) & isfinite(currentDirection);
waveTimestamps = waveTimestamps(validWaveIndices);
waveMeanDirection = waveMeanDirection(validWaveIndices);
currentTimestamps = currentTimestamps(validCurrentIndices);
currentDirection = currentDirection(validCurrentIndices);

% Ensure unique timestamps for interpolation
[uniqueWaveTimestamps, uniqueWaveIndices] = unique(datenum(waveTimestamps));
waveMeanDirection = waveMeanDirection(uniqueWaveIndices);
[uniqueCurrentTimestamps, uniqueCurrentIndices] = unique(datenum(currentTimestamps));
currentDirection = currentDirection(uniqueCurrentIndices);

% Convert wind stress timestamps to datenum
windStressDateNum = datenum(windStressTimestamps);

% Interpolate wave and current directions to match wind and stress data timestamps
waveMeanDirection_interp = interp1(uniqueWaveTimestamps, waveMeanDirection, windStressDateNum, 'linear', 'extrap');
currentDirection_interp = interp1(uniqueCurrentTimestamps, currentDirection, windStressDateNum, 'linear', 'extrap');

% Convert wave and current directions from compass frame to shoreline reference frame
waveMeanDirection_shoreline = waveMeanDirection_interp - Correction_Angle;
currentDirection_shoreline = currentDirection_interp - Correction_Angle;

% Ensure waveMeanDirection_shoreline and currentDirection_shoreline are within -180 to 180
waveMeanDirection_shoreline = mod(waveMeanDirection_shoreline + 180, 360) - 180;
currentDirection_shoreline = mod(currentDirection_shoreline + 180, 360) - 180;

% Define the specific period of interest (September 7th to September 8th)
start_date = datetime(2021, 9, 7);
end_date = datetime(2021, 9, 8, 23, 59, 59);

% Find indices for the specified period
period_indices = find(windStressTimestamps >= start_date & windStressTimestamps <= end_date);

% Extract data for the specified period
period_timestamps = windStressTimestamps(period_indices);
period_Wdir = Wdir_measured_p(period_indices);
period_Sdir = Sdir_measured_p(period_indices);
period_wave = waveMeanDirection_shoreline(period_indices);
period_current = currentDirection_shoreline(period_indices);

% Plot the data for the specified period
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
title('High Stress Period (September 7th to September 8th): Wind, Stress, Wave, and Current Directions');
legend;
grid on;
hold off;