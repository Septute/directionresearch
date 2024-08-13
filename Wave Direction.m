% Load the .mat file containing wind and stress data
load('All_Pier_Data_30_Min_Bins.mat');

% Load the wave and current data from CSV files
waveData = readtable('awac6m_wave_202108to202301_BE.csv');
currentData = readtable('awac6m_current_202108to202301_BE.csv');

% Calculate the correction angle
Correction_Angle = 12;

% Extract relevant variables from the wave and current data
waveMeanDirection = waveData.waveMeanDirection;
currentDirection = currentData.currentDirection;
windStressTimestamps = datetime(Date_p, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
% Plot wind direction
figure;
plot(windStressTimestamps, Wdir_measured_p, 'LineWidth', 1.5);
xlabel('Time');
ylabel('Wind Direction (degrees)');
title('Wind Direction over Time');
grid on;

% Plot stress direction
figure;
plot(windStressTimestamps, Sdir_measured_p, 'LineWidth', 1.5);
xlabel('Time');
ylabel('Stress Direction (degrees)');
title('Stress Direction over Time');
grid on;


% Convert Date_p to datetime and Date_Num_p for interpolation
windStressTimestamps = datetime(Date_p, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
windStressDateNum = Date_Num_p;

% Ensure data is finite and valid for interpolation
validWaveIndices = ~isnan(waveMeanDirection_shoreline) & isfinite(waveMeanDirection_shoreline) & isfinite(datenum(waveTimestamps));
validCurrentIndices = ~isnan(currentDirection_shoreline) & isfinite(currentDirection_shoreline) & isfinite(datenum(currentTimestamps));
waveTimestamps = datenum(waveTimestamps(validWaveIndices));
waveMeanDirection_shoreline = waveMeanDirection_shoreline(validWaveIndices);
currentTimestamps = datenum(currentTimestamps(validCurrentIndices));
currentDirection_shoreline = currentDirection_shoreline(validCurrentIndices);
% Ensure unique timestamps for interpolation
[uniqueWaveTimestamps, uniqueWaveIndices] = unique(waveTimestamps);
waveMeanDirection_shoreline = waveMeanDirection_shoreline(uniqueWaveIndices);
[uniqueCurrentTimestamps, uniqueCurrentIndices] = unique(currentTimestamps);
currentDirection_shoreline = currentDirection_shoreline(uniqueCurrentIndices);
% Interpolate wave and current directions to match wind and stress data timestamps
waveMeanDirection_interp = interp1(uniqueWaveTimestamps, waveMeanDirection_shoreline, windStressDateNum, 'linear', 'extrap');
currentDirection_interp = interp1(uniqueCurrentTimestamps, currentDirection_shoreline, windStressDateNum, 'linear', 'extrap');
% Convert wave and current directions from compass frame to shoreline reference frame
waveMeanDirection_shoreline = waveMeanDirection_interp - Correction_Angle;
currentDirection_shoreline = currentDirection_interp - Correction_Angle;

% Plot wave direction
figure;
plot(windStressTimestamps, waveMeanDirection_shoreline, 'LineWidth', 1.5);
xlabel('Time');
ylabel('Wave Direction (degrees)');
title('Wave Direction over Time');
grid on;
% Plot current direction
figure;
plot(windStressTimestamps, currentDirection_shoreline, 'LineWidth', 1.5);
xlabel('Time');
ylabel('Current Direction (degrees)');
title('Current Direction over Time');
grid on;


figure;
hold on;
% Plot wind direction
plot(windStressTimestamps, Wdir_measured_p, 'DisplayName', 'Wind Direction');
% Plot stress direction
plot(windStressTimestamps, Sdir_measured_p, 'DisplayName', 'Stress Direction');
% Plot wave direction
plot(windStressTimestamps, waveMeanDirection_shoreline, 'DisplayName', 'Wave Direction');
% Plot current direction
plot(windStressTimestamps, currentDirection_shoreline, 'DisplayName', 'Current Direction');
xlabel('Time');
ylabel('Direction (degrees, shoreline reference frame)');
title('Timeline of Wind, Stress, Wave, and Current Directions');
legend;
grid on;
hold off;



%% REMEMBER WHEN YOU COME BACK THAT THE TIME FOR EACH OF THE DATA MAY BE DIFFERET
%% MAKE TIME GRAPH OF EACH ONTOP OF EACH OTHER AND THAN TAKE ONLY DATA FROM TIME THAT IS THE SAME
%% MAYBE DO THIS TO THE GRAPHS BEFORE TOO

