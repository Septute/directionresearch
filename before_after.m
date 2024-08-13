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

% Separate wave and current data into two instances: before April 2022 and after July 2022
beforeApril2022 = (waveTimestamps < datetime(2022, 4, 1));
afterJuly2022 = (waveTimestamps > datetime(2022, 7, 31));

waveTimestamps_before = waveTimestamps(beforeApril2022);
waveMeanDirection_before = waveMeanDirection(beforeApril2022);

currentTimestamps_before = currentTimestamps(currentTimestamps < datetime(2022, 4, 1));
currentDirection_before = currentDirection(currentTimestamps < datetime(2022, 4, 1));

waveTimestamps_after = waveTimestamps(afterJuly2022);
waveMeanDirection_after = waveMeanDirection(afterJuly2022);

currentTimestamps_after = currentTimestamps(currentTimestamps > datetime(2022, 7, 31));
currentDirection_after = currentDirection(currentTimestamps > datetime(2022, 7, 31));

% Separate wind and stress data into two instances: before April 2022 and after July 2022
beforeApril2022_wind = (windStressTimestamps < datetime(2022, 4, 1));
afterJuly2022_wind = (windStressTimestamps > datetime(2022, 7, 31));

windStressTimestamps_before = windStressTimestamps(beforeApril2022_wind);
Sdir_measured_p_before = Sdir_measured_p(beforeApril2022_wind);
Wdir_measured_p_before = Wdir_measured_p(beforeApril2022_wind);

windStressTimestamps_after = windStressTimestamps(afterJuly2022_wind);
Sdir_measured_p_after = Sdir_measured_p(afterJuly2022_wind);
Wdir_measured_p_after = Wdir_measured_p(afterJuly2022_wind);

% Plot the wave and current data before April 2022
figure;
subplot(2, 1, 1);
plot(waveTimestamps_before, waveMeanDirection_before, 'DisplayName', 'Wave Direction Before April 2022');
xlabel('Time');
ylabel('Wave Direction (degrees)');
title('Wave Direction Before April 2022');
legend('show');
grid on;

subplot(2, 1, 2);
plot(currentTimestamps_before, currentDirection_before, 'DisplayName', 'Current Direction Before April 2022');
xlabel('Time');
ylabel('Current Direction (degrees)');
title('Current Direction Before April 2022');
legend('show');
grid on;

% Plot the wave and current data after July 2022
figure;
subplot(2, 1, 1);
plot(waveTimestamps_after, waveMeanDirection_after, 'DisplayName', 'Wave Direction After July 2022');
xlabel('Time');
ylabel('Wave Direction (degrees)');
title('Wave Direction After July 2022');
legend('show');
grid on;

subplot(2, 1, 2);
plot(currentTimestamps_after, currentDirection_after, 'DisplayName', 'Current Direction After July 2022');
xlabel('Time');
ylabel('Current Direction (degrees)');
title('Current Direction After July 2022');
legend('show');
grid on;

% Convert data from 0-360 to -180 to 180
convertTo180 = @(x) mod(x + 180, 360) - 180;
waveMeanDirection_before = convertTo180(waveMeanDirection_before);
currentDirection_before = convertTo180(currentDirection_before);
Sdir_measured_p_before = convertTo180(Sdir_measured_p_before);
Wdir_measured_p_before = convertTo180(Wdir_measured_p_before);

waveMeanDirection_after = convertTo180(waveMeanDirection_after);
currentDirection_after = convertTo180(currentDirection_after);
Sdir_measured_p_after = convertTo180(Sdir_measured_p_after);
Wdir_measured_p_after = convertTo180(Wdir_measured_p_after);

% Ensure data is finite and valid for interpolation
validWaveIndices_before = ~isnan(waveMeanDirection_before) & isfinite(waveMeanDirection_before);
validCurrentIndices_before = ~isnan(currentDirection_before) & isfinite(currentDirection_before);
waveTimestamps_before = waveTimestamps_before(validWaveIndices_before);
waveMeanDirection_before = waveMeanDirection_before(validWaveIndices_before);
currentTimestamps_before = currentTimestamps_before(validCurrentIndices_before);
currentDirection_before = currentDirection_before(validCurrentIndices_before);

validWaveIndices_after = ~isnan(waveMeanDirection_after) & isfinite(waveMeanDirection_after);
validCurrentIndices_after = ~isnan(currentDirection_after) & isfinite(currentDirection_after);
waveTimestamps_after = waveTimestamps_after(validWaveIndices_after);
waveMeanDirection_after = waveMeanDirection_after(validWaveIndices_after);
currentTimestamps_after = currentTimestamps_after(validCurrentIndices_after);
currentDirection_after = currentDirection_after(validCurrentIndices_after);

% Ensure unique timestamps for interpolation
[uniqueWaveTimestamps_before, uniqueWaveIndices_before] = unique(datenum(waveTimestamps_before));
waveMeanDirection_before = waveMeanDirection_before(uniqueWaveIndices_before);
[uniqueCurrentTimestamps_before, uniqueCurrentIndices_before] = unique(datenum(currentTimestamps_before));
currentDirection_before = currentDirection_before(uniqueCurrentIndices_before);

[uniqueWaveTimestamps_after, uniqueWaveIndices_after] = unique(datenum(waveTimestamps_after));
waveMeanDirection_after = waveMeanDirection_after(uniqueWaveIndices_after);
[uniqueCurrentTimestamps_after, uniqueCurrentIndices_after] = unique(datenum(currentTimestamps_after));
currentDirection_after = currentDirection_after(uniqueCurrentIndices_after);

% Convert wind stress timestamps to datenum
windStressDateNum_before = datenum(windStressTimestamps_before);
windStressDateNum_after = datenum(windStressTimestamps_after);

% Interpolate wave and current directions to match wind and stress data timestamps
waveMeanDirection_interp_before = interp1(uniqueWaveTimestamps_before, waveMeanDirection_before, windStressDateNum_before, 'linear', 'extrap');
currentDirection_interp_before = interp1(uniqueCurrentTimestamps_before, currentDirection_before, windStressDateNum_before, 'linear', 'extrap');

waveMeanDirection_interp_after = interp1(uniqueWaveTimestamps_after, waveMeanDirection_after, windStressDateNum_after, 'linear', 'extrap');
currentDirection_interp_after = interp1(uniqueCurrentTimestamps_after, currentDirection_after, windStressDateNum_after, 'linear', 'extrap');

% Convert wave and current directions from compass frame to shoreline reference frame
waveMeanDirection_shoreline_before = waveMeanDirection_interp_before - Correction_Angle;
currentDirection_shoreline_before = currentDirection_interp_before - Correction_Angle;

waveMeanDirection_shoreline_after = waveMeanDirection_interp_after - Correction_Angle;
currentDirection_shoreline_after = currentDirection_interp_after - Correction_Angle;

figure;
plot(windStressTimestamps_before,waveMeanDirection_shoreline_before)