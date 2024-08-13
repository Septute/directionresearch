% Load data
load('All_Pier_Data_30_Min_Bins.mat');
waveData = readtable('awac6m_wave_202108to202301_BE.csv');
currentData = readtable('awac6m_current_202108to202301_BE.csv');

% Correction_Angle
Correction_Angle = Wdir_deg_N_p - Wdir_measured_p;
waveMeanDirection = waveData.waveMeanDirection;
currentDirection = currentData.currentDirection;

% Convert Date_p to datetime and Date_Num_p 
windStressTimestamps = datetime(Date_p, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
windStressDateNum = Date_Num_p;
waveTimestamps = datetime(waveData.year, waveData.month, waveData.day, waveData.hour, waveData.minute, 0);
currentTimestamps = datetime(currentData.year, currentData.month, currentData.day, currentData.hour, currentData.minute, 0);

% Ensure data is finite and valid for interpolation
validWaveIndices = ~isnan(waveMeanDirection) & isfinite(waveMeanDirection) & isfinite(datenum(waveTimestamps));
validCurrentIndices = ~isnan(currentDirection) & isfinite(currentDirection) & isfinite(datenum(currentTimestamps));
waveTimestamps = datenum(waveTimestamps(validWaveIndices));
waveMeanDirection = waveMeanDirection(validWaveIndices);
currentTimestamps = datenum(currentTimestamps(validCurrentIndices));
currentDirection = currentDirection(validCurrentIndices);

% Ensure unique timestamps for interpolation
[uniqueWaveTimestamps, uniqueWaveIndices] = unique(waveTimestamps);
waveMeanDirection = waveMeanDirection(uniqueWaveIndices);
[uniqueCurrentTimestamps, uniqueCurrentIndices] = unique(currentTimestamps);
currentDirection = currentDirection(uniqueCurrentIndices);

% Interpolate wave and current directions to match wind and stress data timestamps
waveMeanDirection_interp = interp1(uniqueWaveTimestamps, waveMeanDirection, windStressDateNum, 'linear', 'extrap');
currentDirection_interp = interp1(uniqueCurrentTimestamps, currentDirection, windStressDateNum, 'linear', 'extrap');

% Apply the correction angle after interpolation
waveMeanDirection_shoreline = waveMeanDirection_interp - Correction_Angle;
currentDirection_shoreline = currentDirection_interp - Correction_Angle;

% Ensure waveMeanDirection_shoreline and currentDirection_shoreline are within -180 to 180 without using mod()
waveMeanDirection_shoreline = waveMeanDirection_shoreline - 360 * floor((waveMeanDirection_shoreline + 180) / 360);
currentDirection_shoreline = currentDirection_shoreline - 360 * floor((currentDirection_shoreline + 180) / 360);

% Plotting the corrected wave and current directions to verify
figure;
subplot(2,1,1);
plot(windStressDateNum, waveMeanDirection_shoreline, 'm');
title('Corrected Wave Mean Direction');
xlabel('Timestamp');
ylabel('Direction (degrees)');
ylim([-180 180]);
grid on;

subplot(2,1,2);
plot(windStressDateNum, currentDirection_shoreline, 'b');
title('Corrected Current Direction');
xlabel('Timestamp');
ylabel('Direction (degrees)');
ylim([-180 180]);
grid on;

% Filter the data for wind coming 90 degrees from onshore
onshore90_idx = onshore90_p == 1;
filtered_Wdir_measured_p = Wdir_measured_p(onshore90_idx);
filtered_Sdir_measured = Sdir_measured_p(onshore90_idx);
filtered_waveMeanDirection = waveMeanDirection_shoreline(onshore90_idx);
filtered_currentDirection = currentDirection_shoreline(onshore90_idx);

% Shift wind direction data
shifted_Wdir_measured_p = filtered_Wdir_measured_p;
shifted_Wdir_measured_p(shifted_Wdir_measured_p < 60) = shifted_Wdir_measured_p(shifted_Wdir_measured_p < 60) + 360;

% Find the mean of the shifted wind direction data
mean_shifted_Wdir = mean(shifted_Wdir_measured_p);

% Determine the shift needed to center data between -90 and 90 degrees
shift_needed = 90 - max(shifted_Wdir_measured_p);
shifted_Wdir_measured_p = shifted_Wdir_measured_p + shift_needed;

% Create scatter plot with wave, stress, and current on the same graph
figure;
hold on;
scatter(shifted_Wdir_measured_p, filtered_Sdir_measured, 20, 'o', 'MarkerFaceColor', 'c', 'DisplayName', 'Stress Direction');
scatter(shifted_Wdir_measured_p, filtered_waveMeanDirection, 20, 'o', 'MarkerFaceColor', 'm', 'DisplayName', 'Wave Direction');
scatter(shifted_Wdir_measured_p, filtered_currentDirection, 20, 'o', 'MarkerFaceColor', 'b', 'DisplayName', 'Current Direction');
xlabel('Wind Direction (degrees)');
ylabel('Direction');
title('Wave, Stress, and Current Directions vs Wind Direction');
legend;
grid on;

% Bin the data and plot error bars for stress, wave, and current directions
bin_edges = -90:3:90;
bin_centers = bin_edges(1:end-1) + 1.5;

% Initialize arrays for binned means and errors
binned_means_stress = zeros(size(bin_centers));
binned_errors_stress = zeros(size(bin_centers));
binned_means_wave = zeros(size(bin_centers));
binned_errors_wave = zeros(size(bin_centers));
binned_means_current = zeros(size(bin_centers));
binned_errors_current = zeros(size(bin_centers));

for i = 1:length(bin_centers)
    bin_idx = shifted_Wdir_measured_p >= bin_edges(i) & shifted_Wdir_measured_p < bin_edges(i+1);
    
    % Stress Direction
    bin_data_stress = filtered_Sdir_measured(bin_idx);
    if ~isempty(bin_data_stress)
        binned_means_stress(i) = mean(bin_data_stress);
        sem_stress = std(bin_data_stress) / sqrt(length(bin_data_stress));
        binned_errors_stress(i) = 1.96 * sem_stress;
    end
    
    % Wave Direction
    bin_data_wave = filtered_waveMeanDirection(bin_idx);
    if ~isempty(bin_data_wave)
        binned_means_wave(i) = mean(bin_data_wave);
        sem_wave = std(bin_data_wave) / sqrt(length(bin_data_wave));
        binned_errors_wave(i) = 1.96 * sem_wave;
    end
    
    % Current Direction
    bin_data_current = filtered_currentDirection(bin_idx);
    if ~isempty(bin_data_current)
        binned_means_current(i) = mean(bin_data_current);
        sem_current = std(bin_data_current) / sqrt(length(bin_data_current));
        binned_errors_current(i) = 1.96 * sem_current;
    end
end

% Plot the binned data with error bars
figure;
hold on;
errorbar(bin_centers, binned_means_stress, binned_errors_stress, 'o', 'MarkerFaceColor', 'c', 'MarkerEdgeColor', 'k', 'Color', 'k', 'CapSize', 4, 'LineWidth', 1, 'DisplayName', 'Stress Direction');
errorbar(bin_centers, binned_means_wave, binned_errors_wave, 'o', 'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'k', 'Color', 'k', 'CapSize', 4, 'LineWidth', 1, 'DisplayName', 'Wave Direction');
errorbar(bin_centers, binned_means_current, binned_errors_current, 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k', 'Color', 'k', 'CapSize', 4, 'LineWidth', 1, 'DisplayName', 'Current Direction');
xlabel('Wind Direction (degrees)');
ylabel('Direction');
title('Binned Wave, Stress, and Current Directions vs Wind Direction');
legend;
grid on;
hold off;


