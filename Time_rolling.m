% Load the .mat file containing wind and stress data
load('All_Pier_Data_30_Min_Bins.mat');
windStressTimestamps = datetime(Date_p, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');

% Apply data quality control
valid_idx = Cz_p <= 0.01; 

mvalid_idx = r2_uw_p >= .9 & valid_idx; 

% Filter the data for wind coming 90 degrees from onshore
onshore90_idx = onshore90_p == 1 & mvalid_idx;
filtered_Wdir_measured_p = Wdir_measured_p(onshore90_idx);
filtered_Sdir_measured = Sdir_measured_p(onshore90_idx);
filtered_Uz_p = Uz_p(onshore90_idx);

% Extract corresponding timestamps
filtered_windStressTimestamps = windStressTimestamps(onshore90_idx);

% Plot stress direction
figure;
scatter(filtered_windStressTimestamps, filtered_Sdir_measured);
xlabel('Time');
ylabel('Stress Direction (degrees)');
title('Stress Direction over Time');
grid on;

% Calculate rolling mean of stress direction
binSize = 3; 
Sdir_measured_p_rolling_mean = movmean(filtered_Sdir_measured, binSize);

% Plot rolling mean of stress direction
figure;
plot(filtered_windStressTimestamps, Sdir_measured_p_rolling_mean, 'r', 'LineWidth', 2);
xlabel('Time');
ylabel('Rolling Mean Stress Direction (degrees)');
title(['Rolling Mean Stress Direction over Time (Bin Size: ' num2str(binSize) ')']);
grid on;