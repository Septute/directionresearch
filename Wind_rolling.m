% Isaiah Sutberry 
% CLEAN ONE
load('All_Pier_Data_30_Min_Bins.mat');

% Apply data quality control
valid_idx = Cz_p <= 0.01; 


mvalid_idx = r2_uw_p >= .9 & valid_idx; 

% Filter the data for wind coming 90 degrees from onshore
onshore90_idx = onshore90_p == 1 & mvalid_idx;
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

% Create the scatter plot
figure;
scatter(shifted_Wdir_measured_p, filtered_Sdir_measured, 20, 'o','MarkerFaceColor','c');
hold on;

% Bin the data into 3-degree bins and calculate mean and 95% confidence intervals
bin_edges = -90:20:90;
bin_centers = bin_edges(1:end-1) + 1;
binned_means = zeros(size(bin_centers));
binned_errors = zeros(size(bin_centers));

y_limits = ylim;
plot([-45 -45], y_limits, '--k', 'LineWidth', 1);
plot([45 45], y_limits, '--k', 'LineWidth', 1);
plot([0 0], y_limits, '--k', 'LineWidth', 1);
text(-67.5, y_limits(1), 'Alongshore', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12, 'Color', 'k', 'FontName', 'Arial');
text(2, y_limits(1), 'Onshore', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12, 'Color', 'k', 'FontName', 'Arial');
text(67.5, y_limits(1), 'Alongshore', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12, 'Color', 'k', 'FontName', 'Arial');

for i = 1:length(bin_centers)
    bin_idx = shifted_Wdir_measured_p >= bin_edges(i) & shifted_Wdir_measured_p < bin_edges(i+1);
    bin_data = filtered_Sdir_measured(bin_idx);
    
    if ~isempty(bin_data)
        binned_means(i) = mean(bin_data);
        sem = std(bin_data) / sqrt(length(bin_data)); 
        binned_errors(i) = 1.96 * sem; 
    end
end

% Plot the binned data with error bars
errorbar(bin_centers, binned_means, binned_errors, 'o', 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k', 'Color', 'k', 'CapSize', 4, 'LineWidth', 1);
title('Stress Direction vs Wind Direction');
xlabel('Wind Direction (degrees)');
ylabel('Stress Direction');
grid on;

% Calculate and plot rolling mean of stress direction over wind direction
rolling_bin_size = 15; % Change this value to your desired bin size
sorted_indices = sortrows([shifted_Wdir_measured_p, filtered_Sdir_measured], 1);
sorted_Wdir_measured_p = sorted_indices(:, 1);
sorted_Sdir_measured = sorted_indices(:, 2);

% Compute rolling mean
rolling_means = movmean(sorted_Sdir_measured, rolling_bin_size);

% Plot rolling mean
figure;
scatter(sorted_Wdir_measured_p, sorted_Sdir_measured, 20, 'o', 'MarkerFaceColor', 'c');
hold on;
plot(sorted_Wdir_measured_p, rolling_means, 'r', 'LineWidth', 2);
xlabel('Wind Direction (degrees)');
ylabel('Stress Direction');
title(['Rolling Mean Stress Direction vs Wind Direction (Bin Size: ' num2str(rolling_bin_size) ')']);
grid on;