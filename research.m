% Isaiah Sutberry
% PART 1 ALL DATA Regression/Data Removal
load('All_Pier_Data_30_Min_Bins.mat')

idx = onshore90_p == 1;
Cz_p_filtered = Cz_p(idx);
Uz_p_filtered = Uz_p(idx);
Ustar_p_filtered = Ustar_p(idx);
Sdir_measured_p_filtered = Sdir_measured_p(idx);
Wdir_measured_p_filtered = Wdir_measured_p(idx);
% IDK if this one is correct or not
% valid_idx = (Cz_p_filtered <= 0.01) & (Uz_p_filtered <= 20) & ...
%            (r2_Tw_p(idx) >= 0.9) & (r2_uw_p(idx) >= 0.9) & (r2_vw_p(idx) >= 0.9);
% 
% Remove data outliers, based on figure 3 in paper
valid_idx = (Cz_p_filtered <= 0.01);
Cz_p_final = Cz_p_filtered(valid_idx);
Uz_p_final = Uz_p_filtered(valid_idx);
Ustar_p_final = Ustar_p(valid_idx);
Sdir_measured_p_final = Sdir_measured_p_filtered(valid_idx);
Wdir_measured_p_final = Wdir_measured_p_filtered(valid_idx);

figure(100);
scatter(Uz_p_final, Cz_p_final);
xlabel('Wind Speed (m/s)');
ylabel('Drag Coefficient');
title('Drag Coefficient Vs Wind Speed');
grid on;


% Linear Regression for Drag COEF and WIND SPEED
hold on
mdl = fitlm(Uz_p_final, Cz_p_final);
disp(mdl);
x_fit = linspace(min(Uz_p_final), max(Uz_p_final), 100);
[y_fit, y_ci] = predict(mdl, x_fit');
plot(x_fit, y_fit, 'r', 'LineWidth', 2);
plot(x_fit, y_ci(:,1), 'r--', 'LineWidth', 1);
plot(x_fit, y_ci(:,2), 'r--', 'LineWidth', 1);
legend('Data', 'Linear Fit', '95% Confidence Interval');
hold off;
disp(['Linear fit equation: Wind Speed = ' num2str(p(1)) ' * Drag Coef + ' num2str(p(2))]);

% Graph of DRAG VS FRIC
figure(200);
scatter(Cz_p_final,Ustar_p_final)
xlabel('Drag Coefficient')
ylabel('Friction Velocity (m/s)')
title('Fric Vel vs Drag Coef')
% Linear Regression for Drag COEF AND FRIC VELOCITY
hold on
mdl = fitlm(Cz_p_final, Ustar_p_final);
disp(mdl);
x_fit = linspace(min(Cz_p_final), max(Cz_p_final), 100);
[y_fit, y_ci] = predict(mdl, x_fit');
plot(x_fit, y_fit, 'r', 'LineWidth', 2);
plot(x_fit, y_ci(:,1), 'r--', 'LineWidth', 1);
plot(x_fit, y_ci(:,2), 'r--', 'LineWidth', 1);
legend('Data', 'Linear Fit', '95% Confidence Interval');
hold off;
 
% Friction versus Wind Speed
X = [ones(length(Uz_p_final), 1) Uz_p_final];
b = X \ Ustar_p_final;
yfit = X * b;
alpha = 0.05;
n = length(Uz_p_final);
p = 2;
dof = n - p;
tval = tinv(1-alpha/2, dof);
residuals = Ustar_p_final - yfit;
s_res = sqrt(sum(residuals.^2) / dof);
X_inv = inv(X' * X);
s_fit = sqrt(sum((X * X_inv) .* X, 2)) * s_res;
ci_upper = yfit + tval * s_fit;
ci_lower = yfit - tval * s_fit;
figure(300);
scatter(Uz_p_final, Ustar_p_final, 'filled');
hold on;
plot(Uz_p_final, yfit, '-r', 'LineWidth', 2);
plot(Uz_p_final, ci_upper, '--r', 'LineWidth', 1);
plot(Uz_p_final, ci_lower, '--r', 'LineWidth', 1);
xlabel('Wind Speed (m/s)');
ylabel('Friction Velocity (m/s)');
title('Linear Regression of Friction Velocity vs Wind Speed');
legend('Data', 'Regression Line', '95% Confidence Interval');
grid on;
hold off;
% PART 2
% ONSHORE AND OFFSHORE
onshore_idx = onshore90_p == 1 & onshore45_p == 1;
offshore_idx = onshore90_p == 1 & onshore45_p == 0;
Cz_p_onshore = Cz_p(onshore_idx);
Uz_p_onshore = Uz_p(onshore_idx);
Cz_p_offshore = Cz_p(offshore_idx);
Uz_p_offshore = Uz_p(offshore_idx);

% Remove data outliers based on figure 3 in paper for onshore data
valid_onshore_idx = (Cz_p_onshore <= 0.01);
Cz_p_onshore_final = Cz_p_onshore(valid_onshore_idx);
Uz_p_onshore_final = Uz_p_onshore(valid_onshore_idx);
valid_offshore_idx = (Cz_p_offshore <= 0.01);
Cz_p_offshore_final = Cz_p_offshore(valid_offshore_idx);
Uz_p_offshore_final = Uz_p_offshore(valid_offshore_idx);

% Plotting
figure(101);
plot(Uz_p_onshore_final, Cz_p_onshore_final, 'o', 'MarkerSize', 3, 'DisplayName', 'Onshore Data');
hold on;
plot(Uz_p_offshore_final, Cz_p_offshore_final, 'x', 'MarkerSize', 3, 'DisplayName', 'Offshore Data');
xlabel('Wind Speed (m/s)');
ylabel('Drag Coefficient');
title('Filtered Cz_p vs. Uz_p (Onshore and Offshore)');
grid on;

% Linear Regression
p_onshore = polyfit(Uz_p_onshore_final, Cz_p_onshore_final, 1); 
yfit_onshore = polyval(p_onshore, Uz_p_onshore_final);
n_onshore = length(Cz_p_onshore_final); 
yresid_onshore = Cz_p_onshore_final - yfit_onshore; 
SSresid_onshore = sum(yresid_onshore.^2); 
SE_onshore = sqrt(SSresid_onshore / (n_onshore - 2));
t_onshore = tinv(0.975, n_onshore - 2); 
ci_onshore = t_onshore * SE_onshore * sqrt(1/n_onshore + (Uz_p_onshore_final - mean(Uz_p_onshore_final)).^2 / ((n_onshore - 1) * var(Uz_p_onshore_final)));
yfit_onshore_upper = yfit_onshore + ci_onshore;
yfit_onshore_lower = yfit_onshore - ci_onshore;

plot(Uz_p_onshore_final, yfit_onshore, '-r', 'LineWidth', 2, 'DisplayName', 'Onshore Linear Fit');
plot(Uz_p_onshore_final, yfit_onshore_upper, '--r', 'LineWidth', 1, 'DisplayName', 'Onshore 95% CI');
plot(Uz_p_onshore_final, yfit_onshore_lower, '--r', 'LineWidth', 1, 'HandleVisibility', 'off');

p_offshore = polyfit(Uz_p_offshore_final, Cz_p_offshore_final, 1);
yfit_offshore = polyval(p_offshore, Uz_p_offshore_final);
n_offshore = length(Cz_p_offshore_final); 
yresid_offshore = Cz_p_offshore_final - yfit_offshore;
SSresid_offshore = sum(yresid_offshore.^2); 
SE_offshore = sqrt(SSresid_offshore / (n_offshore - 2)); 
t_offshore = tinv(0.975, n_offshore - 2); 
ci_offshore = t_offshore * SE_offshore * sqrt(1/n_offshore + (Uz_p_offshore_final - mean(Uz_p_offshore_final)).^2 / ((n_offshore - 1) * var(Uz_p_offshore_final)));
yfit_offshore_upper = yfit_offshore + ci_offshore;
yfit_offshore_lower = yfit_offshore - ci_offshore;

plot(Uz_p_offshore_final, yfit_offshore, '-b', 'LineWidth', 2, 'DisplayName', 'Offshore Linear Fit');
plot(Uz_p_offshore_final, yfit_offshore_upper, '--b', 'LineWidth', 1, 'DisplayName', 'Offshore 95% CI');
plot(Uz_p_offshore_final, yfit_offshore_lower, '--b', 'LineWidth', 1, 'HandleVisibility', 'off');
legend;
hold off;

SStotal_onshore = (n_onshore - 1) * var(Cz_p_onshore_final);
R2_onshore = 1 - SSresid_onshore/SStotal_onshore;
SStotal_offshore = (n_offshore - 1) * var(Cz_p_offshore_final);
R2_offshore = 1 - SSresid_offshore/SStotal_offshore;

disp(['Onshore R^2: ' num2str(R2_onshore)]);
disp(['Offshore R^2: ' num2str(R2_offshore)]);
disp(['Onshore linear fit equation: Cz_p = ' num2str(p_onshore(1)) ' * Uz_p + ' num2str(p_onshore(2))]);
disp(['Offshore linear fit equation: Cz_p = ' num2str(p_offshore(1)) ' * Uz_p + ' num2str(p_offshore(2))]);

% Histogram
figure(102);
hold on;
histogram(Uz_p_onshore_final, 'FaceColor', 'r', 'DisplayName', 'Onshore Uz_p');
histogram(Uz_p_offshore_final, 'FaceColor', 'b', 'DisplayName', 'Offshore Uz_p');
xlabel('Wind Speed (m/s)');
ylabel('Frequency');
title('Histogram of Wind Speed (Onshore and Offshore)');
legend;
hold off;
%
% Playing with directions
%
convert_to_180 = @(angle) mod(angle + 180, 360) - 180;
onshore_idx = onshore90_p == 1 & onshore45_p == 1;
offshore_idx = onshore90_p == 1 & onshore45_p == 0;

% Apply the indices to filter the data
Wdir_measured_p_onshore = Wdir_measured_p(onshore_idx);
Wdir_measured_p_offshore = Wdir_measured_p(offshore_idx);
Sdir_measured_p_offshore = Sdir_measured_p(offshore_idx);

Wdir_measured_p_onshore = arrayfun(convert_to_180, Wdir_measured_p_onshore);
Wdir_measured_p_offshore = arrayfun(convert_to_180, Wdir_measured_p_offshore);

% Wind Direction/STRESS
Wdir_measured_p_final(Wdir_measured_p_final >= 0 & Wdir_measured_p_final <= 60) = Wdir_measured_p_final(Wdir_measured_p_final >= 0 & Wdir_measured_p_final <= 60) + 360;

figure(103);
scatter(Wdir_measured_p_final, Sdir_measured_p_final, 'filled');
xlabel('Wind Direction (degrees)');
ylabel('Stress Direction (degrees)');
title('Scatter plot of Wind Direction vs Stress Direction');
grid on;

% BINNING THE DATA
bin_edges = 240:2:420;
bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;
binned_wind_direction = discretize(Wdir_measured_p_final, bin_edges);

% AVERAGES
bin_avg_wind_direction = arrayfun(@(i) mean(Wdir_measured_p_final(binned_wind_direction == i)), 1:length(bin_centers));
bin_avg_stress_direction = arrayfun(@(i) mean(Sdir_measured_p_final(binned_wind_direction == i)), 1:length(bin_centers));

figure(104);
scatter(bin_avg_wind_direction, bin_avg_stress_direction, 'filled');
xlabel('Binned Wind Direction (degrees)');
ylabel('Average Stress Direction (degrees)');
title('Scatter plot of Binned Wind Direction vs Average Stress Direction');
grid on;

% Linear Regression
valid_bins = ~isnan(bin_avg_wind_direction) & ~isnan(bin_avg_stress_direction);
bin_avg_wind_direction = bin_avg_wind_direction(valid_bins);
bin_avg_stress_direction = bin_avg_stress_direction(valid_bins);
X = [ones(length(bin_avg_wind_direction), 1) bin_avg_wind_direction'];
b = X\bin_avg_stress_direction';
yfit = X * b;
alpha = 0.05;
n = length(bin_avg_wind_direction);
p = 2;
dof = n - p;
tval = tinv(1-alpha/2, dof);
residuals = bin_avg_stress_direction' - yfit;
s_res = sqrt(sum(residuals.^2) / dof);
X_inv = inv(X' * X);
s_fit = sqrt(sum((X * X_inv) .* X, 2)) * s_res;
ci_upper = yfit + tval * s_fit;
ci_lower = yfit - tval * s_fit;

% PLot linear Regression
figure(105);
scatter(bin_avg_wind_direction, bin_avg_stress_direction, 'filled');
hold on;
plot(bin_avg_wind_direction, yfit, '-r', 'LineWidth', 2);
plot(bin_avg_wind_direction, ci_upper, '--r', 'LineWidth', 1);
plot(bin_avg_wind_direction, ci_lower, '--r', 'LineWidth', 1);
xlabel('Binned Wind Direction (degrees)');
ylabel('Stress Direction (degrees)');
title('Linear Regression of Binned Wind Direction vs Stress Direction with 95% Confidence Interval');
legend('Binned Data', 'Linear Fit', '95% Confidence Interval');
grid on;
hold off;

% Histogram of each bin, to show normality, hopefully
for i = 1:length(bin_centers)
    bin_idx = binned_wind_direction == i;
    if any(bin_idx)
        figure;
        
        subplot(2, 1, 1);
        histogram(Wdir_measured_p_final(bin_idx), 'Normalization', 'pdf');
        title(['Histogram of Wind Direction for Bin Centered at ' num2str(bin_centers(i)) ' degrees']);
        xlabel('Wind Direction (degrees)');
        ylabel('Probability Density');
        
        subplot(2, 1, 2);
        histogram(Sdir_measured_p_final(bin_idx), 'Normalization', 'pdf');
        title(['Histogram of Stress Direction for Bin Centered at ' num2str(bin_centers(i)) ' degrees']);
        xlabel('Stress Direction (degrees)');
        ylabel('Probability Density');

    end
end

% Wave Direction is 0 degrees
bin_diff_stress_wave = bin_avg_stress_direction;
bin_wave = bin_diff_stress_wave(valid_bins);
figure(301);
scatter(bin_avg_wind_direction, bin_wave, 'filled');
yline(0, '--k', 'Wave Direction');
xlabel('Binned Wind Direction (degrees)');
ylabel('Stress Direction (degrees)');
title('Graph with Wave Direction');
grid on;

% Stress - Wind
Stress_Wave_Diff = abs(bin_avg_stress_direction); 
figure(302);
scatter(bin_avg_wind_direction, Stress_Wave_Diff, 'filled');
xlabel('Binned Wind Direction (degrees)');
ylabel('Difference between Wave/Stress Direction (degrees)');
title('Distance Stress Direction is from 0');
grid on;




