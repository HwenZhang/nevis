clear; clc;

dx = [1000 500 250 125 62.5 31.25]; % Grid spacing in meters
r = [8000 7500 7.25e3 6.75e3 6.25e3 5.875e3];
figure;
loglog(dx, r, 'k-', 'LineWidth', 2, 'Marker', 'o', 'MarkerFaceColor', 'k');
hold on
yline(6813, 'r--', 'LineWidth', 2); % Reference line at 6813 m

xlabel('Grid spacing (m)');
ylabel('Blister radius (m)');
title('Blister radius vs grid spacing');
grid on; 
set(gca, 'FontSize', 14);
% Fit a power law to the data
log10dx = log10(dx);
log10r = log10(r);
p = polyfit(log10dx, log10r, 1);
% Calculate the fitted values
fitted_r = 10.^(polyval(p, log10(dx)));
% Plot the fitted line
hold on;
loglog(dx, fitted_r, 'r--', 'LineWidth', 2);
legend('Blister radius', 'Fitted power law', 'Location', 'NorthWest');
% Display the power law exponent
exponent = p(1);
text(0.02, 0.5, sprintf('Power law exponent: %.2f', exponent), 'Units', 'normalized', 'FontSize', 14, 'Color', 'r');
% Save the figure
saveas(gcf, 'blister_radius_vs_grid_spacing.png');