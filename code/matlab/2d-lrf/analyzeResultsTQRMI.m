%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2012 by Jerome Maye                                            %
% jerome.maye@gmail.com                                                        %
%                                                                              %
% This program is free software; you can redistribute it and/or modify         %
% it under the terms of the Lesser GNU General Public License as published by  %
% the Free Software Foundation; either version 3 of the License, or            %
% (at your option) any later version.                                          %
%                                                                              %
% This program is distributed in the hope that it will be useful,              %
% but WITHOUT ANY WARRANTY; without even the implied warranty of               %
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                %
% Lesser GNU General Public License for more details.                          %
%                                                                              %
% You should have received a copy of the Lesser GNU General Public License     %
% along with this program. If not, see <http://www.gnu.org/licenses/>.         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script performs the analysis of the results

function analyzeResultsTQRMI(filename)

dx_true = 0.219;
dy_true = 0.1;
psi_true = pi / 4;

load(filename);

dx_ls = zeros(ampResTQRMI(1).rep, cols(ampResTQRMI));
dy_ls = zeros(ampResTQRMI(1).rep, cols(ampResTQRMI));
psi_ls = zeros(ampResTQRMI(1).rep, cols(ampResTQRMI));
dx_ekf = zeros(ampResTQRMI(1).rep, cols(ampResTQRMI));
dy_ekf = zeros(ampResTQRMI(1).rep, cols(ampResTQRMI));
psi_ekf = zeros(ampResTQRMI(1).rep, cols(ampResTQRMI));
dx_tqr = zeros(ampResTQRMI(1).rep, cols(ampResTQRMI));
dy_tqr = zeros(ampResTQRMI(1).rep, cols(ampResTQRMI));
psi_tqr = zeros(ampResTQRMI(1).rep, cols(ampResTQRMI));

dx_box_data = zeros(ampResTQRMI(1).rep, cols(ampResTQRMI) * 3);
dy_box_data = zeros(ampResTQRMI(1).rep, cols(ampResTQRMI) * 3);
psi_box_data = zeros(ampResTQRMI(1).rep, cols(ampResTQRMI) * 3);

dx_ls_rms = zeros(cols(ampResTQRMI), 1);
dy_ls_rms = zeros(cols(ampResTQRMI), 1);
psi_ls_rms = zeros(cols(ampResTQRMI), 1);
dx_tqr_rms = zeros(cols(ampResTQRMI), 1);
dy_tqr_rms = zeros(cols(ampResTQRMI), 1);
psi_tqr_rms = zeros(cols(ampResTQRMI), 1);
dx_ekf_rms = zeros(cols(ampResTQRMI), 1);
dy_ekf_rms = zeros(cols(ampResTQRMI), 1);
psi_ekf_rms = zeros(cols(ampResTQRMI), 1);

data_usage = zeros(cols(ampResTQRMI), 1);

boxIdx = 1;
for i = 1:cols(ampResTQRMI)
  % boxplot
  dx_ls(:, i) = ampResTQRMI(i).Theta_ls(:, 1);
  dy_ls(:, i) = ampResTQRMI(i).Theta_ls(:, 2);
  psi_ls(:, i) = ampResTQRMI(i).Theta_ls(:, 3);
  dx_ekf(:, i) = ampResTQRMI(i).Theta_ekf(:, 1);
  dy_ekf(:, i) = ampResTQRMI(i).Theta_ekf(:, 2);
  psi_ekf(:, i) = ampResTQRMI(i).Theta_ekf(:, 3);
  dx_tqr(:, i) = ampResTQRMI(i).Theta_tqr_mi(:, 1);
  dy_tqr(:, i) = ampResTQRMI(i).Theta_tqr_mi(:, 2);
  psi_tqr(:, i) = ampResTQRMI(i).Theta_tqr_mi(:, 3);
  dx_box_data(:, boxIdx) = ampResTQRMI(i).Theta_ls(:, 1);
  dx_box_data(:, boxIdx + 1) = ampResTQRMI(i).Theta_ekf(:, 1);
  dx_box_data(:, boxIdx + 2) = ampResTQRMI(i).Theta_tqr_mi(:, 1);
  dy_box_data(:, boxIdx) = ampResTQRMI(i).Theta_ls(:, 2);
  dy_box_data(:, boxIdx + 1) = ampResTQRMI(i).Theta_ekf(:, 2);
  dy_box_data(:, boxIdx + 2) = ampResTQRMI(i).Theta_tqr_mi(:, 2);
  psi_box_data(:, boxIdx) = ampResTQRMI(i).Theta_ls(:, 3);
  psi_box_data(:, boxIdx + 1) = ampResTQRMI(i).Theta_ekf(:, 3);
  psi_box_data(:, boxIdx + 2) = ampResTQRMI(i).Theta_tqr_mi(:, 3);
  boxIdx = boxIdx + 3;

  % rms
  dx_ls_rms(i) = sqrt(1 / ampResTQRMI(1).rep * ...
    sum(dx_true - ampResTQRMI(i).Theta_ls(:, 1)).^2);
  dy_ls_rms(i) = sqrt(1 / ampResTQRMI(1).rep * ...
    sum(dy_true - ampResTQRMI(i).Theta_ls(:, 2)).^2);
  psi_ls_rms(i) = sqrt(1 / ampResTQRMI(1).rep * ...
    sum(psi_true - ampResTQRMI(i).Theta_ls(:, 3)).^2);
  dx_tqr_rms(i) = sqrt(1 / ampResTQRMI(1).rep * ...
    sum(dx_true - ampResTQRMI(i).Theta_tqr_mi(:, 1)).^2);
  dy_tqr_rms(i) = sqrt(1 / ampResTQRMI(1).rep * ...
    sum(dy_true - ampResTQRMI(i).Theta_tqr_mi(:, 2)).^2);
  psi_tqr_rms(i) = sqrt(1 / ampResTQRMI(1).rep * ...
    sum(psi_true - ampResTQRMI(i).Theta_tqr_mi(:, 3)).^2);
  dx_ekf_rms(i) = sqrt(1 / ampResTQRMI(1).rep * ...
    sum(dx_true - ampResTQRMI(i).Theta_ekf(:, 1)).^2);
  dy_ekf_rms(i) = sqrt(1 / ampResTQRMI(1).rep * ...
    sum(dy_true - ampResTQRMI(i).Theta_ekf(:, 2)).^2);
  psi_ekf_rms(i) = sqrt(1 / ampResTQRMI(1).rep * ...
    sum(psi_true - ampResTQRMI(i).Theta_ekf(:, 3)).^2);

  % data usage
  data_usage(i) = mean(ampResTQRMI(i).batchSize_rep);

end

dx_fig = figure;
boxplot(dx_box_data, 'labels', [{''}; {'0'}; {''}; {''}; {'0.5'}; {''}; ...
  {''}; {'1.0'}; {''}; {''}; {'1.5'}; {''}; {''}; {'2.0'}; {''}; {''}; ...
  {'2.5'}; {''}; {''}; {'3.0'}; {''}; {''}; {'3.5'}; {''}; {''}; {'4.0'}; ...
  {''}; {''}; {'4.5'}; {''}; {''}; {'5.0'}; {''}], 'plotstyle', 'compact', ...
  'colors', 'rgb', 'labelOrientation', 'horizontal');
set(gca,'XTickLabel',{' '})
x_range = xlim;
line([x_range(1) x_range(2)], [0.219 0.219], 'LineStyle', '--');
line([3.5 3.5], [0.219 - 0.1; 0.219 + 0.1], 'LineStyle', '--', 'Color', 'k');
line([6.5 6.5], [0.219 - 0.1; 0.219 + 0.1], 'LineStyle', '--', 'Color', 'k');
line([9.5 9.5], [0.219 - 0.1; 0.219 + 0.1], 'LineStyle', '--', 'Color', 'k');
line([12.5 12.5], [0.219 - 0.1; 0.219 + 0.1], 'LineStyle', '--', 'Color', 'k');
line([15.5 15.5], [0.219 - 0.1; 0.219 + 0.1], 'LineStyle', '--', 'Color', 'k');
line([18.5 18.5], [0.219 - 0.1; 0.219 + 0.1], 'LineStyle', '--', 'Color', 'k');
line([21.5 21.5], [0.219 - 0.1; 0.219 + 0.1], 'LineStyle', '--', 'Color', 'k');
line([24.5 24.5], [0.219 - 0.1; 0.219 + 0.1], 'LineStyle', '--', 'Color', 'k');
line([27.5 27.5], [0.219 - 0.1; 0.219 + 0.1], 'LineStyle', '--', 'Color', 'k');
line([30.5 30.5], [0.219 - 0.1; 0.219 + 0.1], 'LineStyle', '--', 'Color', 'k');
ylim([0.219 - 0.1; 0.219 + 0.1]);
ylabel('\delta_x [m]');
dy_fig = figure;
boxplot(dy_box_data, 'labels', [{''}; {'0'}; {''}; {''}; {'0.5'}; {''}; ...
  {''}; {'1.0'}; {''}; {''}; {'1.5'}; {''}; {''}; {'2.0'}; {''}; {''}; ...
  {'2.5'}; {''}; {''}; {'3.0'}; {''}; {''}; {'3.5'}; {''}; {''}; {'4.0'}; ...
  {''}; {''}; {'4.5'}; {''}; {''}; {'5.0'}; {''}], 'plotstyle', 'compact', ...
  'colors', 'rgb', 'labelOrientation', 'horizontal');
set(gca,'XTickLabel',{' '})
x_range = xlim;
line([x_range(1) x_range(2)], [0.1 0.1], 'LineStyle', '--');
line([3.5 3.5], [0.1 - 0.1; 0.1 + 0.1], 'LineStyle', '--', 'Color', 'k');
line([6.5 6.5], [0.1 - 0.1; 0.1 + 0.1], 'LineStyle', '--', 'Color', 'k');
line([9.5 9.5], [0.1 - 0.1; 0.1 + 0.1], 'LineStyle', '--', 'Color', 'k');
line([12.5 12.5], [0.1 - 0.1; 0.1 + 0.1], 'LineStyle', '--', 'Color', 'k');
line([15.5 15.5], [0.1 - 0.1; 0.1 + 0.1], 'LineStyle', '--', 'Color', 'k');
line([18.5 18.5], [0.1 - 0.1; 0.1 + 0.1], 'LineStyle', '--', 'Color', 'k');
line([21.5 21.5], [0.1 - 0.1; 0.1 + 0.1], 'LineStyle', '--', 'Color', 'k');
line([24.5 24.5], [0.1 - 0.1; 0.1 + 0.1], 'LineStyle', '--', 'Color', 'k');
line([27.5 27.5], [0.1 - 0.1; 0.1 + 0.1], 'LineStyle', '--', 'Color', 'k');
line([30.5 30.5], [0.1 - 0.1; 0.1 + 0.1], 'LineStyle', '--', 'Color', 'k');
ylim([0.1 - 0.1; 0.1 + 0.1]);
ylabel('\delta_y [m]');
psi_fig = figure;
boxplot(psi_box_data, 'labels', [{''}; {'0'}; {''}; {''}; {'0.5'}; {''}; ...
  {''}; {'1.0'}; {''}; {''}; {'1.5'}; {''}; {''}; {'2.0'}; {''}; {''}; ...
  {'2.5'}; {''}; {''}; {'3.0'}; {''}; {''}; {'3.5'}; {''}; {''}; {'4.0'}; ...
  {''}; {''}; {'4.5'}; {''}; {''}; {'5.0'}; {''}], 'plotstyle', 'compact', ...
  'colors', 'rgb', 'labelOrientation', 'horizontal');
x_range = xlim;
line([x_range(1) x_range(2)], [pi / 4 pi / 4], 'LineStyle', '--');
line([3.5 3.5], [pi / 4 - 0.1; pi / 4 + 0.1], 'LineStyle', '--', 'Color', 'k');
line([6.5 6.5], [pi / 4 - 0.1; pi / 4 + 0.1], 'LineStyle', '--', 'Color', 'k');
line([9.5 9.5], [pi / 4 - 0.1; pi / 4 + 0.1], 'LineStyle', '--', 'Color', 'k');
line([12.5 12.5], [pi / 4 - 0.1; pi / 4 + 0.1], 'LineStyle', '--', ...
  'Color', 'k');
line([15.5 15.5], [pi / 4 - 0.1; pi / 4 + 0.1], 'LineStyle', '--', ...
  'Color', 'k');
line([18.5 18.5], [pi / 4 - 0.1; pi / 4 + 0.1], 'LineStyle', '--', ...
  'Color', 'k');
line([21.5 21.5], [pi / 4 - 0.1; pi / 4 + 0.1], 'LineStyle', '--', ...
  'Color', 'k');
line([24.5 24.5], [pi / 4 - 0.1; pi / 4 + 0.1], 'LineStyle', '--', ...
  'Color', 'k');
line([27.5 27.5], [pi / 4 - 0.1; pi / 4 + 0.1], 'LineStyle', '--', ...
  'Color', 'k');
line([30.5 30.5], [pi / 4 - 0.1; pi / 4 + 0.1], 'LineStyle', '--', ...
  'Color', 'k');
ylim([pi / 4 - 0.1; pi / 4 + 0.1]);
xlabel('Amplitude [m]');
ylabel('\psi [rad]');

amplitudes = [0; 0.5; 1.0; 1.5; 2.0; 2.5; 3.0; 3.5; 4.0; 4.5; 5.0];
rms_dx_fig = figure;
plot(amplitudes, dx_ls_rms, 'r');
hold on;
plot(amplitudes, dx_tqr_rms, 'g');
plot(amplitudes, dx_ekf_rms, 'b');
xlabel('Amplitude [m]');
ylabel('\delta_x [m]');
ylim([0; 0.2]);
rms_dy_fig = figure;
plot(amplitudes, dy_ls_rms, 'r');
hold on;
plot(amplitudes, dy_tqr_rms, 'g');
plot(amplitudes, dy_ekf_rms, 'b');
xlabel('Amplitude [m]');
ylabel('\delta_y [m]');
ylim([0; 0.2]);
rms_psi_fig = figure;
plot(amplitudes, psi_ls_rms, 'r');
hold on;
plot(amplitudes, psi_tqr_rms, 'g');
plot(amplitudes, psi_ekf_rms, 'b');
xlabel('Amplitude [m]');
ylabel('\psi [rad]');
ylim([0; 0.2]);

data_used_fig = figure;
plot(amplitudes, data_usage, 'g');
xlabel('Amplitude [m]');
ylabel('Data usage [steps]');
