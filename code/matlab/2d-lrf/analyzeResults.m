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

function analyzeResults(filename)

load(filename);

dx_ls = zeros(ampResTQRMI(1).rep, cols(ampResTQRMI));
dy_ls = zeros(ampResTQRMI(1).rep, cols(ampResTQRMI));
psi_ls = zeros(ampResTQRMI(1).rep, cols(ampResTQRMI));
dx_ekf = zeros(ampResTQRMI(1).rep, cols(ampResTQRMI));
dy_ekf = zeros(ampResTQRMI(1).rep, cols(ampResTQRMI));
psi_ekf = zeros(ampResTQRMI(1).rep, cols(ampResTQRMI));
dx_tqr_mi = zeros(ampResTQRMI(1).rep, cols(ampResTQRMI));
dy_tqr_mi = zeros(ampResTQRMI(1).rep, cols(ampResTQRMI));
psi_tqr_mi = zeros(ampResTQRMI(1).rep, cols(ampResTQRMI));

for i = 1:cols(ampResTQRMI)
  dx_ls(:, i) = ampResTQRMI(i).Theta_ls(:, 1);
  dy_ls(:, i) = ampResTQRMI(i).Theta_ls(:, 2);
  psi_ls(:, i) = ampResTQRMI(i).Theta_ls(:, 3);
  dx_ekf(:, i) = ampResTQRMI(i).Theta_ekf(:, 1);
  dy_ekf(:, i) = ampResTQRMI(i).Theta_ekf(:, 2);
  psi_ekf(:, i) = ampResTQRMI(i).Theta_ekf(:, 3);
  dx_tqr_mi(:, i) = ampResTQRMI(i).Theta_tqr_mi(:, 1);
  dy_tqr_mi(:, i) = ampResTQRMI(i).Theta_tqr_mi(:, 2);
  psi_tqr_mi(:, i) = ampResTQRMI(i).Theta_tqr_mi(:, 3);
end

save('dx_ls.mat', 'dx_ls');
save('dy_ls.mat', 'dy_ls');
save('psi_ls.mat', 'psi_ls');
save('dx_ekf.mat', 'dx_ekf');
save('dy_ekf.mat', 'dy_ekf');
save('psi_ekf.mat', 'psi_ekf');
save('dx_tqr_mi.mat', 'dx_tqr_mi');
save('dy_tqr_mi.mat', 'dy_tqr_mi');
save('psi_tqr_mi.mat', 'psi_tqr_mi');

boxplot(dx_ls, 'labels', [0; 0.5; 1.0; 1.5; 2.0; 2.5; 3.0; 3.5; 4.0; 4.5; 5.0]);
exportEPSFig(gcf, 'dx_ls.eps');
close;
boxplot(dy_ls, 'labels', [0; 0.5; 1.0; 1.5; 2.0; 2.5; 3.0; 3.5; 4.0; 4.5; 5.0]);
exportEPSFig(gcf, 'dy_ls.eps');
close;
boxplot(psi_ls, 'labels', [0; 0.5; 1.0; 1.5; 2.0; 2.5; 3.0; 3.5; 4.0; 4.5; 5.0]);
exportEPSFig(gcf, 'psi_ls.eps');
close;
boxplot(dx_ekf, 'labels', [0; 0.5; 1.0; 1.5; 2.0; 2.5; 3.0; 3.5; 4.0; 4.5; 5.0]);
exportEPSFig(gcf, 'dx_ekf.eps');
close;
boxplot(dy_ekf, 'labels', [0; 0.5; 1.0; 1.5; 2.0; 2.5; 3.0; 3.5; 4.0; 4.5; 5.0]);
exportEPSFig(gcf, 'dy_ekf.eps');
close;
boxplot(psi_ekf, 'labels', [0; 0.5; 1.0; 1.5; 2.0; 2.5; 3.0; 3.5; 4.0; 4.5; 5.0]);
exportEPSFig(gcf, 'psi_ekf.eps');
close;
boxplot(dx_tqr_mi, 'labels', [0; 0.5; 1.0; 1.5; 2.0; 2.5; 3.0; 3.5; 4.0; 4.5; 5.0]);
exportEPSFig(gcf, 'dx_tqr_mi.eps');
close;
boxplot(dy_tqr_mi, 'labels', [0; 0.5; 1.0; 1.5; 2.0; 2.5; 3.0; 3.5; 4.0; 4.5; 5.0]);
exportEPSFig(gcf, 'dy_tqr_mi.eps');
close;
boxplot(psi_tqr_mi, 'labels', [0; 0.5; 1.0; 1.5; 2.0; 2.5; 3.0; 3.5; 4.0; 4.5; 5.0]);
exportEPSFig(gcf, 'psi_tqr_mi.eps');
close;
