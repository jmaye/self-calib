%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2013 by Jerome Maye                                            %
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

% this script processes a morsel log file and compare with B-spline

% angle normalization between -pi and pi
anglemod = @(x) atan2(sin(x), cos(x));

figure;
plot(data(:, 2),'g');
hold on;
plot(bsplinePoses(:, 2),'r');
ylabel('x');
xnorm = norm(data(:, 2) - bsplinePoses(:, 2))
figure;
plot(data(:, 3),'g');
hold on;
plot(bsplinePoses(:, 3),'r');
ylabel('y');
ynorm = norm(data(:, 3) - bsplinePoses(:, 3))
figure;
plot(data(:, 4),'g');
hold on;
plot(bsplinePoses(:, 4),'r');
ylabel('z');
znorm = norm(data(:, 4) - bsplinePoses(:, 4))
figure;
plot(deg2rad(data(:, 7)),'g');
hold on;
plot(bsplinePoses(:, 5),'r');
ylabel('yaw');
yawnorm = anglemod(norm(data(:, 7) - bsplinePoses(:, 5)))
figure;
plot(deg2rad(data(:, 6)),'g');
hold on;
plot(bsplinePoses(:, 7),'r');
ylabel('pitch');
pitchnorm = anglemod(norm(data(:, 6) - bsplinePoses(:, 7)))
figure;
plot(deg2rad(data(:, 5)),'g');
hold on;
plot(bsplinePoses(:, 6),'r');
ylabel('roll');
rollnorm = anglemod(norm(data(:, 5) - bsplinePoses(:, 6)))
figure;
plot(data(:, 8),'g');
hold on;
plot(bsplinePoses(:, 8),'r');
ylabel('v_x');
vxnorm = norm(data(:, 8) - bsplinePoses(:, 8))
figure;
plot(data(:, 9),'g');
hold on;
plot(bsplinePoses(:, 9),'r');
ylabel('v_y');
vynorm = norm(data(:, 9) - bsplinePoses(:, 9))
figure;
plot(data(:, 10),'g');
hold on;
plot(bsplinePoses(:, 10),'r');
ylabel('v_z');
vznorm = norm(data(:, 10) - bsplinePoses(:, 10))
figure;
plot(deg2rad(data(:, 11)),'g');
hold on;
plot(bsplinePoses(:, 11),'r');
ylabel('\omega_x');
wxnorm = anglemod(norm(data(:, 11) - bsplinePoses(:, 11)))
figure;
plot(deg2rad(data(:, 12)),'g');
hold on;
plot(bsplinePoses(:, 12),'r');
ylabel('\omega_y');
wynorm = anglemod(norm(data(:, 12) - bsplinePoses(:, 12)))
figure;
plot(deg2rad(data(:, 13)),'g');
hold on;
plot(bsplinePoses(:, 13),'r');
ylabel('\omega_z');
wznorm = anglemod(norm(data(:, 13) - bsplinePoses(:, 13)))

v_w = zeros(rows(data), 3);
om_w = zeros(rows(data), 3);
for i = 1:rows(data)
  % rotation from body to world
  C_wb_roll = deg2rad(data(i, 5));
  C_wb_pitch = deg2rad(data(i, 6));
  C_wb_yaw = deg2rad(data(i, 7));
  cx = cos(C_wb_roll);
  sx = sin(C_wb_roll);
  cy = cos(C_wb_pitch);
  sy = sin(C_wb_pitch);
  cz = cos(C_wb_yaw);
  sz = sin(C_wb_yaw);
  C_wb_x = [1, 0, 0; 0, cx, -sx; 0, sx, cx];
  C_wb_y = [cy, 0, sy; 0, 1, 0; -sy, 0, cy];
  C_wb_z = [cz, -sz, 0; sz, cz, 0; 0, 0, 1];
  C_wb = C_wb_z * C_wb_x * C_wb_y;
  v_w(i, :) = C_wb * data(i, 8:10)';
  om_w(i, :) = C_wb * deg2rad(data(i, 11:13))';
end

figure;
plot(v_w(:, 1),'g');
hold on;
plot(bsplinePoses(:, 14),'r');
ylabel('vw_x');
figure;
plot(v_w(:, 2),'g');
hold on;
plot(bsplinePoses(:, 15),'r');
ylabel('vw_y');
figure;
plot(v_w(:, 3),'g');
hold on;
plot(bsplinePoses(:, 16),'r');
ylabel('vw_z');
figure;
plot(om_w(:, 1),'g');
hold on;
plot(bsplinePoses(:, 17),'r');
ylabel('\omega_x_w');
figure;
plot(om_w(:, 2),'g');
hold on;
plot(bsplinePoses(:, 18),'r');
ylabel('\omega_y_w');
figure;
plot(om_w(:, 3),'g');
hold on;
plot(bsplinePoses(:, 19),'r');
ylabel('\omega_z_w');
