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

% This script is performing performing simulation based on Tim Barfoot's example
% with one calibration parameter.

% angle normalization between -pi and pi
normalizeAngle = @(x) atan2(sin(x), cos(x));

% variance on the bearing measurements
b_var = 6.7143e-04;

% variance on the range measurements
r_var = 9.0036e-04;

% variance on the translational speed measurements
v_var = 0.0044;

% variance on the rotational speed measurements
om_var = 0.0082;

% ground truth offset for laser range finder
d = 0.2190;

% frequency of the simulation
dt = 0.1;

% number of landmarks
nl = 17;

% range for landmarks
min_x_l = 0;
min_y_l = 0;
max_x_l = 10;
max_y_l = 10;

% uniform sampling of landmarks
l = zeros(nl, 2);
l(:, 1) = min_x_l + (max_x_l - min_x_l) .* rand(nl, 1);
l(:, 2) = min_y_l + (max_y_l - min_y_l) .* rand(nl, 1);

% initial pose of the robot
x_0 = 1;
y_0 = 1;
th_0 = pi / 4;

% simulation steps to perform
steps = 10000;

% true translational speed of the robot
v_true = zeros(steps, 1);
v_true(1:end) = 0.2;

% true rotational speed of the robot
om_true = zeros(steps, 1);
om_true(1:300) = 0;
om_true(301:end) = pi / 16;

% measured translational speed of the robot
v = zeros(steps, 1);

% measured rotational speed of the robot
om = zeros(steps, 1);

% range measurements
r = zeros(steps, nl);

% bearing measurements
b = zeros(steps, nl);

% true robot trajectory in x
x_true = zeros(steps, 1);

% true robot trajectory in y
y_true = zeros(steps, 1);

% true robot trajectory in theta
th_true = zeros(steps, 1);

% time steps
t = zeros(steps, 1);

x_true(1) = x_0;
y_true(1) = y_0;
th_true(1) = th_0;

% simulation
for i = 2:steps
  % trajectory
  x_true(i) = x_true(i - 1) + dt * cos(th_true(i - 1)) * v_true(i);
  y_true(i) = y_true(i - 1) + dt * sin(th_true(i - 1)) * v_true(i);
  th_true(i) = th_true(i - 1) + dt * om_true(i);
  th_true(i) = normalizeAngle(th_true(i));

  % measurements
  t(i) = t(i - 1) + dt;
  v(i) = v_true(i) + randn * sqrt(v_var);
  om(i) = om_true(i) + randn * sqrt(om_var);
  ct = cos(th_true(i));
  st = sin(th_true(i));
  for j = 1:nl
    a1 = l(j, 1) - x_true(i) - d * ct;
    b1 = l(j, 2) - y_true(i) - d * st;
    r(i, j) = sqrt(a1^2 + b1^2) + randn * sqrt(r_var);
    b(i, j) = atan2(b1, a1) - th_true(i) + randn * sqrt(b_var);
    b(i, j) = normalizeAngle(b(i, j));
  end
end
