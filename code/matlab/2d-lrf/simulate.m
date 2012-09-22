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

% This function simulates a robot run.

function [x_true y_true th_true r b v om t] = ...
  simulate(motion, observation, u, x_0, l, Q, R, Theta, steps, T)

% true robot trajectory in x
x_true = zeros(steps, 1);

% true robot trajectory in y
y_true = zeros(steps, 1);

% true robot trajectory in theta
th_true = zeros(steps, 1);

% range measurements
r = zeros(steps, rows(l));

% bearing measurements
b = zeros(steps, rows(l));

% measured translational speed of the robot
v = zeros(steps, 1);

% measured rotational speed of the robot
om = zeros(steps, 1);

% time steps
t = zeros(steps, 1);

% robot start at initial pose
x_true(1) = x_0(1);
y_true(1) = x_0(2);
th_true(1) = x_0(3);

% simulation
for i = 2:steps
  % trajectory
  x_i = motion([x_true(i - 1); y_true(i - 1); th_true(i - 1)], u(i), T);
  x_true(i) = x_i(1);
  y_true(i) = x_i(2);
  th_true(i) = x_i(3);

  % measurements
  t(i) = t(i - 1) + T;
  v(i) = u(i, 1) + randn * sqrt(Q(1, 1));
  om(i) = u(i, 2) + randn * sqrt(Q(2, 2));
  for j = 1:rows(l)
    y_i_j = observation(x_i, l(j, :)', Theta, ...
      [randn * sqrt(R(1, 1)); randn * sqrt(R(2, 2))]);
    r(i, j) = y_i_j(1);
    b(i, j) = y_i_j(2);
  end
end
