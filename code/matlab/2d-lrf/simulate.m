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

function [x_true y_true th_true r b v om t] =...
  simulate(u, x_0, l, Q, R, Theta, T, maxRange)

% number of steps inferred from motion commands
steps = rows(u);

% true robot trajectory in x
x_true = zeros(steps, 1);

% true robot trajectory in y
y_true = zeros(steps, 1);

% true robot trajectory in theta
th_true = zeros(steps, 1);

% number of landmarks
nl = rows(l);

% range measurements
r = zeros(steps, nl);

% bearing measurements
b = zeros(steps, nl);

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

% angle normalization between -pi and pi
anglemod = @(x) atan2(sin(x), cos(x));

% simulation
for i = 2:steps
  % trajectory
  x_true(i) = x_true(i - 1) + T * cos(th_true(i - 1)) * u(i, 1);
  y_true(i) = y_true(i - 1) + T * sin(th_true(i - 1)) * u(i, 1);
  th_true(i) = anglemod(th_true(i - 1) + T * u(i, 2));

  % measurements
  t(i) = t(i - 1) + T;
  v(i) = u(i, 1) + randn * sqrt(Q(1, 1));
  om(i) = u(i, 2) + randn * sqrt(Q(2, 2));
  % preliminary computations
  ct = cos(th_true(i));
  st = sin(th_true(i));
  numCalib = length(Theta);
  for j = 1:nl
    if numCalib < 3
      aa = l(j, 1) - x_true(i) - Theta(1) * ct;
      bb = l(j, 2) - y_true(i) - Theta(1) * st;
    else
      aa = l(j, 1) - x_true(i) - Theta(1) * ct + Theta(2) * st;
      bb = l(j, 2) - y_true(i) - Theta(1) * st - Theta(2) * ct;
    end
    range = sqrt(aa^2 + bb^2) + randn * sqrt(R(1, 1));
    if nargin < 8 || (nargin == 8 && range < maxRange)
      r(i, j) = range;
      if numCalib < 3
        b(i, j) = anglemod(atan2(bb, aa) - th_true(i) + randn * sqrt(R(2, 2)));
      else
        b(i, j)= anglemod(atan2(bb, aa) - th_true(i) - Theta(3) + randn *...
          sqrt(R(2, 2)));
      end
    end
  end
end
