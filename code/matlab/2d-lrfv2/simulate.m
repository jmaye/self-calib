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

% This function simulates a robot run.

function [x_true r_true b_true u_noisy r_noisy b_noisy t] =...
  simulate(u_true, x_0, l, W, N, theta, T, maxRange)

% number of steps inferred from motion commands
steps = rows(u_true);

% true robot trajectory
x_true = zeros(steps, 3);

% number of landmarks
nl = rows(l);

% range measurements
r_true = zeros(steps, nl);
r_noisy = zeros(steps, nl);

% bearing measurements
b_true = zeros(steps, nl);
b_noisy = zeros(steps, nl);

% robot's speed
u_noisy = zeros(steps, 2);

% time steps
t = zeros(steps, 1);

% robot start at initial pose
x_true(1, :) = x_0;

% angle normalization between -pi and pi
anglemod = @(x) atan2(sin(x), cos(x));

% simulation
for i = 2:steps
  % trajectory
  dTrans = T * u_true(i, 1);
  dRot = T * u_true(i, 2);
  x_true(i, 1) = x_true(i - 1, 1) + dTrans * cos(x_true(i - 1, 3));
  x_true(i, 2) = x_true(i - 1, 2) + dTrans * sin(x_true(i - 1, 3));
  x_true(i, 3) = anglemod(x_true(i - 1, 3) + dRot);

  % measurements
  t(i) = t(i - 1) + T;
  u_noisy(i, 1) = u_true(i, 1) + randn * sqrt(W(1, 1));
  u_noisy(i, 2) = u_true(i, 2) + randn * sqrt(W(3, 3));

  % preliminary computations
  ct = cos(x_true(i, 3));
  st = sin(x_true(i, 3));
  for j = 1:nl
    aa = l(j, 1) - x_true(i, 1) - theta(1) * ct + theta(2) * st;
    bb = l(j, 2) - x_true(i, 2) - theta(1) * st - theta(2) * ct;
    range = sqrt(aa^2 + bb^2);
    if nargin < 8 || (nargin == 8 && range < maxRange)
      r_true(i, j) = range;
      r_noisy(i, j) = range + randn * sqrt(N(1, 1));
      b_true(i, j) = anglemod(atan2(bb, aa) - x_true(i, 3) - theta(3));
      b_noisy(i, j)= b_true(i, j) + randn * sqrt(N(2, 2));
    end
  end
end
