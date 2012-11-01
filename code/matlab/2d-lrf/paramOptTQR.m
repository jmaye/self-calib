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

% Objective function for TQR parameter optimization at a given amplitude

function f = paramOptTQR(rankTol)

% simulation parameters
params;

% number of repetitions
rep = 20;

% sine wave parameters
frequency = 0.01;
amplitude = 0;

% landmarks parameters
nl = 17;
minX = 0;
maxX = 30;
minY = 0;
maxY = 30;

% output for each repetition
Theta_tqr = zeros(rep, 3);

% optimization parameters
optTol = 1e-6;
maxIter = 20;

% initial guess for Theta
Theta_hat = [0.23; 0.11; 0.8];

% steps for the path
steps = 5000;

warning off all;

% perform evaluation
for i = 1:rep
  i
  rankTol

  % setting up the scene
  u = genSinPath(steps, amplitude, frequency, T);
  l = genLandmarks(nl, [minX maxX], [minY maxY]);
  [x_true y_true th_true r b v om t] = simulate(u, x0, l, Q, R, Theta, T);

  % initial guess
  x_odom = odomInt([x_true(1), y_true(1), th_true(1)], [v, om], t);
  l_hat = initLandmarks(x_odom, Theta_hat, r, b);

  % TQR optimization
  disp('TQR');
  [x_est l_est Theta_est Sigma] = ls_slam_calib(x_odom, l_hat, ...
     Theta_hat, [v, om], r, b, t, Q, R, maxIter, optTol, rankTol);
  Theta_tqr(i, :) = Theta_est';
end

f = sum(var(Theta_tqr)) + ...
  norm(Theta_tqr - repmat([Theta_hat(1), Theta_hat(2), Theta_hat(3)], rep, 1));

warning on all;
