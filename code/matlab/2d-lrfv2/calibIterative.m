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

% This function performs least squares SLAM and calibration for a robot with a
% laser range finder.

function [x_est, l_est, theta_est, Sigma, SigmaP, CS, NS] =...
  calibIterative(x_hat, l_hat, theta_hat, u, r, b, t, W, N, maxIter, optTol, ...
  epstol, batchSize, miTol)

% angle normalization between -pi and pi
anglemod = @(x) atan2(sin(x), cos(x));

% default values
if nargin < 10
  maxIter = 30;
end
if nargin < 11
  optTol = 1e-6;
end
if nargin < 12
  epstol = 1e-7;
end
if nargin < 13
  batchSize = 200;
end
if nargin < 14
  miTol = 0.5;
end

dataPointer = 1;
x_est = [];
l_est = l_hat;
theta_est = theta_hat;
u_stored = [];
r_stored = [];
b_stored = [];
t_stored = [];
previousMiSum = [];
previousNRank = 0;
while dataPointer < rows(x_hat)
  % collect new batch of data
  x_est_temp = [x_est; x_hat(dataPointer:dataPointer + batchSize - 1, :)];
  l_est_temp = l_est;
  theta_est_temp = theta_est;
  u_stored_temp = [u_stored; u(dataPointer:dataPointer + batchSize - 1, :)];
  r_stored_temp = [r_stored; r(dataPointer:dataPointer + batchSize - 1, :)];
  b_stored_temp = [b_stored; b(dataPointer:dataPointer + batchSize - 1, :)];
  t_stored_temp = [t_stored; t(dataPointer:dataPointer + batchSize - 1, :)];

  % optimize
  oldRes = 0;
  for s = 1:maxIter
    % print out iteration number
    s

    % compute the Jacobian and error vector based on the current estimates
    [J, e, ns] = computeJacobian(x_est_temp, l_est_temp, theta_est_temp, ...
      t_stored_temp, u_stored_temp, r_stored_temp, b_stored_temp, W, N, ...
      batchSize);

    % perform scaling
    G = normCol(J);

    % convergence check
    res = norm(e)
    if oldRes == 0
      oldRes = res;
    else
      if abs(oldRes - res) < optTol
        break;
      else
        oldRes = res;
      end
    end

    % update estimate
    rankTol = 20 * (rows(J) + cols(J)) * epstol * sqrt(max(diag((J * G)' * ...
      (J * G))));
    update = G * spqr_solve(J * G, -e, struct('tol', rankTol));
    x_est_temp = x_est_temp + [update(1:3:ns * 3) update(2:3:ns * 3) ...
      update(3:3:ns * 3)];
    x_est_temp(:, 3) = anglemod(x_est_temp(:, 3));
    l_est_temp = l_est_temp + [update(ns * 3 + 1:2:end - 3)...
      update(ns * 3 + 2:2:end - 3)];
    theta_est_temp = theta_est_temp + update(end - 3 + 1:end);
    theta_est_temp(3) = anglemod(theta_est_temp(3));
    theta_est_temp
  end

  % compute covariance
  [Omega, Sigma, SigmaP, CS, NS, miSum] = solveMarginal(J, cols(J) - 2);

  % compute current numerical rank
  nRank = cols(CS);

  % information test
  if isempty(previousMiSum) || 0.5 * (miSum - previousMiSum) > miTol || ...
      nRank > previousNRank
    x_est = x_est_temp;
    l_est = l_est_temp;
    theta_est = theta_est_temp;
    u_stored = u_stored_temp;
    r_stored = r_stored_temp;
    b_stored = b_stored_temp;
    t_stored = t_stored_temp;
    disp('batch accepted');
    mi = miSum - previousMiSum
  else
    disp('batch rejected');
    mi = miSum - previousMiSum
  end

  if nRank < previousNRank
    disp('----------------------------WARNING--------------------------------');
    disp('-------------------Numerical rank decreased------------------------');
    disp('----------------------------WARNING--------------------------------');
  end

  previousMiSum = miSum;
  previousNRank = nRank;
  dataPointer = dataPointer + batchSize;
end
