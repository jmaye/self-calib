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
  calibBatch(x_hat, l_hat, theta_hat, u, r, b, t, W, N, maxIter, optTol, ...
  epstol)

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
  epstol = 1e-8;
end

% number of state variables
ns = rows(x_hat);

% non-linear least squares
oldRes = 0;
x_est = x_hat;
l_est = l_hat;
theta_est = theta_hat
for s = 1:maxIter
  % print out iteration number
  s

  % compute the Jacobian and error vector based on the current estimates
  [J, e] = computeJacobian(x_est, l_est, theta_est, t, u, r, b, W, N);

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
  x_est = x_est + [update(1:3:ns * 3) update(2:3:ns * 3) update(3:3:ns * 3)];
  x_est(:, 3) = anglemod(x_est(:, 3));
  l_est = l_est + [update(ns * 3 + 1:2:end - 3)...
    update(ns * 3 + 2:2:end - 3)];
  theta_est = theta_est + update(end - 3 + 1:end);
  theta_est(3) = anglemod(theta_est(3));
  theta_est
end

% compute covariance
[Omega, Sigma, SigmaP, CS, NS, miSum] = solveMarginal(J, cols(J) - 2);
