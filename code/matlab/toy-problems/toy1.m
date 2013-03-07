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

% this script generates the toy problem

% steps to simulate
steps = 10;

% parameter for motion model
kappa = 2.5;

% parameter for motion model estimate
kappa_est = kappa + rand * 0.5

% noise for motion model
sigma_m = 1e-6;

% noise for observation model
sigma_o = 1e-5;

% landmark position
l = 5.0;

% landmark position estimate
l_est = l + randn * 0.1

% robot poses
x = zeros(steps, 1);

% robot poses estimates
x_est = zeros(steps, 1);

% landmark observations
y = zeros(steps, 1);

% motion model input
u_i = rand(steps, 1);

% motion model input observations
u_o = zeros(steps, 1);

% start at 0
x(1) = 0;
x_est(1) = 0;

% simulation loop
for i = 2:steps
  % pose update
  x(i) = x(i - 1) + kappa * u_i(i);

  % generate motion input observation
  u_o(i) = u_i(i) + randn * sqrt(sigma_m);

  % generate landmark observation
  y(i) = l - x(i) + randn * sqrt(sigma_o);

  % generate odometry
  x_est(i) = x_est(i - 1) + kappa_est * u_o(i);
end

% error terms
e = zeros(steps - 1, 1);

% Jacobian of motion model with respect to state variable x_k
Hxk = zeros(1, 1);

% Jacobian of motion model with respect to state variable x_{k-1}
Hxkm1 = zeros(1, 1);

% Jacobian of motion model with respect to kappa
Hk = zeros(1, 1);

% Jacobian of observation model with respect to state variable
Gxk = zeros(1, 1);

% Jacobian of observation model with respect to landmark variable
Gl = zeros(1, 1);

% Jacobian initialization
ii = zeros(5 * (steps - 1), 1);
jj = zeros(5 * (steps - 1), 1);
ss = zeros(5 * (steps - 1), 1);

% non-linear least squares
oldRes = 0;
optTol = 1e-9;
for s = 1:20
  % print out iteration number
  s

  % update Jacobian and error term
  row = 1;
  nzcount = 1;

  for i = 2:steps
    % update Jacobians
    Hxk = 1 / kappa_est;
    Hxkm1 = -1 / kappa_est;
    Hk = -(x(i) - x(i - 1)) / kappa_est / kappa_est;
    Gxk = -1;
    Gl = 1;

    % error terms
    e(row) = sqrt(1 / sigma_m) * (u_o(i) - ...
      (x_est(i) - x_est(i - 1)) / kappa_est);
    e(row + 1) = sqrt(1 / sigma_o) * (y(i) - (l_est - x_est(i)));

    % update sparse matrix filling
    ii(nzcount) = row;
    jj(nzcount) = i - 1;
    ss(nzcount) = -sqrt(1 / sigma_m) * Hxkm1;
    nzcount = nzcount + 1;
    ii(nzcount) = row;
    jj(nzcount) = i;
    ss(nzcount) = -sqrt(1 / sigma_m) * Hxk;
    nzcount = nzcount + 1;
    ii(nzcount) = row;
    jj(nzcount) = steps + 1;
    ss(nzcount) = -sqrt(1 / sigma_m) * Hk;
    nzcount = nzcount + 1;
    ii(nzcount) = row + 1;
    jj(nzcount) = i;
    ss(nzcount) = -sqrt(1 / sigma_o) * Gxk;
    nzcount = nzcount + 1;
    ii(nzcount) = row + 1;
    jj(nzcount) = steps + 2;
    ss(nzcount) = -sqrt(1 / sigma_o) * Gl;
    nzcount = nzcount + 1;

    % update row counter
    row = row + 2;
  end

  % build sparse matrix
  H = sparse(ii, jj, ss, (steps - 1) * 2, steps + 2, (steps - 1) * 5);
  norms = colNorm(H); % could be included in the above loop for speedup
  G = spdiags(1 ./ norms, 0, cols(H), cols(H));

  % convergence check
  res = norm(e);
  if oldRes == 0
    oldRes = res;
  else
    if abs(oldRes - res) < optTol
      break;
    else
      oldRes = res;
    end
  end

  % output residual
  res

  % update estimate
  update = G * spqr_solve(H * G, -e, struct('tol', -1));
  x_est = x_est + update(1:steps);
  kappa_est = kappa_est + update(steps + 1)
  l_est = l_est + update(steps + 2)
end
