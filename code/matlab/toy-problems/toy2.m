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
steps = 100;

% parameter for motion model
kappa1 = 2.5;

% parameter for motion model estimate
kappa1_est = kappa1 + rand * 0.1

% noise for motion model
sigma_m = 1e-6;

% noise for observation model
sigma_o = 1e-5;

% parameter for observation model
kappa2 = 5;

% parameter for observation model estimate
kappa2_est = kappa2 + rand * 0.1

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
  x(i) = x(i - 1) + kappa1 * u_i(i);

  % generate motion input observation
  u_o(i) = u_i(i) + randn * sqrt(sigma_m);

  % generate landmark observation
  y(i) = 1 / kappa2 * x(i) + randn * sqrt(sigma_o);

  % generate odometry
  x_est(i) = x_est(i - 1) + kappa1_est * u_o(i);
end

% error terms
e = zeros(steps - 1, 1);

% Jacobian of motion model with respect to state variable x_k
Hxk = zeros(1, 1);

% Jacobian of motion model with respect to state variable x_{k-1}
Hxkm1 = zeros(1, 1);

% Jacobian of motion model with respect to kappa1
Hk1 = zeros(1, 1);

% Jacobian of observation model with respect to state variable x_k
Gxk = zeros(1, 1);

% Jacobian of observation model with respect to kappa2
Gk2 = zeros(1, 1);

% Jacobian initialization
ii = zeros(5 * (steps - 1), 1);
jj = zeros(5 * (steps - 1), 1);
ss = zeros(5 * (steps - 1), 1);

% non-linear least squares
oldRes = 0;
optTol = 1e-6;
for s = 1:20
  % print out iteration number
  s

  % update Jacobian and error term
  row = 1;
  nzcount = 1;

  for i = 2:steps
    % update Jacobians
    Hxk = 1;
    Hxkm1 = -1;
    Hk1 = -u_o(i);
    Gxk = -1 / kappa2_est;
    Gk2 = x_est(i) / kappa2_est / kappa2_est;

    % error terms
    e(row) = 1 / kappa1_est * sqrt(1 / sigma_m) * (x_est(i) - x_est(i - 1) - ...
      kappa1_est * u_o(i));
    e(row + 1) = sqrt(1 / sigma_o) * (y(i) - 1 / kappa2_est * x_est(i));

    % update sparse matrix filling
    ii(nzcount) = row;
    jj(nzcount) = i - 1;
    ss(nzcount) = 1 / kappa1_est * sqrt(1 / sigma_m) * Hxkm1;
    nzcount = nzcount + 1;
    ii(nzcount) = row;
    jj(nzcount) = i;
    ss(nzcount) = 1 / kappa1_est * sqrt(1 / sigma_m) * Hxk;
    nzcount = nzcount + 1;
    ii(nzcount) = row;
    jj(nzcount) = steps + 1;
    ss(nzcount) = 1 / kappa1_est * sqrt(1 / sigma_m) * Hk1;
    nzcount = nzcount + 1;
    ii(nzcount) = row + 1;
    jj(nzcount) = i;
    ss(nzcount) = sqrt(1 / sigma_o) * Gxk;
    nzcount = nzcount + 1;
    ii(nzcount) = row + 1;
    jj(nzcount) = steps + 2;
    ss(nzcount) = sqrt(1 / sigma_o) * Gk2;
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
  kappa1_est = kappa1_est + update(steps + 1)
  kappa2_est = kappa2_est + update(steps + 2)
end
