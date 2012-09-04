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

% This script is performing self-calibration of a laser range finder with 3
% parameters based on Tim Barfoot's example. For the sake of simplicity,
% landmark's positions are not estimated here.

% motion model (without noise)
h = @(x, u, T) [(x(:, 1) + T .* cos(x(:, 3)) .* u(:, 1)), (x(:, 2) +...
  T .* sin(x(:, 3)) .* u(:, 1)), (x(:, 3) + T .* u(:, 2))];

% angle normalization between -pi and pi
normalizeAngle = @(x) atan2(sin(x), cos(x));

% initial odometry
x_odom = zeros(length(x_true), 3);
x_odom(1, :) = [x_true(1), y_true(1), th_true(1)];

% odometry measurements
u = [v, om];

% integrate odometry
for i = 2:length(x_odom)
  x_odom(i, :) = h(x_odom(i - 1, :), u(i, :), t(i) - t(i - 1));
  x_odom(i, 3) = normalizeAngle(x_odom(i, 3));
end

% number of steps to estimate
steps = 1000;

% initial estimates
x_est = x_odom(1:1 + steps, :);
d_x_est = d_x + randn * 0.5
d_y_est = d_y + randn * 0.5
phi_est = phi + randn * 0.5

% count the number of observation rows the matrices have
rowCount = 0;
for i = 2:1 + steps
  for j = 1:rows(l)
    if r(i, j) > 0
      rowCount = rowCount + 1;
    end
  end
end

% matrices allocation
% jacobian matrix
H = sparse(3 * steps + rowCount * 2, 3 * steps + 3 + 3);

% error terms
e = zeros(3 * steps + rowCount * 2, 1);

% motion model covariance
Q = [v_var, 0; 0, om_var];

% observation model covariance
R = [r_var, 0; 0, b_var];

% non-linear least square
tol = 1e-9;
maxNumIter = 200;
oldll = 0;
for s = 1:maxNumIter
  s

  % build matrices
  row = 1;
  col = 1;
  tic;
  for i = 2:1 + steps
    % some pre-computations
    stm1 = sin(x_est(i - 1, 3));
    ctm1 = cos(x_est(i - 1, 3));

    % jacobian of motion model with respect to noise
    H_w = zeros(3, 2);
    H_w(1, 1) = (t(i) - t(i - 1)) * ctm1;
    H_w(2, 1) = (t(i) - t(i - 1)) * stm1;
    H_w(3, 2) = (t(i) - t(i - 1));

    % covariance matrix
    W = H_w * Q * H_w';

    % inverted covariance matrix (simplification because of numerical issues)
    invW = diag(1 ./ diag(W));

    % Cholesky factor of covariance matrix
    invW_sqrt = chol(invW);

    % jacobian of motion model with respect to state variables
    H_x = eye(cols(x_est), cols(x_est));
    H_x(1, 3) = -(t(i) - t(i - 1)) * stm1 * u(i, 1);
    H_x(2, 3) = (t(i) - t(i - 1)) * ctm1 * u(i, 1);
    H_x_cov = -invW_sqrt * H_x;

    % setting everything into H and e
    H(row, col) = H_x_cov(1, 1);
    H(row, col + 2) = H_x_cov(1, 3);
    H(row + 1, col + 1) = H_x_cov(2, 2);
    H(row + 1, col + 2) = H_x_cov(2, 3);
    H(row + 2, col + 2) = H_x_cov(3, 3);
    id_cov = invW_sqrt * eye(cols(x_est), cols(x_est));
    H(row, col + 3) = id_cov(1, 1);
    H(row + 1, col + 4) = id_cov(2, 2);
    H(row + 2, col + 5) = id_cov(3, 3);
    e(row:row + 2) = x_est(i, :)' - h(x_est(i - 1, :), u(i, :), ...
      t(i) - t(i - 1))';
    e(row + 2) = normalizeAngle(e(row + 2));
    e(row:row + 2) = invW_sqrt * e(row:row + 2);
    row = row + cols(x_est);
    col = col + cols(x_est);

    % some pre-computations
    st1 = sin(x_est(i, 3));
    ct1 = cos(x_est(i, 3));

    % covariance matrix of observation model
    N = R;

    % inverted covariance matrix
    invN = diag(1 ./ diag(R));

    % Cholesky factor of covariance matrix
    invN_sqrt = chol(invN);

    % loop over the observations
    for j = 1:rows(l)
      if r(i, j) > 0
        % some pre-computations
        a1 = l(j, 1) - x_est(i, 1) - d_x_est * ct1 + d_y_est * st1;
        b1 = l(j, 2) - x_est(i, 2) - d_x_est * st1 - d_y_est * ct1;
        temp1 = a1^2 + b1^2;
        temp2 = sqrt(temp1);

        % jacobian of observation model with respect to state variables
        G_x = zeros(2, 3);
        G_x(1, 1) = -a1 / temp2;
        G_x(1, 2) = -b1 / temp2;
        G_x(1, 3) = (a1 * (d_x_est * st1 + d_y_est * ct1) + b1 *...
          (-d_x_est * ct1 + d_y_est * st1)) / temp2;
        G_x(2, 1) = b1 / temp1;
        G_x(2, 2) = -a1 / temp1;
        G_x(2, 3) = (a1 * (-d_x_est * ct1 + d_y_est * st1) - b1 *...
          (d_x_est * st1 + d_y_est * ct1)) / temp1 - 1;
        G_x_cov = -invN_sqrt * G_x;

        % jacobian of observation model with respect to calibration parameters
        G_d = zeros(2, 3);
        G_d(1, 1) = -(a1 * ct1 + b1 * st1) / temp2;
        G_d(1, 2) = (a1 * st1 - b1 * ct1) / temp2;
        G_d(2, 1) = (-a1 * st1 + b1 * ct1) / temp1;
        G_d(2, 2) = -(a1 * ct1 + b1 * st1) / temp1;
        G_d(2, 3) = -1;
        G_d_cov = -invN_sqrt * G_d;

        % setting everything into H and e
        H(row, col) = G_x_cov(1, 1);
        H(row, col + 1) = G_x_cov(1, 2);
        H(row, col + 2) = G_x_cov(1, 3);
        H(row + 1, col) = G_x_cov(2, 1);
        H(row + 1, col + 1) = G_x_cov(2, 2);
        H(row + 1, col + 2) = G_x_cov(2, 3);
        H(row, end - 2) = G_d_cov(1, 1);
        H(row, end - 1) = G_d_cov(1, 2);
        H(row + 1, end -2) = G_d_cov(2, 1);
        H(row + 1, end - 1) = G_d_cov(2, 2);
        H(row + 1, end) = G_d_cov(2, 3);
        e(row) = r(i, j) - temp2;
        e(row + 1) = b(i, j) - (atan2(b1, a1) - x_est(i, 3) - phi_est);
        e(row + 1) = normalizeAngle(e(row + 1));
        e(row:row + 1) = invN_sqrt * e(row:row + 1);
        row = row + 2;
      end
    end
  end
  toc;

  % convergence check
  ll = norm(e);
  if oldll == 0
    oldll = ll;
  else
    if abs(oldll - ll) < tol
      break;
    else
      oldll = ll;
    end
  end

  ll

  % compute update
  tic;
  dx = spqr_solve(H, -e, struct('tol', 5));
  toc;

  % update estimate
  x_update = [dx(1:3:3 * steps + 3) dx(2:3:3 * steps + 3)...
    dx(3:3:3 * steps + 3)];
  x_est = x_est + x_update;
  x_est(:, 3) = normalizeAngle(x_est(:, 3));
  d_x_est = d_x_est + dx(end - 2);
  d_y_est = d_y_est + dx(end - 1);
  phi_est = phi_est + dx(end);
  phi_est = normalizeAngle(phi_est);
end

% covariance of the calibration parameters
[C1, R1, P1] = spqr(H, -e, struct('permutation', 'matrix', 'econ', cols(H)));
R1 = P1 * R1 * P1';
covariance = zeros(3, 3);
cov_k = cols(covariance);
for k = cols(H):-1:cols(H) - 2
  temp1 = 0;
  cov_j = cov_k + 1;
  for j = k + 1:cols(H)
    if R1(k, j) ~= 0
      temp1 = temp1 + R1(k, j) * covariance(cov_j, cov_k);
    end
    cov_j = cov_j + 1;
  end
  covariance(cov_k, cov_k) = 1 / R1(k, k) * (1 / R1(k, k) - temp1);
  cov_i = cov_k - 1;
  for i = k - 1:-1:cols(H) - 2
    temp1 = 0;
    temp2 = 0;
    cov_j = cov_i + 1;
    for j = i + 1:k
      if R1(i, j) ~= 0
        temp1 = temp1 + R1(i, j) * covariance(cov_j, cov_k);
      end
      cov_j = cov_j + 1;
    end
    cov_j = cov_k + 1;
    for j = k + 1:cols(H)
      if R1(i, j) ~= 0
        temp2 = temp2 + R1(i, j) * covariance(cov_k, cov_j);
      end
      cov_j = cov_j + 1;
    end
    covariance(cov_i, cov_k) = 1 / R1(i, i) * (-temp1 - temp2);
    covariance(cov_k, cov_i) = covariance(cov_i, cov_k);
    cov_i = cov_i - 1;
  end
  cov_k = cov_k - 1;
end

disp('Final calibration');
d_x_est
d_y_est
phi_est
covariance
