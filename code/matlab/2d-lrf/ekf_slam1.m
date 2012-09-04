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

% This script is performing SLAM based on Tim Barfoot's example using
% an Extended Kalman Filter (EKF). Laser range finder
% offset is considered as given.

% motion model (without noise)
h = @(x, u, T) [(x(:, 1) + T .* cos(x(:, 3)) .* u(:, 1)), (x(:, 2) +...
  T .* sin(x(:, 3)) .* u(:, 1)), (x(:, 3) + T .* u(:, 2))];

% angle normalization between -pi and pi
normalizeAngle = @(x) atan2(sin(x), cos(x));

% odometry measurements
u = [v, om];

% number of steps to estimate
steps = 5000;

% map initialization
l_init = l + randn * 0.1;
l_init = zeros(rows(l), 2);
x_map_init = h([x_true(1), y_true(1), th_true(1)], u(2, :), t(2) - t(1));
for i = 1:rows(l)
  l_init(i, 1) = x_map_init(1) + d * cos(x_map_init(3)) + r(2, i) *...
    cos(x_map_init(3) + b(2, i));
  l_init(i, 2) = x_map_init(2) + d * sin(x_map_init(3)) + r(2, i) *...
    sin(x_map_init(3) + b(2, i));
end

% odometry
x_odom = zeros(steps, 3);
x_odom(1, :) = [x_true(1), y_true(1), th_true(1)];

% filtered state
x = zeros(steps, 3 + rows(l) * 2);
x(1, :) = [x_true(1), y_true(1), th_true(1), ...
  reshape(l_init', rows(l) * cols(l), 1)'];

% filtered covariance
P = zeros(3 + rows(l) * 2, 3 + rows(l) * 2, steps);
P(1:3, 1:3, 1) = diag([1, 1, 0.1]);
P(4:end, 4:end, 1) = eye(rows(l) * 2, rows(l) * 2) * 0.1;

% motion model covariance
Q = [v_var, 0; 0, om_var];

% observation model covariance
R = [r_var, 0; 0, b_var];

% matrices allocation
H_x = eye(3, 3);
H_w = zeros(3, 2);

% perform filtering
for k = 2:steps
  % some pre-computations
  stm1 = sin(x(k - 1, 3));
  ctm1 = cos(x(k - 1, 3));

  % jacobian of motion model with respect to state variables
  H_x(1, 3) = -(t(k) - t(k - 1)) * stm1 * u(k, 1);
  H_x(2, 3) = (t(k) - t(k - 1)) * ctm1 * u(k, 1);

  % jacobian of motion model with respect to noise
  H_w(1, 1) = (t(k) - t(k - 1)) * ctm1;
  H_w(2, 1) = (t(k) - t(k - 1)) * stm1;
  H_w(3, 2) = (t(k) - t(k - 1));

  % predictor equations
  Pm_kRR = H_x * P(1:3, 1:3, k - 1) * H_x' + H_w * Q * H_w';
  Pm_kRM = H_x * P(1:3, 4:end, k - 1);
  Pm_k = P(:, :, k - 1);
  Pm_k(1:3, 1:3) = Pm_kRR;
  Pm_k(1:3, 4:end) = Pm_kRM;
  Pm_k(4:end, 1:3) = Pm_kRM';
  xm_kR = h(x(k - 1, 1:3), u(k, :), t(k) - t(k - 1));
  xm_k = x(k - 1, :);
  xm_k(1:3) = xm_kR;
  x_odom(k, :) = h(x_odom(k - 1, :), u(k, :), t(k) - t(k - 1));

  % number of observations
  L_k = nnz(r(k, :));

  % loop over observations if any
  if L_k > 0
    % matrices allocation
    y_k = zeros(2 * L_k, 1);
    g_k = zeros(2 * L_k, 1);
    R_k = zeros(2 * L_k, 2 * L_k);
    G_xk = zeros(2 * L_k, 3 + rows(l) * 2);
    innovation = zeros(2 * L_k, 1);
    y = zeros(2, 1);
    g = zeros(2, 1);
    G_x = zeros(2, 3);
    G_l = zeros(2, 2);

    % some pre-computations
    st1 = sin(xm_k(3));
    ct1 = cos(xm_k(3));

    idx = 1;
    for j = 1:rows(l)
      if r(k, j) > 0
        % some pre-computations
        a1 = xm_k(4 + (j - 1) * 2) - xm_k(1) - d * ct1;
        b1 = xm_k(4 + (j - 1) * 2 + 1) - xm_k(2) - d * st1;
        temp1 = a1^2 + b1^2;
        temp2 = sqrt(temp1);

        % update observation vector
        y(1) = r(k, j);
        y(2) = b(k, j);
        y_k(idx:idx + 1) = y;

        % update prediction vector
        g(1) = temp2;
        g(2) = atan2(b1, a1) - xm_k(3);
        g_k(idx:idx + 1) = g;

        % jacobian of observation model with respect to state variables
        G_x(1, 1) = -a1 / temp2;
        G_x(1, 2) = -b1 / temp2;
        G_x(1, 3) = (a1 * d * st1 - b1 * d * ct1) / temp2;
        G_x(2, 1) = b1 / temp1;
        G_x(2, 2) = -a1 / temp1;
        G_x(2, 3) = -(a1 * d * ct1 + b1 * d * st1) / temp1 - 1;
        G_xk(idx:idx + 1, 1:3) = G_x;

        % jacobian of observation model with respect to landmark positions
        G_l(1, 1) = a1 / temp2;
        G_l(1, 2) = b1 / temp2;
        G_l(2, 1) = -b1 / temp1;
        G_l(2, 2) = a1 / temp1;
        G_xk(idx:idx + 1, 4 + (j - 1) * 2:4 + (j - 1) * 2 + 1) = G_l;

        % compute innovation
        innovation(idx:idx + 1) = y - g;
        innovation(idx + 1) = normalizeAngle(innovation(idx + 1));

        % observation covariance
        R_k(idx:idx + 1, idx:idx + 1) = R;

        % update count
        idx = idx + 2;
      end
    end

    % kalman gain
    K_k = Pm_k * G_xk' * inv(G_xk * Pm_k * G_xk' + R_k);

    % corrector equations
    P(:, :, k) = Pm_k - K_k * G_xk * Pm_k;
    x(k, :) = xm_k + (K_k * innovation)';
    x(k, 3) = normalizeAngle(x(k, 3));
  else % no corrector equations
    P(:, :, k) = Pm_k;
    x(k, 1:3) = xm_kR;
    x(k, 4:end) = x(k - 1, 4:end);
    x(k, 3) = normalizeAngle(x(k, 3));
  end
end
