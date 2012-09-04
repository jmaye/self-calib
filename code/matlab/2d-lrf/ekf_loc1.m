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

% This script is performing localization based on Tim Barfoot's example using
% an Extended Kalman Filter (EKF). Landmark's positions and laser range finder
% offset are considered as given.

% motion model (without noise)
h = @(x, u, T) [(x(:, 1) + T .* cos(x(:, 3)) .* u(:, 1)), (x(:, 2) +...
  T .* sin(x(:, 3)) .* u(:, 1)), (x(:, 3) + T .* u(:, 2))];

% angle normalization between -pi and pi
normalizeAngle = @(x) atan2(sin(x), cos(x));

% odometry measurements
u = [v, om];

% number of steps to estimate
steps = length(x_true);

% odometry
x_odom = zeros(steps, 3);
x_odom(1, :) = [x_true(1), y_true(1), th_true(1)];

% filtered state
x = zeros(steps, 3);
x(1, :) = [x_true(1), y_true(1), th_true(1)];

% filtered covariance
P = zeros(3, 3, steps);
P(:, :, 1) = diag([1, 1, 0.1]);

% motion model covariance
Q = [v_var, 0; 0, om_var];

% observation model covariance
R = [r_var, 0; 0, b_var];

% matrices allocation
H_x = eye(cols(x), cols(x));
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
  Pm_k = H_x * P(:, :, k - 1) * H_x' + H_w * Q * H_w';
  xm_k = h(x(k - 1, :), u(k, :), t(k) - t(k - 1));
  x_odom(k, :) = h(x_odom(k - 1, :), u(k, :), t(k) - t(k - 1));

  % number of observations
  L_k = nnz(r(k, :));

  % loop over observations if any
  if L_k > 0
    % matrices allocation
    y_k = zeros(2 * L_k, 1);
    g_k = zeros(2 * L_k, 1);
    R_k = zeros(2 * L_k, 2 * L_k);
    G_xk = zeros(2 * L_k, 3);
    innovation = zeros(2 * L_k, 1);
    y = zeros(2, 1);
    g = zeros(2, 1);
    G_x = zeros(2, 3);

    % some pre-computations
    st1 = sin(xm_k(3));
    ct1 = cos(xm_k(3));

    idx = 1;
    for j = 1:rows(l)
      if r(k, j) > 0
        % some pre-computations
        a1 = l(j, 1) - xm_k(1) - d * ct1;
        b1 = l(j, 2) - xm_k(2) - d * st1;
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
        G_xk(idx:idx + 1, :) = G_x;

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
    x(k, :) = xm_k;
    x(k, 3) = normalizeAngle(x(k, 3));
  end
end
