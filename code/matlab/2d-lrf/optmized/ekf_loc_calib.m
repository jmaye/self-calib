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

% This function performs Extended Kalman Filter (EKF) localization and
% calibration for a robot with a laser range finder with known landmark
% positions.

function [x P] = ekf_loc_calib(x_0, P_0, Theta_0, Sigma_0, l, u, r, b, t, Q, R)

% angle normalization between -pi and pi
anglemod = @(x) atan2(sin(x), cos(x));

% number of calibration parameters
numCalib = length(Theta_0);

% timesteps to evaluate
steps = rows(r);

% pointer to robot
rr = 1:3;

% pointer to calibration parameter
tt = 4:4 + numCalib - 1;

% init filtered output state
x = zeros(steps, 3 + numCalib);
x(1, rr) = x_0;
x(1, tt) = Theta_0;

% init filtered output covariance
P = zeros(3 + numCalib, 3 + numCalib, steps);
P(rr, rr, 1) = P_0;
P(tt, tt, 1) = Sigma_0;

% Jacobian of motion model with respect to state variable
Hx = eye(3, 3);

% Jacobian of motion model with respect to noise variable
Hw = zeros(3, 2);

% Jacobian of observation model with respect to state variable
Gx = zeros(2, 3);

% Jacobian of observation model with respect to calibration parameter
Gt = zeros(2, numCalib);

% observation
y = zeros(2, 1);

% result of observation model
g = zeros(2, 1);

% perform filtering
for k = 2:steps
  % some pre-computations
  stm1 = sin(x(k - 1, 3));
  ctm1 = cos(x(k - 1, 3));

  % update Jacobian of motion model with respect to state variables
  Hx(1, 3) = -(t(k) - t(k - 1)) * stm1 * u(k, 1);
  Hx(2, 3) = (t(k) - t(k - 1)) * ctm1 * u(k, 1);

  % update Jacobian of motion model with respect to noise
  Hw(1, 1) = (t(k) - t(k - 1)) * ctm1;
  Hw(2, 1) = (t(k) - t(k - 1)) * stm1;
  Hw(3, 2) = (t(k) - t(k - 1));

  % predictor equations
  T = t(k) - t(k - 1);
  P(rr, rr, k) = Hx * P(rr, rr, k - 1) * Hx' + Hw * Q * Hw';
  P(rr, tt, k) = Hx * P(rr, tt, k - 1);
  P(tt, rr, k) = P(rr, tt, k)';
  P(tt, tt, k) = P(tt, tt, k - 1);
  x(k, 1) = x(k - 1, 1) + T * ctm1 * u(k, 1);
  x(k, 2) = x(k - 1, 2) + T * stm1 * u(k, 1);
  x(k, 3) = anglemod(x(k - 1, 3) + T * u(k, 2));
  x(k, tt) = x(k - 1, tt);

  % number of observations
  L_k = nnz(r(k, :));

  % loop over observations if any
  if L_k > 0
    % matrices allocation
    y_k = zeros(2 * L_k, 1);
    g_k = zeros(2 * L_k, 1);
    R_k = zeros(2 * L_k, 2 * L_k);
    G_xk = zeros(2 * L_k, 3 + numCalib);
    innovation = zeros(2 * L_k, 1);

    % some pre-computations
    st1 = sin(x(k, 3));
    ct1 = cos(x(k, 3));

    idx = 1;
    for j = 1:rows(l)
      if r(k, j) > 0
        % some pre-computations
        if numCalib < 3
          dct = x(k, tt) * ct1;
          dst = x(k, tt) * st1;
          aa = l(j, 1) - x(k, 1) - dct;
          bb = l(j, 2) - x(k, 2) - dst;
        else
          dxct = x(k, 4) * ct1;
          dxst = x(k, 4) * st1;
          dyct = x(k, 5) * ct1;
          dyst = x(k, 5) * st1;
          aa = l(j, 1) - x(k, 1) - dxct + dyst;
          bb = l(j, 2) - x(k, 2) - dxst - dyct;
        end
        temp1 = aa^2 + bb^2;
        temp2 = sqrt(temp1);

        % update observation vector
        y(1) = r(k, j);
        y(2) = b(k, j);
        y_k(idx:idx + 1) = y;

        % update prediction vector
        g(1) = temp2;
        if numCalib < 3
          g(2) = anglemod(atan2(bb, aa) - x(k, 3));
        else
          g(2) = anglemod(atan2(bb, aa) - x(k, 3) - x(k, 6));
        end
        g_k(idx:idx + 1) = g;

        % update Jacobian of observation model with respect to state variable
        Gx(1, 1) = -aa / temp2;
        Gx(1, 2) = -bb / temp2;
        Gx(2, 1) = bb / temp1;
        Gx(2, 2) = -aa / temp1;
        if numCalib < 3
          Gx(1, 3) = (aa * dst - bb * dct) / temp2;
          Gx(2, 3) = -(aa * dct + bb * dst) / temp1 - 1;
        else
          Gx(1, 3) = (aa * (dxst + dyct) + bb * (-dxct + dyst)) / temp2;
          Gx(2, 3) = (aa * (-dxct + dyst) - bb * (dxst + dyct)) / temp1 - 1;
        end
        G_xk(idx:idx + 1, rr) = Gx;

        % update Jacobian of observation model with respect to calib. variable
        if numCalib < 3
          Gt(1, 1) = -(aa * ct1 + bb * st1) / temp2;
          Gt(2, 1) = (-aa * st1 + bb * ct1) / temp1;
        else
          Gt(1, 1) = -(aa * ct1 + bb * st1) / temp2;
          Gt(1, 2) = (aa * st1 - bb * ct1) / temp2;
          Gt(2, 1) = (-aa * st1 + bb * ct1) / temp1;
          Gt(2, 2) = -(aa * ct1 + bb * st1) / temp1;
          Gt(2, 3) = -1;
        end
        G_xk(idx:idx + 1, tt) = Gt;

        % compute innovation
        innovation(idx:idx + 1) = y - g;
        innovation(idx + 1) = anglemod(innovation(idx + 1));

        % observation covariance
        R_k(idx:idx + 1, idx:idx + 1) = R;

        % update count
        idx = idx + 2;
      end
    end

    % Kalman gain
    K_k = P(:, :, k) * G_xk' * inv(G_xk * P(:, :, k) * G_xk' + R_k);

    % corrector equations
    P(:, :, k) = P(:, :, k) - K_k * G_xk * P(:, :, k);
    x(k, :) = x(k, :)' + K_k * innovation;
    x(k, 3) = anglemod(x(k, 3));
    if numCalib == 3
      x(k, 6) = anglemod(x(k, 6));
    end
  end
end
