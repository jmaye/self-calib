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

% This function performs Extended Kalman Filter (EKF) localization for a robot
% with a laser range finder with known landmark positions and calibration
% parameters.

function [x P] = ekf_loc(x_0, P_0, l, Theta, u, r, b, t, Q, R)

% angle normalization between -pi and pi
anglemod = @(x) atan2(sin(x), cos(x));

% number of calibration parameters
numCalib = length(Theta);

% timesteps to evaluate
steps = rows(r);

% init filtered output state
x = zeros(steps, 3);
x(1, :) = x_0;

% init filtered output covariance
P = zeros(3, 3, steps);
P(:, :, 1) = P_0;

% Jacobian of motion model with respect to state variable
Hx = eye(3, 3);

% Jacobian of motion model with respect to noise variable
Hw = zeros(3, 2);

% Jacobian of observation model with respect to state variable
Gx = zeros(2, 3);

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
  P(:, :, k) = Hx * P(:, :, k - 1) * Hx' + Hw * Q * Hw';
  x(k, 1) = x(k - 1, 1) + T * ctm1 * u(k, 1);
  x(k, 2) = x(k - 1, 2) + T * stm1 * u(k, 1);
  x(k, 3) = anglemod(x(k - 1, 3) + T * u(k, 2));

  % number of observations
  L_k = nnz(r(k, :) > 0);

  % loop over observations if any
  if L_k > 0
    % matrices allocation
    y_k = zeros(2 * L_k, 1);
    g_k = zeros(2 * L_k, 1);
    R_k = zeros(2 * L_k, 2 * L_k);
    G_xk = zeros(2 * L_k, 3);
    innovation = zeros(2 * L_k, 1);

    % some pre-computations
    st1 = sin(x(k, 3));
    ct1 = cos(x(k, 3));

    idx = 1;
    for j = 1:rows(l)
      if r(k, j) > 0
        % some pre-computations
        if numCalib < 3
          dct = Theta(1) * ct1;
          dst = Theta(1) * st1;
          aa = l(j, 1) - x(k, 1) - dct;
          bb = l(j, 2) - x(k, 2) - dst;
        else
          dxct = Theta(1) * ct1;
          dxst = Theta(1) * st1;
          dyct = Theta(2) * ct1;
          dyst = Theta(2) * st1;
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
          g(2) = anglemod(atan2(bb, aa) - x(k, 3) - Theta(3));
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
        G_xk(idx:idx + 1, :) = Gx;

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
  end
end
