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

% This function performs Extended Kalman Filter (EKF) SLAM for a robot with a
% laser range finder and known calibration parameters.

function [x P] = ekf_slam(x_0, P_0, Theta, u, r, b, t, Q, R)

% angle normalization between -pi and pi
anglemod = @(x) atan2(sin(x), cos(x));

% number of calibration parameters
numCalib = length(Theta);

% timesteps to evaluate
steps = rows(r);

% number of landmarks
nl = cols(r);

% pointer to robot
rr = 1:3;

% pointer to map
mm = 4:4 + nl * 2 - 1;

% init filtered output state
x = zeros(steps, 3 + nl * 2);
x(1, rr) = x_0;

% init filtered output covariance
P = zeros(3 + nl * 2, 3 + nl * 2, steps);
P(rr, rr, 1) = P_0;

% init landmarks
Fx = zeros(2, 3);
Fy = zeros(2, 2);
for i = 1:nl
  Fx(1, 1) = 1;
  Fx(2, 2) = 1;
  if numCalib < 3
    x(1, 4 + (i - 1) * 2) = x(1, 1) + Theta(1) * cos(x(1, 3)) + r(2, i) *...
      cos(b(2, i) + x(1, 3));
    x(1, 4 + (i - 1) * 2 + 1) = x(1, 2) + Theta(1) * sin(x(1, 3)) + r(2, i) *...
      sin(b(2, i) + x(1, 3));
    Fx(1, 3) = -Theta(1) * sin(x(1, 3)) - r(2, i) * sin(x(1, 3) + b(2, i));
    Fx(2, 3) = Theta(1) * cos(x(1, 3)) + r(2, i) * cos(x(1, 3) + b(2, i));
    Fy(1, 1) = cos(x(1, 3) + b(2, i));
    Fy(1, 2) = -r(2, i) * sin(b(2, i) + x(1, 3));
    Fy(2, 1) = sin(b(2, i) + x(1, 3));
    Fy(2, 2) = -r(2, i) * cos(b(2, i) + x(1, 3));
  else
    x(1, 4 + (i - 1) * 2)  = x(1, 1) + Theta(1) * cos(x(1, 3)) - Theta(2) *...
      sin(x(1, 3)) + r(2, i) * cos(b(2, i) + Theta(3) + x(1, 3));
    x(1, 4 + (i - 1) * 2 + 1) = x(1, 2) + Theta(1) * sin(x(1, 3)) +...
      Theta(2) * cos(x(1, 3)) + r(2, i) * sin(b(2, i) + Theta(3) + x(1, 3));
    Fx(1, 3) = -Theta(1) * sin(x(1, 3)) - Theta(2) * cos(x(1, 3)) -...
      r(2, i) * sin(x(1, 3) + b(2, i) + Theta(3));
    Fx(2, 3) = Theta(1) * cos(x(1, 3)) - Theta(2) * sin(x(1, 3)) +...
      r(2, i) * cos(x(1, 3) + b(2, i) + Theta(3));
    Fy(1, 1) = cos(x(1, 3) + b(2, i) + Theta(3));
    Fy(1, 2) = -r(2, i) * sin(b(2, i) + x(1, 3) + Theta(3));
    Fy(2, 1) = sin(b(2, i) + x(1, 3) + Theta(3));
    Fy(2, 2) = -r(2, i) * cos(b(2, i) + x(1, 3) + Theta(3));
  end
  P(4 + (i - 1) * 2:4 + (i - 1) * 2 + 1, ...
    4 + (i - 1) * 2:4 + (i - 1) * 2 + 1, 1) =...
    Fx * P(rr, rr, 1) * Fx' + Fy * R * Fy';
end

% Jacobian of motion model with respect to state variable
Hx = eye(3, 3);

% Jacobian of motion model with respect to noise variable
Hw = zeros(3, 2);

% Jacobian of observation model with respect to state variable
Gx = zeros(2, 3);

% Jacobian of observation model with respect to landmark variable
Gl = zeros(2, 2);

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
  P(rr, mm, k) = Hx * P(rr, mm, k - 1);
  P(mm, rr, k) = P(rr, mm, k)';
  P(mm, mm, k) = P(mm, mm, k - 1);
  x(k, 1) = x(k - 1, 1) + T * ctm1 * u(k, 1);
  x(k, 2) = x(k - 1, 2) + T * stm1 * u(k, 1);
  x(k, 3) = anglemod(x(k - 1, 3) + T * u(k, 2));
  x(k, mm) = x(k - 1, mm);

  % number of observations
  L_k = nnz(r(k, :) > 0);

  % loop over observations if any
  if L_k > 0
    % matrices allocation
    y_k = zeros(2 * L_k, 1);
    g_k = zeros(2 * L_k, 1);
    R_k = zeros(2 * L_k, 2 * L_k);
    G_xk = zeros(2 * L_k, 3 + nl * 2);
    innovation = zeros(2 * L_k, 1);

    % some pre-computations
    st1 = sin(x(k, 3));
    ct1 = cos(x(k, 3));

    idx = 1;
    for j = 1:nl
      if r(k, j) > 0
        % some pre-computations
        if numCalib < 3
          dct = Theta(1) * ct1;
          dst = Theta(1) * st1;
          aa = x(k, 4 + (j - 1) * 2) - x(k, 1) - dct;
          bb = x(k, 4 + (j - 1) * 2 + 1) - x(k, 2) - dst;
        else
          dxct = Theta(1) * ct1;
          dxst = Theta(1) * st1;
          dyct = Theta(2) * ct1;
          dyst = Theta(2) * st1;
          aa = x(k, 4 + (j - 1) * 2) - x(k, 1) - dxct + dyst;
          bb = x(k, 4 + (j - 1) * 2 + 1) - x(k, 2) - dxst - dyct;
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
        G_xk(idx:idx + 1, rr) = Gx;

        % update Jacobian of observation model with respect to landmark variable
        Gl(1, 1) = aa / temp2;
        Gl(1, 2) = bb / temp2;
        Gl(2, 1) = -bb / temp1;
        Gl(2, 2) = aa / temp1;
        G_xk(idx:idx + 1, 4 + (j - 1) * 2:4 + (j - 1) * 2 + 1) = Gl;

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
