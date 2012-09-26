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

% This function performs least squares SLAM for a robot with a laser
% range finder with known calibration parameters.

function [x_est l_est] = ls_slam(Theta, x_hat, l_hat, u, r, b, t, Q, R, ...
  maxIter, optTol, rankTol)

% default values
if nargin < 10
  maxIter = 100;
end
if nargin < 11
  optTol = 1e-6;
end

% timesteps to evaluate
steps = rows(r);

% number of landmarks
nl = rows(l_hat);

% number of state variables
ns = rows(x_hat);

% number of calibration parameters
numCalib = length(Theta);

% angle normalization between -pi and pi
anglemod = @(x) atan2(sin(x), cos(x));

% number of lrf observations
numObs = nnz(r(2:end, :) > 0) * 2;

% number of variables to estimate
numVar = ns * 3 + nl * 2;

% number of non-zero entries in the Jacobian
nzmax = (steps - 1) * 8 + numObs / 2 * 10;

% Jacobian initialization
ii = zeros(nzmax, 1);
jj = zeros(nzmax, 1);
ss = zeros(nzmax, 1);

% error term
e = zeros((ns - 1) * 3 + numObs, 1);

% transformed covariance matrix of observation model
N = R;

% Cholesky factor of the inverted transformed observation covariance
invNChol = sqrt(diag(1 ./ diag(R)));

% Jacobian of motion model with respect to state variable
Hx = eye(3, 3);

% Jacobian of motion model with respect to noise
Hw = zeros(3, 2);

% Jacobian of observation model with respect to state variables
Gx = zeros(2, 3);

% Jacobian of observation model with respect to landmark variables
Gl = zeros(2, 2);

% non-linear least squares
oldRes = 0;
x_est = x_hat;
l_est = l_hat;
for s = 1:maxIter
  % print out iteration number
  s

  % update Jacobian and error term
  row = 1;
  col = 1;
  nzcount = 1;

  for i = 2:steps
    % some pre-computations
    stm1 = sin(x_est(i - 1, 3));
    ctm1 = cos(x_est(i - 1, 3));

    % update Jacobian of the motion model with respect to noise
    Hw(1, 1) = (t(i) - t(i - 1)) * ctm1;
    Hw(2, 1) = (t(i) - t(i - 1)) * stm1;
    Hw(3, 2) = (t(i) - t(i - 1));

    % transform odometry covariance
    W = Hw * Q * Hw';

    % Cholesky factor of the inverted transformed motion covariance
    invWChol = sqrt(diag(1 ./ diag(W)));

    % update Jacobian of motion model with respect to state variable
    Hx(1, 3) = -(t(i) - t(i - 1)) * stm1 * u(i, 1);
    Hx(2, 3) = (t(i) - t(i - 1)) * ctm1 * u(i, 1);

    % update sparse matrix filling
    ii(nzcount) = row;
    jj(nzcount) = col;
    ss(nzcount) = -invWChol(1, 1);
    nzcount = nzcount + 1;
    ii(nzcount) = row;
    jj(nzcount) = col + 2;
    ss(nzcount) = -invWChol(1, 1) * Hx(1, 3);
    nzcount = nzcount + 1;
    ii(nzcount) = row + 1;
    jj(nzcount) = col + 1;
    ss(nzcount) = -invWChol(2, 2);
    nzcount = nzcount + 1;
    ii(nzcount) = row + 1;
    jj(nzcount) = col + 2;
    ss(nzcount) = -invWChol(2, 2) * Hx(2, 3);
    nzcount = nzcount + 1;
    ii(nzcount) = row + 2;
    jj(nzcount) = col + 2;
    ss(nzcount) = -invWChol(3, 3);
    nzcount = nzcount + 1;
    ii(nzcount) = row;
    jj(nzcount) = col + 3;
    ss(nzcount) = invWChol(1, 1);
    nzcount = nzcount + 1;
    ii(nzcount) = row + 1;
    jj(nzcount) = col + 4;
    ss(nzcount) = invWChol(2, 2);
    nzcount = nzcount + 1;
    ii(nzcount) = row + 2;
    jj(nzcount) = col + 5;
    ss(nzcount) = invWChol(3, 3);
    nzcount = nzcount + 1;

    % update error term
    e(row) = invWChol(1, 1) * (x_est(i, 1) -...
      (x_est(i - 1, 1) + (t(i) - t(i - 1)) * ctm1 * u(i, 1)));
    e(row + 1) = invWChol(2, 2) *...
      (x_est(i, 2) - (x_est(i - 1, 2) + (t(i) - t(i - 1)) * stm1 * u(i, 1)));
    e(row + 2) = anglemod(x_est(i, 3) -...
      (x_est(i - 1, 3) + (t(i) - t(i - 1)) * u(i, 2)));
    e(row + 2) = invWChol(3, 3) * e(row + 2);

    % update row/col counter
    row = row + 3;
    col = col + 3;

    % some pre-computations
    st1 = sin(x_est(i, 3));
    ct1 = cos(x_est(i, 3));

    % loop over the observations
    for j = 1:nl
      if r(i, j) > 0
        % some pre-computations
        if numCalib < 3
          dct = Theta(1) * ct1;
          dst = Theta(1) * st1;
          aa = l_est(j, 1) - x_est(i, 1) - dct;
          bb = l_est(j, 2) - x_est(i, 2) - dst;
        else
          dxct = Theta(1) * ct1;
          dxst = Theta(1) * st1;
          dyct = Theta(2) * ct1;
          dyst = Theta(2) * st1;
          aa = l_est(j, 1) - x_est(i, 1) - dxct + dyst;
          bb = l_est(j, 2) - x_est(i, 2) - dxst - dyct;
        end
        temp1 = aa^2 + bb^2;
        temp2 = sqrt(temp1);

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

        % update sparse matrix filling
        ii(nzcount) = row;
        jj(nzcount) = col;
        ss(nzcount) = -invNChol(1, 1) * Gx(1, 1);
        nzcount = nzcount + 1;
        ii(nzcount) = row;
        jj(nzcount) = col + 1;
        ss(nzcount) = -invNChol(1, 1) * Gx(1, 2);
        nzcount = nzcount + 1;
        ii(nzcount) = row;
        jj(nzcount) = col + 2;
        ss(nzcount) = -invNChol(1, 1) * Gx(1, 3);
        nzcount = nzcount + 1;
        ii(nzcount) = row + 1;
        jj(nzcount) = col;
        ss(nzcount) = -invNChol(2, 2) * Gx(2, 1);
        nzcount = nzcount + 1;
        ii(nzcount) = row + 1;
        jj(nzcount) = col + 1;
        ss(nzcount) = -invNChol(2, 2) * Gx(2, 2);
        nzcount = nzcount + 1;
        ii(nzcount) = row + 1;
        jj(nzcount) = col + 2;
        ss(nzcount) = -invNChol(2, 2) * Gx(2, 3);
        nzcount = nzcount + 1;

        % update Jacobian of observation model with respect to landmark variable
        Gl(1, 1) = aa / temp2;
        Gl(1, 2) = bb / temp2;
        Gl(2, 1) = -bb / temp1;
        Gl(2, 2) = aa / temp1;

        % update sparse matrix filling
        temp3 = ns * 3 + 1 + (j - 1) * 2;
        ii(nzcount) = row;
        jj(nzcount) = temp3;
        ss(nzcount) = -invNChol(1, 1) * Gl(1, 1);
        nzcount = nzcount + 1;
        ii(nzcount) = row;
        jj(nzcount) = temp3 + 1;
        ss(nzcount) = -invNChol(1, 1) * Gl(1, 2);
        nzcount = nzcount + 1;
        ii(nzcount) = row + 1;
        jj(nzcount) = temp3;
        ss(nzcount) = -invNChol(2, 2) * Gl(2, 1);
        nzcount = nzcount + 1;
        ii(nzcount) = row + 1;
        jj(nzcount) = temp3 + 1;
        ss(nzcount) = -invNChol(2, 2) * Gl(2, 2);
        nzcount = nzcount + 1;

        % update error term
        e(row) = invNChol(1, 1) * (r(i, j) - temp2);
        if numCalib < 3
          e(row + 1) = invNChol(2, 2) *...
            anglemod(b(i, j) - (atan2(bb, aa) - x_est(i, 3)));
        else
          e(row + 1) = invNChol(2, 2) *...
            anglemod(b(i, j) - (atan2(bb, aa) - x_est(i, 3) - Theta(3)));
        end
        row = row + 2;
      end
    end
  end
  H = sparse(ii, jj, ss, (ns - 1) * 3 + numObs, numVar, nzmax);
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
  if nargin < 12
    update = G * spqr_solve(H * G, -e);
  else
    update = G * spqr_solve(H * G, -e, struct('tol', rankTol));
  end
  x_est = x_est + [update(1:3:ns * 3) update(2:3:ns * 3) update(3:3:ns * 3)];
  x_est(:, 3) = anglemod(x_est(:, 3));
  l_est = l_est + [update(ns * 3 + 1:2:end) update(ns * 3 + 2:2:end)];
end
