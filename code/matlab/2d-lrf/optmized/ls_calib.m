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

% This function performs least squares calibration for a laser range finder with
% known poses and landmark positions.

function [est Sigma] = ls_calib(x, l, Theta_hat, r, b, R, maxIter, optTol, ...
  rankTol)

% default values
if nargin < 7
  maxIter = 100;
end
if nargin < 8
  optTol = 1e-6;
end

% number of calibration parameters
numCalib = length(Theta_hat);

% angle normalization between -pi and pi
anglemod = @(x) atan2(sin(x), cos(x));

% number of lrf observations
numObs = nnz(r(2:end, :)) * 2;

% number of non-zero entries in the Jacobian
if numCalib < 3
  nzmax = numObs;
else
  nzmax = numObs / 2 * 5;
end

% Jacobian initialization
ii = zeros(nzmax, 1);
jj = zeros(nzmax, 1);
ss = zeros(nzmax, 1);

% error term
e = zeros(numObs, 1);

% transformed covariance matrix of observation model
N = R;

% Cholesky factor of the inverted transformed observation covariance
invNChol = sqrt(diag(1 ./ diag(R)));

% timesteps to evaluate
steps = rows(r);

% number of landmarks
nl = rows(l);

% Jacobian of observation model with respect to calibration parameters
Gt = zeros(2, numCalib);

% non-linear least squares
oldRes = 0;
est = Theta_hat;
for s = 1:maxIter
  % print out iteration number and current estimate
  s
  est

  % update Jacobian and error term
  row = 1;
  nzcount = 1;
  for i = 2:steps
    % some pre-computations
    st = sin(x(i, 3));
    ct = cos(x(i, 3));

    % loop over the observations
    for j = 1:nl
      if r(i, j) > 0
        % some pre-computations
        if numCalib < 3
          aa = l(j, 1) - x(i, 1) - est(1) * ct;
          bb = l(j, 2) - x(i, 2) - est(1) * st;
        else
          aa = l(j, 1) - x(i, 1) - est(1) * ct + est(2) * st;
          bb = l(j, 2) - x(i, 2) - est(1) * st - est(2) * ct;
        end
        temp1 = aa^2 + bb^2;
        temp2 = sqrt(temp1);

        % update Jacobian
        if numCalib < 3
          Gt(1, 1) = -(aa * ct + bb * st) / temp2;
          Gt(2, 1) = (-aa * st + bb * ct) / temp1;
        else
          Gt(1, 1) = -(aa * ct + bb * st) / temp2;
          Gt(1, 2) = (aa * st - bb * ct) / temp2;
          Gt(2, 1) = (-aa * st + bb * ct) / temp1;
          Gt(2, 2) = -(aa * ct + bb * st) / temp1;
          Gt(2, 3) = -1;
        end

        % update sparse matrix filling
        if numCalib < 3
          ii(nzcount) = row;
          jj(nzcount) = 1;
          ss(nzcount) = invNChol(1, 1) * Gt(1, 1);
          nzcount = nzcount + 1;
          ii(nzcount) = row + 1;
          jj(nzcount) = 1;
          ss(nzcount) = invNChol(2, 2) * Gt(2, 1);
          nzcount = nzcount + 1;
        else
          ii(nzcount) = row;
          jj(nzcount) = 1;
          ss(nzcount) = invNChol(1, 1) * Gt(1, 1);
          nzcount = nzcount + 1;
          ii(nzcount) = row;
          jj(nzcount) = 2;
          ss(nzcount) = invNChol(1, 1) * Gt(1, 2);
          nzcount = nzcount + 1;
          ii(nzcount) = row + 1;
          jj(nzcount) = 1;
          ss(nzcount) = invNChol(2, 2) * Gt(2, 1);
          nzcount = nzcount + 1;
          ii(nzcount) = row + 1;
          jj(nzcount) = 2;
          ss(nzcount) = invNChol(2, 2) * Gt(2, 2);
          nzcount = nzcount + 1;
          ii(nzcount) = row + 1;
          jj(nzcount) = 3;
          ss(nzcount) = invNChol(2, 2) * Gt(2, 3);
          nzcount = nzcount + 1;
        end

        % update error term
        e(row) = invNChol(1, 1) * (r(i, j) - temp2);
        if numCalib < 3
          e(row + 1) = invNChol(2, 2) *...
            anglemod(b(i, j) - (atan2(bb, aa) - x(i, 3)));
        else
          e(row + 1) = invNChol(2, 2) *...
            anglemod(b(i, j) - (atan2(bb, aa) - x(i, 3) - est(3)));
        end
        row = row + 2;
      end
    end
  end
  H = sparse(ii, jj, ss, numObs, numCalib, nzmax);

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
  if nargin < 9
    est = est + spqr_solve(H, e);
  else
    est = est + spqr_solve(H, e, struct('tol', rankTol));
  end
  if numCalib == 3
    est(3) = anglemod(est(3));
  end
end

% compute covariance
Sigma = computeCov(H, e, numCalib);
