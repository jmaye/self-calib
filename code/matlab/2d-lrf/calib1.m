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

% This script is performing self-calibration of a laser range finder offset
% based on Tim Barfoot's example. For the sake of simplicity, robot's poses
% and landmark' positions are not estimated here.

% angle normalization between -pi and pi
normalizeAngle = @(x) atan2(sin(x), cos(x));

% number of steps to estimate
steps = 100;

% initial estimate
d_est = d + randn * 0.5

% count the number of observation rows the matrices have
rowCount = 0;
for i = 1:steps
  for j = 1:rows(l)
    if r(i, j) > 0
      rowCount = rowCount + 1;
    end
  end
end

% jacobian matrix
H = sparse(rowCount * 2, 1);

% error terms
e = zeros(rowCount * 2, 1);

% observation model covariance
R = [r_var, 0; 0, b_var];

% covariance matrix of observation model
N = R;

% inverted covariance matrix
invN = diag(1 ./ diag(R));

% Cholesky factor of covariance matrix
invN_sqrt = chol(invN);

% non-linear least square
tol = 1e-9;
maxNumIter = 200;
oldll = 0;
for s = 1:maxNumIter
  s

  % build matrices
  row = 1;
  tic;
  for i = 1:steps
    % some pre-computations
    st1 = sin(th_true(i));
    ct1 = cos(th_true(i));

    % loop over the observations
    for j = 1:rows(l)
      if r(i, j) > 0
        % some pre-computations
        a1 = l(j, 1) - x_true(i) - d_est * ct1;
        b1 = l(j, 2) - y_true(i) - d_est * st1;
        temp1 = a1^2 + b1^2;
        temp2 = sqrt(temp1);

        % jacobian of observation model with respect to calibration param.
        G_d = zeros(2, 1);
        G_d(1, 1) = -(a1 * ct1 + b1 * st1) / temp2;
        G_d(2, 1) = (-a1 * st1 + b1 * ct1) / temp1;
        G_d_cov = invN_sqrt * G_d;

        % setting everything into H and e
        H(row, 1) = G_d_cov(1, 1);
        H(row + 1, 1) = G_d_cov(2, 1);
        e(row) = r(i, j) - temp2;
        e(row + 1) = b(i, j) - (atan2(b1, a1) - th_true(i));
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
  dx = spqr_solve(H, e, struct('tol', 5));
  toc;

  % update estimate
  d_est = d_est + dx

end

% variance on the calibration parameter
[C1, R1, P1] = spqr(H, -e, struct('permutation', 'matrix', 'econ', cols(H)));
R1 = P1 * R1 * P1';
variance = 1 / R1(end, end)^2;

disp('Final calibration');
d_est
variance
