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
% based on Tim Barfoot's example. This is the definitive version of the
% iterative algorithm we need, but with known landmark's poses and robot's poses

% angle normalization between -pi and pi
normalizeAngle = @(x) atan2(sin(x), cos(x));

% observation model covariance
R = [r_var, 0; 0, b_var];

% batch size
batchSize = 100;

% batch indices saved for optimization
batchIndices = [];

% mutual information threshold in bits
miThreshold = 0.5;

% Estimates
d_est = d + randn * 0.1

varianceRecord = [];

% loop over dataset
%for i = 2:length(x_odom)
for i = 2:1000
  % process batch
  if mod(i, batchSize) == 0
    i

    % initial estimates
    d_est_temp = d_est;

    % temporarily add this batch
    batchIndices = [batchIndices; i];

    % count how many laser range measurements to be considered
    obsCount = 0;
    for j = 1:length(batchIndices)
      for k = batchIndices(j) - batchSize + 2:batchIndices(j)
        for m = 1:rows(l)
          if r(k, m) > 0
            obsCount = obsCount + 1;
          end
        end
      end
    end

    % jacobian matrix allocation
    H = sparse(obsCount * 2, 1);

    % error terms allocation
    e = zeros(obsCount * 2, 1);

    % optimization using old and new dataset
    tol = 1e-6;
    maxNumIter = 200;
    oldll = 0;
    for s = 1:maxNumIter
      % emit iteration number
      s;

      % build matrices
      row = 1;
      tic;
      n = 1;
      for j = 1:length(batchIndices)
        n = n + 1;
        for k = batchIndices(j) - batchSize + 2:batchIndices(j)
          % some pre-computations
          st1 = sin(th_true(n));
          ct1 = cos(th_true(n));

          % covariance matrix of observation model
          N = R;

          % inverted covariance matrix
          invN = diag(1 ./ diag(R));

          % Cholesky factor of covariance matrix
          invN_sqrt = chol(invN);

          % loop over the observations
          for m = 1:rows(l)
            if r(k, m) > 0
              % some pre-computations
              a1 = l(m, 1) - x_true(n) - d_est_temp * ct1;
              b1 = l(m, 2) - y_true(n) - d_est_temp * st1;
              temp1 = a1^2 + b1^2;
              temp2 = sqrt(temp1);

              % jacobian of observation model with respect to calibration param.
              G_d = zeros(2, 1);
              G_d(1, 1) = -(a1 * ct1 + b1 * st1) / temp2;
              G_d(2, 1) = (-a1 * st1 + b1 * ct1) / temp1;
              G_d_cov = -invN_sqrt * G_d;

              % setting everything into H and e
              H(row, 1) = G_d_cov(1, 1);
              H(row + 1, 1) = G_d_cov(2, 1);
              e(row) = r(k, m) - temp2;
              e(row + 1) = b(k, m) - (atan2(b1, a1) - th_true(n));
              e(row + 1) = normalizeAngle(e(row + 1));
              e(row:row + 1) = invN_sqrt * e(row:row + 1);
              row = row + 2;
            end
          end
          n = n + 1;
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

      ll;

      % compute update
      tic;
      dx = spqr_solve(H, -e, struct('tol', 5));
      toc;

      % update estimate
      d_update = dx(end);
      d_est_temp = d_est_temp + d_update;
    end

    % variance on the calibration parameter
    [C1, R1, P1] = spqr(H, -e, struct('permutation', 'matrix', 'econ', ...
      cols(H)));
    R1 = P1 * R1 * P1';
    variance = 1 / R1(end, end)^2;

    % check if we need this batch
    if length(batchIndices) == 1 % first batch always taken
      d_est = d_est_temp;
      varianceRecord = variance;
    else
      % compute mutual information
      mi = 0.5 * log2(varianceRecord / variance);

      % add batch if needed
      if mi > miThreshold
        d_est = d_est_temp;
        varianceRecord = variance;
      else
        % kick out the batch
        batchIndices(end) = [];
      end
    end
  end
end
