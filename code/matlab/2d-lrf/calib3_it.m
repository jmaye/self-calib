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
% parameters based on Tim Barfoot's example. This is the definitive version of
% the iterative algorithm we need, but with known landmark's poses and robot's
% poses

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

% estimates
d_x_est = d_x + randn * 0.5
d_y_est = d_y + randn * 0.5
phi_est = phi + randn * 0.5

covarianceDetRecord = [];

% loop over dataset
%for i = 2:length(x_odom)
for i = 2:1000
  % process batch
  if mod(i, batchSize) == 0
    i

    % initial estimates
    d_x_est_temp = d_x_est;
    d_y_est_temp = d_y_est;
    phi_est_temp = phi_est;

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
    H = sparse(obsCount * 2, 3);

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
              a1 = l(m, 1) - x_true(n) - d_x_est_temp * ct1 + d_y_est_temp *...
                st1;
              b1 = l(m, 2) - y_true(n) - d_x_est_temp * st1 - d_y_est_temp *...
                ct1;
              temp1 = a1^2 + b1^2;
              temp2 = sqrt(temp1);

              % jacobian of observation model with respect to calibration param.
              G_d = zeros(2, 3);
              G_d(1, 1) = -(a1 * ct1 + b1 * st1) / temp2;
              G_d(1, 2) = (a1 * st1 - b1 * ct1) / temp2;
              G_d(2, 1) = (-a1 * st1 + b1 * ct1) / temp1;
              G_d(2, 2) = -(a1 * ct1 + b1 * st1) / temp1;
              G_d(2, 3) = -1;
              G_d_cov = -invN_sqrt * G_d;

              % setting everything into H and e
              H(row, 1) = G_d_cov(1, 1);
              H(row, 2) = G_d_cov(1, 2);
              H(row + 1, 1) = G_d_cov(2, 1);
              H(row + 1, 2) = G_d_cov(2, 2);
              H(row + 1, 3) = G_d_cov(2, 3);
              e(row) = r(k, m) - temp2;
              e(row + 1) = b(k, m) - (atan2(b1, a1) - th_true(n) -...
                phi_est_temp);
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
      d_x_est_temp = d_x_est_temp + dx(1);
      d_y_est_temp = d_y_est_temp + dx(2);
      phi_est_temp = phi_est_temp + dx(3);
      phi_est_temp = normalizeAngle(phi_est_temp);
    end

    % covariance of the calibration parameters
    [C1, R1, P1] = spqr(H, -e, struct('permutation', 'matrix', 'econ', ...
      cols(H)));
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
      cov_m = cov_k - 1;
      for m = k - 1:-1:cols(H) - 2
        temp1 = 0;
        temp2 = 0;
        cov_j = cov_m + 1;
        for j = m + 1:k
          if R1(m, j) ~= 0
            temp1 = temp1 + R1(m, j) * covariance(cov_j, cov_k);
          end
          cov_j = cov_j + 1;
        end
        cov_j = cov_k + 1;
        for j = k + 1:cols(H)
          if R1(m, j) ~= 0
            temp2 = temp2 + R1(m, j) * covariance(cov_k, cov_j);
          end
          cov_j = cov_j + 1;
        end
        covariance(cov_m, cov_k) = 1 / R1(m, m) * (-temp1 - temp2);
        covariance(cov_k, cov_m) = covariance(cov_m, cov_k);
        cov_m = cov_m - 1;
      end
      cov_k = cov_k - 1;
    end

    % check if we need this batch
    if length(batchIndices) == 1 % first batch always taken
      d_x_est = d_x_est_temp;
      d_y_est = d_y_est_temp;
      phi_est = phi_est_temp;
      covarianceDetRecord = det(covariance);
    else
      % compute covariance determinant
      covarianceDet = det(covariance);

      % compute mutual information
      mi = 0.5 * log2(covarianceDetRecord / covarianceDet);

      % add batch if needed
      if mi > miThreshold
        d_x_est = d_x_est_temp;
        d_y_est = d_y_est_temp;
        phi_est = phi_est_temp;
        covarianceDetRecord = covarianceDet;
      else
        % kick out the batch
        batchIndices(end) = [];
      end
    end
  end
end
