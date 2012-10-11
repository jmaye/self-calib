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
% known poses and landmark positions in an iterative fashion using mutual
% information.

function [est Sigma batchIdx] = ls_calib_it(x, l, Theta_hat, r, b, R, ...
  maxIter, optTol, batchSize, miTol, rankTol)

% default values
if nargin < 7
  maxIter = 100;
end
if nargin < 8
  optTol = 1e-6;
end
if nargin < 9
  batchSize = 100;
end
if nargin < 10
  miTol = 0.5;
end

% number of calibration parameters
numCalib = length(Theta_hat);

% angle normalization between -pi and pi
anglemod = @(x) atan2(sin(x), cos(x));

% batch indices saved for optimization
batchIdx = [];

% number of batches used
numBatches = 0;

% record for covariances from previous iterations
sigmaRecord = [];

% record for covariances determinatn from previous iterations
sigmaDetRecord = [];

% timesteps
timesteps = rows(x);

% number of landmarks
nl = rows(l);

% init estimate
est = Theta_hat;

% transformed covariance matrix of observation model
N = R;

% Cholesky factor of the inverted transformed observation covariance
invNChol = sqrt(diag(1 ./ diag(R)));

% Jacobian of observation model with respect to calibration parameters
Gt = zeros(2, numCalib);

% loop over dataset
for i = 1:timesteps
  % temporarily add this data point into our data pool
  batchIdx = [batchIdx; i];

  % process batch
  if mod(i, batchSize) == 0
    % increment number of batches
    numBatches = numBatches + 1;

    % debugging print out
    i

    % number of lrf observations
    numObs = 0;
    for j = 1:numBatches
      % start of batch, don't count first point
      batchStart = (j - 1) * batchSize + 2;

      % increment observations
      numObs = numObs +...
        nnz(r(batchIdx(batchStart:batchStart + batchSize - 2), :) > 0) * 2;
    end

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

    % init temporary estimate
    est_temp = est;

    % non-linear least squares
    oldRes = 0;
    for s = 1:maxIter
      % update Jacobian and error term
      row = 1;
      nzcount = 1;
      for j = 1:numBatches
        % start of batch, don't count first point
        batchStart = (j - 1) * batchSize + 2;

        % process batch
        for k = batchIdx(batchStart):batchIdx(batchStart + batchSize - 2)
          % some pre-computations
          st = sin(x(k, 3));
          ct = cos(x(k, 3));

          % loop over the observations
          for ll = 1:nl
            if r(k, ll) > 0
              % some pre-computations
              if numCalib < 3
                aa = l(ll, 1) - x(k, 1) - est_temp(1) * ct;
                bb = l(ll, 2) - x(k, 2) - est_temp(1) * st;
              else
                aa = l(ll, 1) - x(k, 1) - est_temp(1) * ct + est_temp(2) * st;
                bb = l(ll, 2) - x(k, 2) - est_temp(1) * st - est_temp(2) * ct;
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
              e(row) = invNChol(1, 1) * (r(k, ll) - temp2);
              if numCalib < 3
                e(row + 1) = invNChol(2, 2) *...
                  anglemod(b(k, ll) - (atan2(bb, aa) - x(k, 3)));
              else
                e(row + 1) = invNChol(2, 2) *...
                  anglemod(b(k, ll) - (atan2(bb, aa) - x(k, 3) - est_temp(3)));
              end
              row = row + 2;
            end
          end
        end
      end

      H = sparse(ii, jj, ss, numObs, numCalib, nzmax);
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

      % update estimate
      if nargin < 11
        est_temp = est_temp + G * spqr_solve(H * G, e);
      else
        est_temp = est_temp + G * spqr_solve(H * G, e, struct('tol', rankTol));
      end
      if numCalib == 3
        est_temp(3) = anglemod(est_temp(3));
      end
    end

    % compute covariance
    Sigma = computeCov(H, e, numCalib);

    % compute covariance determinant
    sigmaDet = det(Sigma);

    % check if we need this batch
    if numBatches == 1 % first batch always taken
      est = est_temp;
      sigmaRecord = Sigma;
      sigmaDetRecord = sigmaDet;
    else
      % compute mutual information
      mi = 0.5 * log2(sigmaDetRecord / sigmaDet);

      % add batch if needed
      if mi > miTol
        est = est_temp;
        sigmaRecord = Sigma;
        sigmaDetRecord = sigmaDet;
      else
        % kick out the batch
        batchIdx(end - batchSize + 1:end) = [];
        numBatches = numBatches - 1;
      end
    end
  end
end

Sigma = sigmaRecord;
