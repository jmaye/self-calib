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

% This function performs least squares SLAM and calibration for a robot
% with a laser range finder in an iterative fashion using mutual information.

function [x_est l_est Theta_est Sigma batchIdx] =...
  ls_slam_calib_it_auto(x_hat, l_hat, Theta_hat, u, r, b, t, Q, R, maxIter, ...
  optTol, batchSize, miTol, rankGap)

% default values
if nargin < 10
  maxIter = 20;
end
if nargin < 11
  optTol = 1e-6;
end
if nargin < 12
  batchSize = 200;
end
if nargin < 13
  miTol = 1.5;
end
if nargin < 14
  rankGap = 0.02;
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
timesteps = rows(x_hat);

% number of landmarks
nl = rows(l_hat);

% init estimates
x_est = [];
l_est = l_hat;
Theta_est = Theta_hat;

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

    % number of state variables
    ns = batchSize * numBatches;

    % number of variables to estimate
    numVar = ns * 3 + nl * 2 + numCalib;

    % number of non-zero entries in the Jacobian
    if numCalib < 3
      nzmax = (ns - numBatches) * 8 + numObs / 2 * 12;
    else
      nzmax = (ns - numBatches) * 8 + numObs / 2 * 15;
    end

    % Jacobian initialization
    ii = zeros(nzmax, 1);
    jj = zeros(nzmax, 1);
    ss = zeros(nzmax, 1);

    % error term
    e = zeros((ns - numBatches) * 3 + numObs, 1);

    % init temporary estimates
    x_est_temp = [x_est; x_hat(batchIdx(end - batchSize  + 1:end), :)];
    l_est_temp = l_est;
    Theta_est_temp = Theta_est;

    % non-linear least squares
    oldRes = 0;
    for s = 1:maxIter
      % update Jacobian and error term
      row = 1;
      col = 1;
      nzcount = 1;
      n = 1;
      for j = 1:numBatches
        % start of batch, don't count first point
        batchStart = (j - 1) * batchSize + 2;
        n = n + 1;

        % process batch
        for k = batchIdx(batchStart):batchIdx(batchStart + batchSize - 2)
          % some pre-computations
          stm1 = sin(x_est_temp(n - 1, 3));
          ctm1 = cos(x_est_temp(n - 1, 3));

          % update Jacobian of the motion model with respect to noise
          Hw(1, 1) = (t(k) - t(k - 1)) * ctm1;
          Hw(2, 1) = (t(k) - t(k - 1)) * stm1;
          Hw(3, 2) = (t(k) - t(k - 1));

          % transform odometry covariance
          W = Hw * Q * Hw';

          % Cholesky factor of the inverted transformed motion covariance
          invWChol = sqrt(diag(1 ./ diag(W)));

          % update Jacobian of motion model with respect to state variable
          Hx(1, 3) = -(t(k) - t(k - 1)) * stm1 * u(k, 1);
          Hx(2, 3) = (t(k) - t(k - 1)) * ctm1 * u(k, 1);

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
          e(row) = invWChol(1, 1) * (x_est_temp(n, 1) -...
            (x_est_temp(n - 1, 1) + (t(k) - t(k - 1)) * ctm1 * u(k, 1)));
          e(row + 1) = invWChol(2, 2) *...
            (x_est_temp(n, 2) - (x_est_temp(n - 1, 2) + (t(k) - t(k - 1)) *...
            stm1 * u(k, 1)));
          e(row + 2) = invWChol(3, 3) * anglemod(x_est_temp(n, 3) -...
            (x_est_temp(n - 1, 3) + (t(k) - t(k - 1)) * u(k, 2)));

          % update row/col counter
          row = row + 3;
          col = col + 3;

          % some pre-computations
          st1 = sin(x_est_temp(n, 3));
          ct1 = cos(x_est_temp(n, 3));

          % loop over the observations
          for ll = 1:nl
            if r(k, ll) > 0
              % some pre-computations
              if numCalib < 3
                dct = Theta_est_temp(1) * ct1;
                dst = Theta_est_temp(1) * st1;
                aa = l_est_temp(ll, 1) - x_est_temp(n, 1) - dct;
                bb = l_est_temp(ll, 2) - x_est_temp(n, 2) - dst;
              else
                dxct = Theta_est_temp(1) * ct1;
                dxst = Theta_est_temp(1) * st1;
                dyct = Theta_est_temp(2) * ct1;
                dyst = Theta_est_temp(2) * st1;
                aa = l_est_temp(ll, 1) - x_est_temp(n, 1) - dxct + dyst;
                bb = l_est_temp(ll, 2) - x_est_temp(n, 2) - dxst - dyct;
              end
              temp1 = aa^2 + bb^2;
              temp2 = sqrt(temp1);

              % update Jacobian of observation model with respect to state var.
              Gx(1, 1) = -aa / temp2;
              Gx(1, 2) = -bb / temp2;
              Gx(2, 1) = bb / temp1;
              Gx(2, 2) = -aa / temp1;
              if numCalib < 3
                Gx(1, 3) = (aa * dst - bb * dct) / temp2;
                Gx(2, 3) = -(aa * dct + bb * dst) / temp1 - 1;
              else
                Gx(1, 3) = (aa * (dxst + dyct) + bb * (-dxct + dyst)) / temp2;
                Gx(2, 3) = (aa * (-dxct + dyst) -...
                  bb * (dxst + dyct)) / temp1 - 1;
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

              % update Jacobian of observation model with respect to land. var.
              Gl(1, 1) = aa / temp2;
              Gl(1, 2) = bb / temp2;
              Gl(2, 1) = -bb / temp1;
              Gl(2, 2) = aa / temp1;

              % update sparse matrix filling
              temp3 = ns * 3 + 1 + (ll - 1) * 2;
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

              % update Jacobian of observation model with respect to calib. var.
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

              % update sparse matrix filling
              if numCalib < 3
                ii(nzcount) = row;
                jj(nzcount) = numVar;
                ss(nzcount) = -invNChol(1, 1) * Gt(1, 1);
                nzcount = nzcount + 1;
                ii(nzcount) = row + 1;
                jj(nzcount) = numVar;
                ss(nzcount) = -invNChol(2, 2) * Gt(2, 1);
                nzcount = nzcount + 1;
              else
                ii(nzcount) = row;
                jj(nzcount) = numVar - 2;
                ss(nzcount) = -invNChol(1, 1) * Gt(1, 1);
                nzcount = nzcount + 1;
                ii(nzcount) = row;
                jj(nzcount) = numVar - 1;
                ss(nzcount) = -invNChol(1, 1) * Gt(1, 2);
                nzcount = nzcount + 1;
                ii(nzcount) = row + 1;
                jj(nzcount) = numVar - 2;
                ss(nzcount) = -invNChol(2, 2) * Gt(2, 1);
                nzcount = nzcount + 1;
                ii(nzcount) = row + 1;
                jj(nzcount) = numVar - 1;
                ss(nzcount) = -invNChol(2, 2) * Gt(2, 2);
                nzcount = nzcount + 1;
                ii(nzcount) = row + 1;
                jj(nzcount) = numVar;
                ss(nzcount) = -invNChol(2, 2) * Gt(2, 3);
                nzcount = nzcount + 1;
              end

              % update error term
              e(row) = invNChol(1, 1) * (r(k, ll) - temp2);
              if numCalib < 3
                e(row + 1) = invNChol(2, 2) *...
                  anglemod(b(k, ll) - (atan2(bb, aa) - x_est_temp(n, 3)));
              else
                e(row + 1) = invNChol(2, 2) *...
                  anglemod(b(k, ll) - (atan2(bb, aa) - x_est_temp(n, 3) -...
                  Theta_est_temp(3)));
              end
              row = row + 2;
            end
          end
          n = n + 1;
        end
        col = col + 3;
      end

      H = sparse(ii, jj, ss, (ns - numBatches) * 3 + numObs, numVar, nzmax);
      norms = colNorm(H); % could be included in the above loop for speedup
      G = spdiags(1 ./ norms, 0, cols(H), cols(H));

      % rank inference
      [C1, R1, P1] = spqr(H * G, -e, struct('permutation', 'matrix', ...
        'econ', cols(H)));
      sortR1 = sort(abs(diag(R1)), 'ascend');
      rankTol = 0;
      for rankIdx = 2:min(length(sortR1), 10)
        if sortR1(rankIdx) - sortR1(rankIdx - 1) > rankGap
          rankTol = sortR1(rankIdx);
          break;
        end
      end

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
      if rankTol == 0
        break;
      else
        update = G * spqr_solve(H * G, -e, struct('tol', rankTol));
      end
      x_est_temp = x_est_temp + [update(1:3:ns * 3) update(2:3:ns * 3) ...
        update(3:3:ns * 3)];
      x_est_temp(:, 3) = anglemod(x_est_temp(:, 3));
      l_est_temp = l_est_temp + [update(ns * 3 + 1:2:end - numCalib)...
        update(ns * 3 + 2:2:end - numCalib)];
      % deal with invisible landmarks that will become NaN
      for ll = 1:nl
        if isnan(l_est_temp(ll, 1)) || isnan(l_est_temp(ll, 2))
          l_est_temp(ll, :) = l_est(ll, :);
        end
      end
      Theta_est_temp = Theta_est_temp + update(end - numCalib + 1:end);
      if numCalib == 3
        Theta_est_temp(3) = anglemod(Theta_est_temp(3));
      end
    end

    s
    res

    % compute covariance
    Sigma = computeCov(H, e, numCalib);

    % compute covariance determinant
    sigmaDet = det(Sigma);

    % check if we need this batch
    if numBatches == 1 % first batch always taken
      x_est = x_est_temp;
      l_est = l_est_temp;
      Theta_est = Theta_est_temp;
      sigmaRecord = Sigma;
      sigmaDetRecord = sigmaDet;
    else
      % compute mutual information
      mi = 0.5 * log2(sigmaDetRecord / sigmaDet);

      % add batch if needed
      if mi > miTol
        x_est = x_est_temp;
        l_est = l_est_temp;
        Theta_est = Theta_est_temp;
        sigmaRecord = Sigma;
        sigmaDetRecord = sigmaDet;
      else
        % kick out the batch
        batchIdx(end - batchSize + 1:end) = [];
        numBatches = numBatches - 1;
      end
    end
    Theta_est
  end
end

Sigma = sigmaRecord;
