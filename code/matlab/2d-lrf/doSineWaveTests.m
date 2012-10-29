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

% This script performs systematic evaluation on the sine wave path with
% varying amplitudes.

% number of repetitions for each amplitude
rep = 1;

% sine wave frequency
frequency = 0.01;

% landmarks parameters
nl = 17;
minX = 0;
maxX = 25;
minY = 0;
maxY = 25;

% amplitude parameters
minAmplitude = 0;
maxAmplitude = 2;
amplitudeStep = 1.0;

% output for each repetition with a given amplitude
Theta_hat_rep = zeros(rep, 3);
l_rep = zeros(nl, 2, rep);
Theta_tqr_mi = zeros(rep, 3);
Sigma_tqr_mi = zeros(3, 3, rep);
batchSize = zeros(rep, 1);
Theta_ls = zeros(rep, 3);
Sigma_ls = zeros(3, 3, rep);
Theta_ekf = zeros(rep, 3);
Sigma_ekf = zeros(3, 3, rep);

% amplitude index
ampIdx = 1;

% optimization parameters
optTol = 1e-9;
maxIter = 20;
batchSize = 100;
miThreshold = 0.5;
rankTol = 0.018;

% perform evaluation
for amplitude = minAmplitude:amplitudeStep:maxAmplitude
  amplitude
  for i = 1:rep
    i

    % setting up the scene
    params;
    u = genSinPath(steps, amplitude, frequency, T);
    l = genLandmarks(nl, [minX maxX], [minY maxY]);
    [x_true y_true th_true r b v om t] = simulate(u, x0, l, Q, R, Theta, T);

    % initial guess
    x_odom = odomInt([x_true(1), y_true(1), th_true(1)], [v, om], t);
    l_hat = initLandmarks(x_odom, Theta_hat, r, b);
    Theta_hat_rep(i, :) = Theta_hat';
    l_rep(:, :, i) = l;
    Theta_hat

    % TQR-MI optimization
    disp('TQR-MI');
    [x_est l_est Theta_est Sigma batchIdx] = ls_slam_calib_it(x_odom, l_hat, ...
      Theta_hat, u, r, b, t, Q, R, maxIter, optTol, batchSize, miThreshold, ...
      rankTol);
    Theta_tqr_mi(i, :) = Theta_est';
    Sigma_tqr_mi(:, :, i) = Sigma;
    Theta_est
    Sigma
    batchSize(i) = length(batchIdx);

    % Non-linear least squares without regularization
    disp('LS');
    [x_est l_est Theta_est Sigma] = ls_slam_calib_prior(x_odom, l_hat, ...
      Theta_hat, u, r, b, t, Q, R, maxIter, optTol, x0, diag([1e-6;1e-6;1e-6]));
    Theta_ls(i, :) = Theta_est';
    Sigma_ls(:, :, i) = Sigma;
    Theta_est
    Sigma

    % EKF
    disp('EKF');
    [x_est, P] = ekf_slam_calib(x0, diag([1e-6; 1e-6; 1e-6]), Theta_hat, ...
      diag([1e-3; 1e-3; 1e-3]), u, r, b, t, Q, R);
    Theta_ekf(i, :) = x_est(end, end - 2:end);
    Sigma_ekf(:, :, i) = P(end - 2:end, end - 2:end, end);
    x_est(end, end - 2:end)
    P(end - 2:end, end - 2:end, end)
  end

  % save results for this amplitude
  ampRes(ampIdx).amplitude = amplitude;
  ampRes(ampIdx).rep = rep;
  ampRes(ampIdx).Theta_hat_rep = Theta_hat_rep;
  ampRes(ampIdx).Theta_tqr_mi = Theta_tqr_mi;
  ampRes(ampIdx).Sigma_tqr_mi = Sigma_tqr_mi;
  ampRes(ampIdx).batchSize = batchSize;
  ampRes(ampIdx).Theta_ls = Theta_ls;
  ampRes(ampIdx).Sigma_ls = Sigma_ls;
  ampRes(ampIdx).Theta_ekf = Theta_ekf;
  ampRes(ampIdx).Sigma_ekf = Sigma_ekf;

  % increment amplitude index
  ampIdx = ampIdx + 1;
end
