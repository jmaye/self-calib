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

clear all;

% simulation parameters
params;

% number of repetitions for each amplitude
rep = 50;

% sine wave frequency
frequency = 0.01;

% landmarks parameters
nl = 17;
minX = 0;
maxX = 15;
minY = 0;
maxY = 15;

% amplitude parameters
minAmplitude = 0;
maxAmplitude = 3;
amplitudeStep = 0.5;

% output for each repetition with a given amplitude
l_rep = zeros(nl, 2, rep);
Theta_tqr = zeros(rep, 3);
Sigma_tqr = zeros(3, 3, rep);
Theta_ls = zeros(rep, 3);
Sigma_ls = zeros(3, 3, rep);
Theta_ekf = zeros(rep, 3);
Sigma_ekf = zeros(3, 3, rep);
R1_rep = zeros(rep, 100);
R2_rep = zeros(rep, 100);

% amplitude index
ampIdx = 1;

% optimization parameters
optTol = 1e-6;
maxIter = 20;
rankTol = 0.02;

% initial guess for Theta
Theta_hat = [0.22; 0.12; 0.8];

warning off all;

% perform evaluation
for amplitude = minAmplitude:amplitudeStep:maxAmplitude
  amplitude

  if amplitude < 1
    rankTol = 0.028;
  end
  if amplitude >= 1 && amplitude < 3
    rankTol = 0.021;
  end
  if amplitude >= 3
    rankTol = 0.013;
  end

  for i = 1:rep
    i

    % setting up the scene
    u = genSinPath(steps, amplitude, frequency, T);
    l = genLandmarks(nl, [minX maxX], [minY maxY]);
    [x_true y_true th_true r b v om t] = simulate(u, x0, l, Q, R, Theta, T);

    % initial guess
    x_odom = odomInt([x_true(1), y_true(1), th_true(1)], [v, om], t);
    l_hat = initLandmarks(x_odom, Theta_hat, r, b);
    l_rep(:, :, i) = l;

    % TQR optimization
    disp('TQR');
    [x_est l_est Theta_est Sigma] = ls_slam_calib(x_odom, l_hat, ...
       Theta_hat, u, r, b, t, Q, R, maxIter, optTol, rankTol);
    Theta_tqr(i, :) = Theta_est';
    Sigma_tqr(:, :, i) = Sigma;
    Theta_est
    Sigma

    % Non-linear least squares without regularization
    disp('LS');
    [x_est l_est Theta_est Sigma R1 R2] = ls_slam_calib_prior(x_odom, l_hat, ...
      Theta_hat, u, r, b, t, Q, R, maxIter, optTol, x0, diag([1e-6;1e-6;1e-6]));
    Theta_ls(i, :) = Theta_est';
    Sigma_ls(:, :, i) = Sigma;
    R1_rep(i, :) = R1(end - 99:end);
    R2_rep(i, :) = R2(end - 99:end);
    Theta_est
    Sigma

    % EKF
    disp('EKF');
    [x_est, P] = ekf_slam_calib(x0, diag([1e-6; 1e-6; 1e-6]), Theta_hat, ...
      diag([1e-1; 1e-1; 1e-1]), u, r, b, t, Q, R);
    Theta_ekf(i, :) = x_est(end, end - 2:end);
    Sigma_ekf(:, :, i) = P(end - 2:end, end - 2:end, end);
    x_est(end, end - 2:end)
    P(end - 2:end, end - 2:end, end)
  end

  % save results for this amplitude
  ampRes(ampIdx).amplitude = amplitude;
  ampRes(ampIdx).rep = rep;
  ampRes(ampIdx).l_rep = l_rep;
  ampRes(ampIdx).R1_rep = R1_rep;
  ampRes(ampIdx).R2_rep = R2_rep;
  ampRes(ampIdx).Theta_tqr = Theta_tqr;
  ampRes(ampIdx).Sigma_tqr = Sigma_tqr;
  ampRes(ampIdx).Theta_ls = Theta_ls;
  ampRes(ampIdx).Sigma_ls = Sigma_ls;
  ampRes(ampIdx).Theta_ekf = Theta_ekf;
  ampRes(ampIdx).Sigma_ekf = Sigma_ekf;

  save('ampRes.mat', 'ampRes');

  % increment amplitude index
  ampIdx = ampIdx + 1;
end

warning on all;
