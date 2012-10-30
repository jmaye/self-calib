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
% varying amplitudes and varying parameters.

clear all;

% simulation parameters
params;

% number of repetitions for each amplitude and parameter
rep = 10;

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
amplitudeStep = 1.0;

% rank tolerance parameters
minRankTol = 0;
maxRankTol = 0.02;
rankTolStep = 0.001;

% output for each repetition with a given amplitude
Theta_tqr = zeros(rep, 3);
Sigma_tqr = zeros(3, 3, rep);

%% amplitude index
ampIdx = 1;

% optimization parameters
optTol = 1e-6;
maxIter = 15;

% initial guess for Theta
Theta_hat = [0.22; 0.12; 0.8];

warning off all;

% perform evaluation
for amplitude = minAmplitude:amplitudeStep:maxAmplitude
  amplitude

  for rankTol = minRankTol:rankTolStep:maxRankTol
    rankTol

    for i = 1:rep
      i

      % setting up the scene
      u = genSinPath(steps, amplitude, frequency, T);
      l = genLandmarks(nl, [minX maxX], [minY maxY]);
      [x_true y_true th_true r b v om t] = simulate(u, x0, l, Q, R, Theta, T);

      % initial guess
      x_odom = odomInt([x_true(1), y_true(1), th_true(1)], [v, om], t);
      l_hat = initLandmarks(x_odom, Theta_hat, r, b);

      % TQR optimization
      disp('TQR');
      [x_est l_est Theta_est Sigma] = ls_slam_calib(x_odom, l_hat, ...
         Theta_hat, [v, om], r, b, t, Q, R, maxIter, optTol, rankTol);
      Theta_tqr(i, :) = Theta_est';
      Sigma_tqr(:, :, i) = Sigma;
      Theta_est
      Sigma
    end

    % save results for this amplitude and parameter
    paramRes(ampIdx).amplitude = amplitude;
    paramRes(ampIdx).rankTol = rankTol;
    paramRes(ampIdx).rep = rep;
    paramRes(ampIdx).Theta_tqr = Theta_tqr;
    paramRes(ampIdx).Sigma_tqr = Sigma_tqr;

    % save file to disk
    save('paramRes.mat', 'paramRes');

    % increment amplitude index
    ampIdx = ampIdx + 1;
  end
end

warning on all;
