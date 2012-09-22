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

% This file contains simulation and estimation parameters.

% translational speed measurement noise
v_var = 0.0044;

% rotational speed measurement noise
om_var = 0.0082;

% odometry covariance matrix
Q = diag([v_var; om_var]);

% range measurement noise
r_var = 9.0036e-04;

% bearing measurement noise
b_var = 6.7143e-04;

% lrf covariance matrix
R = diag([r_var; b_var]);

% true calibration parameters
Theta = [0.219; 0.1; pi / 4];

% calibration steps
steps = 2000;

% sampling time
T = 0.1;

% initial pose
x0 = [1; 1; pi / 4];

% true translational speed of the robot
v_true = zeros(steps, 1);
v_true(1:end) = 0.2;

% true rotational speed of the robot
om_true = zeros(steps, 1);
om_true(1:300) = 0;
om_true(301:end) = pi / 16;

% motion commands
u = [v_true, om_true];

% calibration parameters guess
Theta_hat = Theta + 0.1 .* randn(3, 1);

% maximum number of optimization iterations
maxIter = 100;

% tolerance for optimization
optTol = 1e-6;

% tolerance for TQR
rankTol = 10;
