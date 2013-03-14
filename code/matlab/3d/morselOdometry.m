%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2013 by Jerome Maye                                            %
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

% this script processes a morsel log file

% angle normalization between -pi and pi
anglemod = @(x) atan2(sin(x), cos(x));

% translation from odometry to IMU center
t_io = [0.687484145; -0.255547822; -1.149330497]; % [m]

% rotation from odometry to IMU center
C_io_roll = deg2rad(15.619976997); % [rad]
C_io_pitch = deg2rad(-41.629974365); % [rad]
C_io_yaw = deg2rad(-46.374362946); % [rad]
cx = cos(C_io_roll);
sx = sin(C_io_roll);
cy = cos(C_io_pitch);
sy = sin(C_io_pitch);
cz = cos(C_io_yaw);
sz = sin(C_io_yaw);
C_io_x = [1, 0, 0; 0, cx, -sx; 0, sx, cx];
C_io_y = [cy, 0, sy; 0, 1, 0; -sy, 0, cy];
C_io_z = [cz, -sz, 0; sz, cz, 0; 0, 0, 1];
C_io = C_io_z * C_io_x * C_io_y;

% rear-front axle distance
L = 1.810069084; % [m]

% half-width rear axle
e_R = 1.315000057 / 2; % [m]

% half-width front axle
e_F = 1.284999967 / 2; % [m]

% rear right wheel radius
r_RR = 0.285004228; % [m]

% rear left wheel radius
r_RL = 0.285004228; % [m]

% front right wheel radius
r_FR = 0.280002475; % [m]

% front left wheel radius
r_FL = 0.280002475; % [m]

% odometry measurements and their predictions
v_fl = zeros(rows(data), 1);
v_fr = zeros(rows(data), 1);
v_rl = zeros(rows(data), 1);
v_rr = zeros(rows(data), 1);
phi = zeros(rows(data), 1);
v_fl_pred = zeros(rows(data), 1);
v_fr_pred = zeros(rows(data), 1);
v_rl_pred = zeros(rows(data), 1);
v_rr_pred = zeros(rows(data), 1);
phi_pred = zeros(rows(data), 1);

% Mahalanobis distances
errors = zeros(rows(data), 1);

% process data sequentially
for i = 1:rows(data)
  % translational velocity of IMU center in its own reference frame
  v_ii_c = [data(i, 8); data(i, 9); data(i, 10)];

  % rotational velocity of IMU center in its own reference frame
  om_ii_c = [deg2rad(data(i, 11)); deg2rad(data(i, 12)); deg2rad(data(i, 13))];

  % translational velocity of odometry center in its own reference frame
  v_oo = C_io' * (v_ii_c + cross(om_ii_c, t_io));

  % rotational velocity of odometry center in its own reference frame
  om_oo = C_io' * om_ii_c;

  % translational speed measurable by odometry
  v_oo_x = v_oo(1);

  % rotational speed measurable by odometry
  om_oo_z = om_oo(3);

  % velocity of the front left wheel
  v_fl(i) = r_FL * deg2rad(data(i, 14));

  % velocity of the front right wheel
  v_fr(i) = r_FR * deg2rad(data(i, 15));

  % velocity of the rear left wheel
  v_rl(i) = r_RL * deg2rad(data(i, 16));

  % velocity of the rear right wheel
  v_rr(i) = r_RR * deg2rad(data(i, 17));

  % steering wheel angle
  phi(i) = deg2rad(data(i, 18));

  % predicted steering angle of the virtual middle front wheel
  if v_oo_x == 0
    phi_pred(i) = phi(i);
    v_rl_pred(i) = v_rl(i);
    v_rr_pred(i) = v_rr(i);
    v_fl_pred(i) = v_fl(i);
    v_fr_pred(i) = v_fr(i);
    continue;
  end

  phi_pred(i) = atan(L * om_oo_z / v_oo_x);

  % predicted steering angle of the front left wheel
  phi_L_pred =  atan(L * om_oo_z / (v_oo_x - e_F * om_oo_z));

  % predicted steering angle of the front right wheel
  phi_R_pred = atan(L * om_oo_z / (v_oo_x + e_F * om_oo_z));

  % predicted velocity of the rear left wheel
  v_rl_pred(i) = v_oo_x - e_R * om_oo_z;

  % predicted velocity of the rear right wheel
  v_rr_pred(i) = v_oo_x + e_R * om_oo_z;

  % predicted velocity of the front left wheel
  v_fl_pred(i) = (v_oo_x - e_F * om_oo_z) / cos(phi_L_pred);

  % predicted velocity of the front right wheel
  v_fr_pred(i) = (v_oo_x + e_F * om_oo_z) / cos(phi_R_pred);

  % covariance matrix of the error
  R(1, 1) = 1e-3;
  R(2, 2) = 1e-3;
  R(3, 3) = 1e-3;
  R(4, 4) = 1e-3;
  R(5, 5) = 1e-3;

  % prediction error
  error = zeros(5, 1);
  error(1, 1) = anglemod(phi(i) - phi_pred(i));
  error(2, 1) = v_rl(i) - v_rl_pred(i);
  error(3, 1) = v_rr(i) - v_rr_pred(i);
  error(4, 1) = v_fl(i) - v_fl_pred(i);
  error(5, 1) = v_fr(i) - v_fr_pred(i);

  error = error + mvnrnd(zeros(cols(R), 1), R)';

  % Mahalanobis distance
  errors(i) = error' * inv(R) * error;

end

plot(v_rl, 'g');
hold on;
plot(v_rl_pred, 'r');
figure;
plot(v_rr, 'g');
hold on;
plot(v_rr_pred, 'r');
figure
plot(v_fl, 'g');
hold on;
plot(v_fl_pred, 'r');
figure;
plot(v_fr, 'g');
hold on;
plot(v_fr_pred, 'r');
figure;
plot(phi, 'g');
hold on;
plot(phi_pred, 'r');
figure;
[binCount, binPos] = hist(errors, 10000);
bar(binPos, binCount / trapz(binPos, binCount));
hold on;
x = 0:0.01:30;
plot(x, chi2pdf(x, cols(R)),'r');
