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
t_io = [0.195761055; -0.692404866; -1.158081532]; % [m]

% rotation from odometry to IMU center
C_io_roll = deg2rad(-10.676511765); % [rad]
C_io_pitch = deg2rad(-42.901313782); % [rad]
C_io_yaw = deg2rad(-56.549423218); % [rad]
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

errors = zeros(rows(data), 1);

% process data sequentially
for i = 1:rows(data)
  % translation from IMU body to world
  t_wi_b = [data(i, 19) data(i, 20) data(i, 21)]';

  % translational velocity of IMU body in world reference frame
  v_iw_b = [data(i, 8); data(i, 9); data(i, 10)];

  % rotational velocity of IMU body in world reference frame
  om_iw_b = [deg2rad(data(i, 11)); deg2rad(data(i, 12)); deg2rad(data(i, 13))];

  % translation from IMU center to world
  t_wi_c = [data(i, 2); data(i, 3); data(i, 4)];

  % rotation from IMU center to world
  C_wi_c_roll = deg2rad(data(i, 5));
  C_wi_c_pitch = deg2rad(data(i, 6));
  C_wi_c_yaw = deg2rad(data(i, 7));
  cx = cos(C_wi_c_roll);
  sx = sin(C_wi_c_roll);
  cy = cos(C_wi_c_pitch);
  sy = sin(C_wi_c_pitch);
  cz = cos(C_wi_c_yaw);
  sz = sin(C_wi_c_yaw);
  C_wi_c_x = [1, 0, 0; 0, cx, -sx; 0, sx, cx];
  C_wi_c_y = [cy, 0, sy; 0, 1, 0; -sy, 0, cy];
  C_wi_c_z = [cz, -sz, 0; sz, cz, 0; 0, 0, 1];
  C_wi_c = C_wi_c_z * C_wi_c_x * C_wi_c_y;

  % translational velocity of IMU center in world reference frame
  v_iw_c = v_iw_b + cross(om_iw_b, t_wi_c - t_wi_b);

  % rotational velocity of IMU center in world reference frame
  om_iw_c = om_iw_b;

  % translational velocity in IMU center reference frame
  v_ii_c = C_wi_c' * v_iw_c;

  % rotational velocity in IMU center reference frame
  om_ii_c = C_wi_c' * om_iw_c;

  % translational velocity in odometry reference frame
  v_oo = C_io' * (v_ii_c + cross(om_ii_c, t_io));

  % rotational velocity in odometry reference frame
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
  error(1) = anglemod(phi(i) - phi_pred(i));
  error(2) = v_rl(i) - v_rl_pred(i);
  error(3) = v_rr(i) - v_rr_pred(i);
  error(4) = v_fl(i) - v_fl_pred(i);
  error(5) = v_fr(i) - v_fr_pred(i);

  error = error + mvnrnd(zeros(1, cols(R)), R);

  % Mahalanobis distance
  errors(i) = error * inv(R) * error';

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
