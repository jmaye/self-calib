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

% rear-front axle distance
L = 1.81007; % [m]

% half-width rear axle
e_R = 0.6575; % [m]

% half-width front axle
e_F = 0.6425; % [m]

% rear right wheel radius
r_RR = 0.285004; % [m]

% rear left wheel radius
r_RL = 0.285004; % [m]

% front right wheel radius
r_FR = 0.280002; % [m]

% front left wheel radius
r_FL = 0.280002; % [m]

% translation from odometry to IMU
t_io = [0.261236; -0.0220428; -1.33797]; % [m]

% rotation from odometry to IMU
C_io_roll = deg2rad(7.26519); % [rad]
C_io_pitch = deg2rad(-21.1075); % [rad]
C_io_yaw = deg2rad(-45.4123); % [rad]
cx = cos(C_io_roll);
sx = sin(C_io_roll);
cy = cos(C_io_pitch);
sy = sin(C_io_pitch);
cz = cos(C_io_yaw);
sz = sin(C_io_yaw);
C_io = [cz * cy, -sz * cx + cz * sy * sx, sz * sx + cz * sy * cx;
        sz * cy, cz * cx + sz * sy * sx, -cz * sx + sz * sy * cx;
        -sy, cy * sx, cy * cx];

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

% process data sequentially
for i = 1:rows(data)
  % translation from IMU to world
  t_wi = [data(i, 2); data(i, 3); data(i, 4)];

  % rotation from IMU to world
  C_wi_roll = deg2rad(data(i, 5));
  C_wi_pitch = deg2rad(data(i, 6));
  C_wi_yaw = deg2rad(data(i, 7));
  cx = cos(C_wi_roll);
  sx = sin(C_wi_roll);
  cy = cos(C_wi_pitch);
  sy = sin(C_wi_pitch);
  cz = cos(C_wi_yaw);
  sz = sin(C_wi_yaw);
  C_wi = [cz * cy, -sz * cx + cz * sy * sx, sz * sx + cz * sy * cx;
          sz * cy, cz * cx + sz * sy * sx, -cz * sx + sz * sy * cx;
          -sy, cy * sx, cy * cx];

  % translational velocity of IMU in world reference frame
  v_iw = [data(i, 8); data(i, 9); data(i, 10)];

  % rotational velocity of IMU in world reference frame
  om_iw = [deg2rad(data(i, 11)); deg2rad(data(i, 12)); deg2rad(data(i, 13))];

  % translational velocity in IMU reference frame
  v_ii = C_wi' * v_iw;

  % rotational velocity in IMU reference frame
  om_ii = C_wi' * om_iw;

  % translational velocity in odometry reference frame
  v_oo = C_io' * (v_ii + cross(om_ii, t_io));

  % rotational velocity in odometry reference frame
  om_oo = C_io' * om_ii;

  % predicted translational speed measurable by odometry
  v_oo_x = v_oo(1);

  % predicted rotational speed measurable by odometry
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

  phi_pred(i) = atan2(L * om_oo_z, v_oo_x);

  % predicted steering angle of the front left wheel
  phi_L_pred =  atan2(L * om_oo_z, v_oo_x - e_F * om_oo_z);

  % predicted steering angle of the front right wheel
  phi_R_pred = atan2(L * om_oo_z, v_oo_x + e_F * om_oo_z);

  % predicted velocity of the rear left wheel
  v_rl_pred(i) = v_oo_x - e_R * om_oo_z;

  % predicted velocity of the rear right wheel
  v_rr_pred(i) = v_oo_x + e_R * om_oo_z;

  % predicted velocity of the front left wheel
  v_fl_pred(i) = (v_oo_x - e_F * om_oo_z) / cos(phi_L_pred);

  % predicted velocity of the front right wheel
  v_fr_pred(i) = (v_oo_x + e_F * om_oo_z) / cos(phi_R_pred);

end

norm(v_rl - v_rl_pred) / rows(data)
norm(v_rr - v_rr_pred) / rows(data)
norm(v_fl - v_fl_pred) / rows(data)
norm(v_fr - v_fr_pred) / rows(data)
anglemod = @(x) atan2(sin(x), cos(x));
norm(anglemod(phi - phi_pred)) / rows(data)

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
