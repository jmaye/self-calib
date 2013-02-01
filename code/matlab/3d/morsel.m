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

% translation from odometry to vehicle
t_ov = [-0.0939582; -0.0220428; -1.36]; % [m]

% rotation from odometry to vehicle
R_ov_roll = deg2rad(-10.6197); % [rad]
R_ov_pitch = deg2rad(-3.47954); % [rad]
R_ov_yaw = deg2rad(-45.7612); % [rad]
R_ov_x = [1 0 0;
          0 cos(R_ov_roll) -sin(R_ov_roll);
          0 sin(R_ov_roll) cos(R_ov_roll)];
R_ov_y = [cos(R_ov_pitch) 0 sin(R_ov_pitch);
          0 1 0;
          -sin(R_ov_pitch) 0 cos(R_ov_pitch)];
R_ov_z = [cos(R_ov_yaw) -sin(R_ov_yaw) 0;
          sin(R_ov_yaw) cos(R_ov_yaw) 0;
          0 0 1];
R_ov = R_ov_z * R_ov_y * R_ov_x;

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
  % translation from vehicle to world
  t_vw = [data(i, 2); data(i, 3); data(i, 4)];

  % rotation from vehicle to world
  R_vw_roll = deg2rad(data(i, 5));
  R_vw_pitch = deg2rad(data(i, 6));
  R_vw_yaw = deg2rad(data(i, 7));
  R_vw_x = [1 0 0;
            0 cos(R_vw_roll) -sin(R_vw_roll);
            0 sin(R_vw_roll) cos(R_vw_roll)];
  R_vw_y = [cos(R_vw_pitch) 0 sin(R_vw_pitch);
            0 1 0;
            -sin(R_vw_pitch) 0 cos(R_vw_pitch)];
  R_vw_z = [cos(R_vw_yaw) -sin(R_vw_yaw) 0;
            sin(R_vw_yaw) cos(R_vw_yaw) 0;
            0 0 1];
  R_vw = R_vw_z * R_vw_y * R_vw_x;

  % translational velocity of vehicle in world reference frame
  v_vw = [data(i, 8); data(i, 9); data(i, 10)];

  % rotational velocity of vehicle in world reference frame
  om_vw = [deg2rad(data(i, 11)); deg2rad(data(i, 12)); deg2rad(data(i, 13))];

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

  % translational velocity in vehicle reference frame
  v_vv = R_vw' * v_vw;

  % rotational velocity in vehicle reference frame
  om_vv = R_vw' * om_vw;

  % translational velocity in odometry reference frame
  v_oo = R_ov' * (v_vv + cross(om_vv, t_ov));

  % rotational velocity in odometry reference frame
  om_oo = R_ov' * om_vv;

  % predicted translational speed measurable by odometry
  v_x = v_oo(1);

  % predicted rotational speed measurable by odometry
  om_z = om_oo(3);

  % predicted steering angle of the virtual middle front wheel
  if v_x == 0
    phi_pred(i) = phi(i);
    v_rl_pred(i) = v_rl(i);
    v_rr_pred(i) = v_rr(i);
    v_fl_pred(i) = v_fl(i);
    v_fr_pred(i) = v_fr(i);
    continue;
  end

  tan_phi_pred = L * om_z / v_x;
  phi_pred(i) = atan(tan_phi_pred);

  % predicted steering angle of the front left wheel
  phi_L_pred = atan(L * om_z / v_x / (1 - e_F * om_z / v_x));

  % predicted steering angle of the front right wheel
  phi_R_pred = atan(L * om_z / v_x / (1 + e_F * om_z / v_x));

  % predicted velocity of the rear left wheel
  v_rl_pred(i) = v_x - e_R * om_z;

  % predicted velocity of the rear right wheel
  v_rr_pred(i) = v_x + e_R * om_z;

  % predicted velocity of the front left wheel
  v_fl_pred(i) = (v_x - e_F * om_z) / cos(phi_L_pred);

  % predicted velocity of the front right wheel
  v_fr_pred(i) = (v_x + e_F * om_z) / cos(phi_R_pred);

end

norm(v_rl - v_rl_pred) / rows(data)
norm(v_rr - v_rr_pred) / rows(data)
norm(v_fl - v_fl_pred) / rows(data)
norm(v_fr - v_fr_pred) / rows(data)
