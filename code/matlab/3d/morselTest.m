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

% translation from odometry to IMU
t_io = [0.261235952; -0.022042794; -1.337972522]; % [m]

% rotation from odometry to IMU
C_io_roll = deg2rad(7.265191078); % [rad]
C_io_pitch = deg2rad(-21.107488632); % [rad]
C_io_yaw = deg2rad(-45.412288666); % [rad]
cx = cos(C_io_roll);
sx = sin(C_io_roll);
cy = cos(C_io_pitch);
sy = sin(C_io_pitch);
cz = cos(C_io_yaw);
sz = sin(C_io_yaw);
C_io = [cz * cy, -sz * cx + cz * sy * sx, sz * sx + cz * sy * cx;
        sz * cy, cz * cx + sz * sy * sx, -cz * sx + sz * sy * cx;
        -sy, cy * sx, cy * cx];

% process data sequentially
for i = 450:rows(data)
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

  % translation from odometry to world
  t_wo = [data(i, 19); data(i, 20); data(i, 21)];

  % rotation from odometry to world
  C_wo_roll = deg2rad(data(i, 22));
  C_wo_pitch = deg2rad(data(i, 23));
  C_wo_yaw = deg2rad(data(i, 24));
  cx = cos(C_wo_roll);
  sx = sin(C_wo_roll);
  cy = cos(C_wo_pitch);
  sy = sin(C_wo_pitch);
  cz = cos(C_wo_yaw);
  sz = sin(C_wo_yaw);
  C_wo = [cz * cy, -sz * cx + cz * sy * sx, sz * sx + cz * sy * cx;
          sz * cy, cz * cx + sz * sy * sx, -cz * sx + sz * sy * cx;
          -sy, cy * sx, cy * cx];

  % translational velocity of odometry in world reference frame
  v_ow = [data(i, 25); data(i, 26); data(i, 27)];

  % rotational velocity of odometry in world reference frame
  om_ow = [deg2rad(data(i, 28)); deg2rad(data(i, 29)); deg2rad(data(i, 30))];

  % translational velocity in odometry reference frame
  v_oo_2 = C_wo' * v_ow;

  % rotational velocity in odometry reference frame
  om_oo_2 = C_wo' * om_ow;

  t_wi_2 = [data(i,31) data(i,32) data(i,33)];
  C_wi_2 = [data(i,34) data(i,35) data(i,36);data(i,37) data(i,38) data(i,39); data(i,40) data(i,41) data(i,42)];
  t_wo_2 = [data(i,43) data(i,44) data(i,45)];
  C_wo_2 = [data(i,46) data(i,47) data(i,48);data(i,49) data(i,50) data(i,51); data(i,52) data(i,53) data(i,54)];

end
