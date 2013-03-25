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

% This script performs computations for Euler angles

syms xdot ydot zdot cx sx cy sy cz sz;
Rx = [1, 0, 0; 0, cx, -sx; 0, sx, cx];
Cx = [1, 0, 0; 0, cx, sx; 0, -sx, cx];
Rxt = [1, 0, 0; 0, cx, sx; 0, -sx, cx];
Rxdot = [0, 0, 0; 0, -sx, -cx; 0, cx, -sx] * xdot;
Ry = [cy, 0, sy; 0, 1, 0; -sy, 0, cy];
Cy = [cy, 0, -sy; 0, 1, 0; sy, 0, cy];
Ryt = [cy, 0, -sy; 0, 1, 0; sy, 0, cy];
Rydot = [-sy, 0, cy; 0, 0, 0; -cy, 0, -sy] * ydot;
Rz = [cz, -sz, 0; sz, cz, 0; 0, 0, 1];
Cz = [cz, sz, 0; -sz, cz, 0; 0, 0, 1];
Rzt = [cz, sz, 0; -sz, cz, 0; 0, 0, 1];
Rzdot = [-sz, -cz, 0; cz, -sz, 0; 0, 0, 0] * zdot;

% Z-X-Y convention
wskew = Rzdot * Rzt + Rz * Rxdot * Rxt * Rzt + ...
  Rz * Rx * Rydot * Ryt * Rxt * Rzt;

% Z-Y-X convention
%wskew = Rzdot * Rzt + Rz * Rydot * Ryt * Rzt + ...
%  Rz * Ry * Rxdot * Rxt * Ryt * Rzt;
