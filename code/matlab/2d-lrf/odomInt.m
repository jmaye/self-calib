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

% This function integrates odometry from an initial point with motion commands
% measurements.

function x_odom = odomInt(x_0, u, t)

% number of steps
steps = rows(u);

% init odometry
x_odom = zeros(steps, 3);
x_odom(1, :) = x_0;

% angle normalization between -pi and pi
anglemod = @(x) atan2(sin(x), cos(x));

for i = 2:steps
  x_odom(i, 1) = x_odom(i - 1, 1) + (t(i) - t(i - 1)) *...
    cos(x_odom(i - 1, 3)) * u(i, 1);
  x_odom(i, 2) = x_odom(i - 1, 2) + (t(i) - t(i - 1)) *...
    sin(x_odom(i - 1, 3)) * u(i, 1);
  x_odom(i, 3) = anglemod(x_odom(i - 1, 3) + (t(i) - t(i - 1)) * u(i, 2));
end
