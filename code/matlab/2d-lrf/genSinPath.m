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

% This function generates a sinus path.

function u = genSinPath(steps, amplitude, frequency, T)

% true translational speed of the robot
v = zeros(steps, 1);

% true rotational speed of the robot
om = zeros(steps, 1);

% generate sin wave
t = 0:T:steps * T;
sint = zeros(length(t), 1);
cost = zeros(length(t), 1);
for i = 1:length(t)
  sint(i) = amplitude * sin(2 * pi * frequency * t(i));
  cost(i) = amplitude * cos(2 * pi * frequency * t(i));
end

% sample path and generate commands
theta = 0;
for i = 2:length(t)
  theta_new = atan(cost(i));
  v(i) = 2 * pi * frequency *...
    (sqrt(1 + cost(i)^2) + sqrt(1 + cost(i - 1)^2)) / 2;
  om(i) = (theta_new - theta) / T;
  theta = theta_new;
end

% motion commands
u = [v, om];
