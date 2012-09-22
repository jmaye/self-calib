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

% This function implements an observation model for a laser range finder with
% calibration parameter observing a landmark. Furthermore, a noise variable
% can be added as input. Optionally, Jacobians with respect to state variable,
% landmark variable, calibration variable, and noise variable can be outputted.

function [y_k_l Gx Gl Gt Gn] = observation(x_k, x_l, Theta, n_k_l)

% if noise not specified, set it to zero
if nargin < 4
  n_k_l = zeros(2, 1);
end

% preliminary computations
theta_k = x_k(3);
ct = cos(theta_k);
st = sin(theta_k);
if length(Theta) < 3
  d = Theta;
  dst = d * st;
  dct = d * ct;
  a = x_l(1) - x_k(1) - dct;
  b = x_l(2) - x_k(2) - dst;
else
  dx = Theta(1);
  dy = Theta(2);
  psi = Theta(3);
  dxct = dx * ct;
  dxst = dx * st;
  dyct = dy * ct;
  dyst = dy * st;
  a = x_l(1) - x_k(1) - dxct + dyst;
  b = x_l(2) - x_k(2) - dxst - dyct;
end
temp1 = a^2 + b^2;
temp2 = sqrt(temp1);

if length(Theta) < 3
  y_k_l = [temp2; atan2(b, a) - theta_k] + n_k_l;
else
  y_k_l = [temp2; atan2(b, a) - theta_k - psi] + n_k_l;
end
y_k_l(2) = anglemod(y_k_l(2));

if nargout >= 2
  % Calculate the Jacobian with respect to state variable
  Gx = zeros(2, 3);
  Gx(1, 1) = -a / temp2;
  Gx(1, 2) = -b / temp2;
  if length(Theta) < 3
    Gx(1, 3) = (a * dst - b * dct) / temp2;
  else
    Gx(1, 3) = (a * (dxst + dyct) + b * (-dxct + dyst)) / temp2;
  end
  Gx(2, 1) = b / temp1;
  Gx(2, 2) = -a / temp1;
  if length(Theta) < 3
    Gx(2, 3) = -(a * dct + b * dst) / temp1 - 1;
  else
    Gx(2, 3) = (a * (-dxct + dyst) - b * (dxst + dyct)) / temp1 - 1;
  end
end

if nargout >= 3
  % Calculate the Jacobian with respect to landmark variable
  Gl = zeros(2, 2);
  Gl(1, 1) = a / temp2;
  Gl(1, 2) = b / temp2;
  Gl(2, 1) = -b / temp1;
  Gl(2, 2) = a / temp1;
end

if nargout >= 4
  % Calculate the Jacobian with respect to calibration variable
  if length(Theta) < 3
    Gt = zeros(2, 1);
    Gt(1, 1) = -(a * ct + b * st) / temp2;
    Gt(2, 1) = (-a * st + b * ct) / temp1;
  else
    Gt = zeros(2, 3);
    Gt(1, 1) = -(a * ct + b * st) / temp2;
    Gt(1, 2) = (a * st - b * ct) / temp2;
    Gt(2, 1) = (-a * st + b * ct) / temp1;
    Gt(2, 2) = -(a * ct + b * st) / temp1;
    Gt(2, 3) = -1;
  end
end

if nargout >= 5
  % Calculate the Jacobian with respect to noise variable
  Gn = eye(2, 2);
end
