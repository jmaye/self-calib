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

% This function implements a 2d motion model integrating the equations of
% motion using a first-order method. Furthermore, a noise variable can be added
% as input. Optionally, Jacobians with respect to state variable and noise
% variable can be outputted.

function [x_k Hx Hw] = motion(x_km1, u_k, T, w_k)

% if noise not specified, set it to zero
if nargin == 3
  w_k = zeros(2, 1);
end

% preliminary computations
theta_km1 = x_km1(3);
ct = cos(theta_km1);
st = sin(theta_km1);
C = [ct, 0;
     st, 0;
     0, 1];

% state update
x_k = x_km1 + T * C * (u_k + w_k);
x_k(3) = anglemod(x_k(3));

if nargout >= 2
  % Calculate the Jacobian with respect to state variable
  vk = u_k(1);
  Hx = eye(3);
  Hx(1, 3) = -T * vk * st;
  Hx(2, 3) =  T * vk * ct;
end

if nargout == 3
  % Calculate the Jacobian with respect to noise variable
  Hw = zeros(3, 2);
  Hw(1, 1) = T * ct;
  Hw(2, 1) = T * st;
  Hw(3, 2) = T;
end
