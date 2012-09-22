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

% This function implements an inverse observation model for a laser range finder
% with calibration parameter observing a landmark. Optionally, Jacobians with
% respect to state variable, landmark measurement, and calibration variable can
% be outputted.

function [x_l Fx Fy Ft] = invObservation(x_k, y_k_l, Theta)

% preliminary computations
theta_k = x_k(3);
ct = cos(theta_k);
st = sin(theta_k);
r = y_k_l(1);
phi = y_k_l(2);
if length(Theta) < 3
  d = Theta;
  dst = d * st;
  dct = d * ct;
else
  dx = Theta(1);
  dy = Theta(2);
  psi = Theta(3);
  dxct = dx * ct;
  dxst = dx * st;
  dyct = dy * ct;
  dyst = dy * st;
end

if length(Theta) < 3
  x_l = x_k(1:2) + [dct; dst] + r * [cos(phi + theta_k); sin(phi + theta_k)];
else
  x_l = x_k(1:2) + [dxct - dyst; dxst + dyct] +...
    r * [cos(phi + psi + theta_k); sin(phi + psi + theta_k)];
  disp('using');
end

if nargout >= 2
  % Calculate the Jacobian with respect to state variable
  Fx = zeros(2, 3);
  Fx(1, 1) = 1;
  Fx(2, 2) = 1;
  if length(Theta) < 3
    Fx(1, 3) = -dst - r * sin(phi + theta_k);
    Fx(2, 3) = dct + r * cos(phi + theta_k);
  else
    Fx(1, 3) = -dxst - dyct - r * sin(phi + psi + theta_k);
    Fx(2, 3) = dxct - dyst + r * cos(phi + psi + theta_k);
  end
end

if nargout >= 3
  % Calculate the Jacobian with respect to observation variable
  Fy = zeros(2, 2);
  if length(Theta) < 3
    Fy(1, 1) = cos(phi + theta_k);
    Fy(1, 2) = -r * sin(phi + theta_k);
    Fy(2, 1) = sin(phi + theta_k);
    Fy(2, 2) = r * cos(phi + theta_k);
  else
    Fy(1, 1) = cos(phi + psi + theta_k);
    Fy(1, 2) = -r * sin(phi + psi + theta_k);
    Fy(2, 1) = sin(phi + psi + theta_k);
    Fy(2, 2) = r * cos(phi + psi + theta_k);
  end
end

if nargout >= 4
  % Calculate the Jacobian with respect to calibration variable
  if length(Theta) < 3
    Ft = zeros(2, 1);
    Ft(1, 1) = ct;
    Ft(2, 1) = st;
  else
    Ft = zeros(2, 3);
    Ft(1, 1) = ct;
    Ft(1, 2) = -st;
    Ft(1, 3) = r * cos(phi + psi + theta_k);
    Ft(2, 1) = st;
    Ft(2, 2) = ct;
    Ft(2, 3) = r * sin(phi + psi + theta_k);
  end
end
