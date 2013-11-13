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

% This function implements a linear solver for incremental calibration
function x = linearSolver(J, e, i)

% extract the part corresponding to the state/landmarks/...
J_x = J(:, 1:i - 1);
G_x = normCol(J_x);

% extract the part corresponding to the calibration parameters
J_theta = J(:, i:end);
G_theta = normCol(J_theta);

% compute the thin QR factorization of J_x
[Q, R, P, info] = spqr(J_x * G_x, struct('Q','Householder', 'econ', ...
  cols(J_x), 'permutation', 'vector'));

% compute the Jacobian of the reduced system
J_thetatQ = spqr_qmult(Q, (J_theta * G_theta)', 3);
J_thetatQ = J_thetatQ(:, 1:cols(J_x));
OmegaTheta = full((J_theta * G_theta)' * J_theta * G_theta - ...
  J_thetatQ * J_thetatQ');

% compute the right-hand side of the reduced system
Qte = spqr_qmult(Q, e, 0);
Qte = Qte(1:cols(J_x), :);
b_theta = (J_theta * G_theta)' * e - J_thetatQ * Qte;

% compute the thin SVD of OmegaTheta
[U, S, V] = svd(OmegaTheta, 0);

% compute the numerical rank
nrank = cols(S);

% compute the least-squares solution for x_theta
x_theta = zeros(cols(J_theta), 1);
for i = 1:nrank
  x_theta = x_theta + U(:, i)' * b_theta * V(:, i) / S(i, i);
end

% compute the spqr solution for the rest of the state
x_x = spqr_solve((J_x * G_x)' * J_x * G_x, (J_x * G_x)' * e - ...
  (J_x * G_x)' * J_theta * G_theta * x_theta);
x_x = G_x * x_x;
x_theta = G_theta * x_theta;

x = [x_x; x_theta];
