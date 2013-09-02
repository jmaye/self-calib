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

% This function computes the marginal covariance, column space, null space out
% of a sparse Jacobian.

function [Omega, Sigma, SigmaP, CS, NS, miSum] = solveMarginal(J, i, epstol)

% default arguments
if nargin < 3
  epstol = 1e-4;
end

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
OmegaScaled = full((J_theta * G_theta)' * J_theta * G_theta - ...
  J_thetatQ * J_thetatQ');

% compute the thin SVD of Omega
[U, S, V] = svd(OmegaScaled, 0);

% compute the numerical rank
nrank = cols(S);
tol = rows(OmegaScaled) * S(1, 1) * epstol;
for i = cols(S):-1:1
  if S(i, i) > tol
    break;
  else
    nrank = nrank - 1;
  end
end

% compute the thin QR factorization of J_x
[Q, R, P, info] = spqr(J_x, struct('Q','Householder', 'econ', cols(J_x), ...
  'permutation', 'vector'));

% compute the Jacobian of the reduced system
J_thetatQ = spqr_qmult(Q, (J_theta)', 3);
J_thetatQ = J_thetatQ(:, 1:cols(J_x));
Omega = full(J_theta' * J_theta - J_thetatQ * J_thetatQ');

% compute the thin SVD of Omega
[U, S, V] = svd(Omega, 0);

% compute the numerical column space
CS = V(:, 1:nrank);

% compute the numerical null space
NS = V(:, end - (cols(V) - nrank) + 1:end);

% compute the projected covariance matrix
invS = zeros(cols(S), cols(S));
SigmaP = zeros(nrank, nrank);
miSum = 0;
for i = 1:nrank
  SigmaP(i, i) = 1.0 / S(i, i);
  miSum = miSum + log2(S(i, i));
  invS(i, i) = SigmaP(i, i);
end

Sigma = V * invS * U';
