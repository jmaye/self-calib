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

% This function computes the marginal covariance of a sparse Jacobian.

function Sigma = margCovSparseJacobian(J, i, epstol)

% default arguments
if nargin < 3
  epstol = 1e-9;
end

% extract the part corresponding to the state/landmarks/...
J_x = J(:, 1:i - 1);

% extract the part corresponding to the calibration parameters
J_theta = J(:, i:end);

% compute the thin QR factorization of J_x
[Q, R, P, info] = spqr(J_x, struct('Q','Householder', 'econ', cols(J_x), ...
  'permutation', 'vector'));

% compute the intermediate Jacobian
J_thetatQ = spqr_qmult(Q, J_theta', 3);
J_thetatQ = J_thetatQ(:, 1:cols(J_x));
JThetaTrans = J_theta' * J_theta - J_thetatQ * J_thetatQ';

% compute the QR decomposition of JThetaTrans
[Q, R, P, info] = spqr(JThetaTrans, struct('Q','Householder', 'econ', ...
  cols(JThetaTrans), 'permutation', 'matrix'));

% compute the inverse of R
n = cols(R);
invR = zeros(n, n);
for k = n:-1:1
  if abs(R(k, k)) < epstol
    continue;
  end
  invR(k, k) = 1.0 / R(k, k);
  for i = k - 1:-1:1
    temp = 0;
    for j = i + 1:k
      temp = temp + R(i, j) * invR(j, k);
    end
    invR(i, k) = -1.0 / R(i, i) * temp;
  end
end

% compute the marginal covariance
Sigma = P * spqr_qmult(Q, invR, 2);
