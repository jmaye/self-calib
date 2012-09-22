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

% This function computes the covariance from the Jacobian of a least square
% estimate.

function Sigma = computeCov(H, e, n)

% pre-computations
colsH = cols(H);

% compute QR decomposition
[C, R, P] = spqr(H, -e, struct('permutation', 'matrix', 'econ', colsH));

% apply permutation to recover original order
R = P * R * P';

% init covariance
Sigma = zeros(n, n);

% computation loop, start from last element
Sigma_l = cols(Sigma);
for l = colsH:-1:colsH - n + 1
  temp1 = 0;
  Sigma_j = Sigma_l + 1;
  for j = l + 1:colsH
    if R(l, j) ~= 0
      temp1 = temp1 + R(l, j) * Sigma(Sigma_j, Sigma_l);
    end
    Sigma_j = Sigma_j + 1;
  end
  Sigma(Sigma_l, Sigma_l) = 1 / R(l, l) * (1 / R(l, l) - temp1);
  Sigma_i = Sigma_l - 1;
  for i = l - 1:-1:colsH - n + 1
    temp1 = 0;
    temp2 = 0;
    Sigma_j = Sigma_i + 1;
    for j = i + 1:l
      if R(i, j) ~= 0
        temp1 = temp1 + R(i, j) * Sigma(Sigma_j, Sigma_l);
      end
      Sigma_j = Sigma_j + 1;
    end
    Sigma_j = Sigma_l + 1;
    for j = l + 1:colsH
      if R(i, j) ~= 0
        temp2 = temp2 + R(i, j) * Sigma(Sigma_l, Sigma_j);
      end
      Sigma_j = Sigma_j + 1;
    end
    Sigma(Sigma_i, Sigma_l) = 1 / R(i, i) * (-temp1 - temp2);
    Sigma(Sigma_l, Sigma_i) = Sigma(Sigma_i, Sigma_l);
    Sigma_i = Sigma_i - 1;
  end
  Sigma_l = Sigma_l - 1;
end
