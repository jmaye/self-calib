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

% This function computes the Jacobian and the error vector

function [J, e] = computeJacobian(x, l, theta, t, u, r, b, W, N, batchSize)

% angle normalization between -pi and pi
anglemod = @(x) atan2(sin(x), cos(x));

% noise-free case
if nargin < 9
  W = eye(3, 3);
  N = eye(2, 2);
end

% number of landmarks
nl = rows(l);

% number of state variables
ns = rows(x);
if nargin == 10
  if mod(ns, batchSize) == 1
    ns = ns - 1;
  end
end

% number of variables to estimate
numVar = ns * 3 + nl * 2 + length(theta);

% number of lrf observations
if nargin == 10
  numObs = 0;
  numBatches = 0;
  for i = 1:batchSize:rows(r) - 1
    if i + batchSize - 1 > rows(r)
      numObs = numObs + nnz(r(i + 1:end, :) > 0) * 2;
    else
      numObs = numObs + nnz(r(i + 1:i + batchSize - 1, :) > 0) * 2;
    end
    numBatches = numBatches + 1;
  end
else
  numObs = nnz(r(2:end, :) > 0) * 2;
  numBatches = 1;
end

% number of non-zero entries in the Jacobian
nzmax = (ns - numBatches) * 12 + numObs / 2 * 15;

% Jacobian initialization
ii = zeros(nzmax, 1);
jj = zeros(nzmax, 1);
ss = zeros(nzmax, 1);

% error term
e = zeros((ns - numBatches) * 3 + numObs, 1);

% row, column, and non-zeros counters
row = 1;
col = 1;
nzcount = 1;

% inverse squared root covariances
invWChol = sqrt(diag(1 ./ diag(W)));
invNChol = sqrt(diag(1 ./ diag(N)));

% intermediate matrices instantation
Hxk = zeros(3, 3);
Hxkm1 = zeros(3, 3);
Gxk = zeros(2, 3);
Gl = zeros(2, 2);
Gt = zeros(2, 3);

% build Jacobian and error vector
i = 2;
while i <= ns
  if nargin == 10 && mod(i - 1, batchSize) == 0
    i = i + 1;
    if i > ns
      break;
    end
    col = col + 3;
  end

  % some pre-computations
  stm1 = sin(x(i - 1, 3));
  ctm1 = cos(x(i - 1, 3));
  T = t(i) - t(i - 1);

  % Jacobian with respect to x_k
  Hxk(1, 1) = ctm1;
  Hxk(1, 2) = stm1;
  Hxk(2, 1) = -stm1;
  Hxk(2, 2) = ctm1;
  Hxk(3, 3) = 1;
  Hxk = -Hxk / T;

  % Jacobian with respect to x_k-1
  dx = x(i, :) - x(i - 1, :);
  Hxkm1(1, 1) = -ctm1;
  Hxkm1(1, 2) = -stm1;
  Hxkm1(1, 3) = -stm1 * dx(1) + ctm1 * dx(2);
  Hxkm1(2, 1) = stm1;
  Hxkm1(2, 2) = -ctm1;
  Hxkm1(2, 3) = -ctm1 * dx(1) - stm1 * dx(2);
  Hxkm1(3, 3) = -1;
  Hxkm1 = -Hxkm1 / T;

  % error term
  e(row) = invWChol(1, 1) * (u(i, 1) - 1.0 / T * (ctm1 * dx(1) + stm1 * dx(2)));
  e(row + 1) = invWChol(2, 2) * -1.0 / T * (-stm1 * dx(1) + ctm1 * dx(2));
  e(row + 2) = invWChol(3, 3) * (u(i, 2) - 1.0 / T * dx(3));

  % update sparse matrix filling
  ii(nzcount) = row;
  jj(nzcount) = col;
  ss(nzcount) = invWChol(1, 1) * Hxkm1(1, 1);
  nzcount = nzcount + 1;
  ii(nzcount) = row;
  jj(nzcount) = col + 1;
  ss(nzcount) = invWChol(1, 1) * Hxkm1(1, 2);
  nzcount = nzcount + 1;
  ii(nzcount) = row;
  jj(nzcount) = col + 2;
  ss(nzcount) = invWChol(1, 1) * Hxkm1(1, 3);
  nzcount = nzcount + 1;
  ii(nzcount) = row + 1;
  jj(nzcount) = col;
  ss(nzcount) = invWChol(2, 2) * Hxkm1(2, 1);
  nzcount = nzcount + 1;
  ii(nzcount) = row + 1;
  jj(nzcount) = col + 1;
  ss(nzcount) = invWChol(2, 2) * Hxkm1(2, 2);
  nzcount = nzcount + 1;
  ii(nzcount) = row + 1;
  jj(nzcount) = col + 2;
  ss(nzcount) = invWChol(2, 2) * Hxkm1(2, 3);
  nzcount = nzcount + 1;
  ii(nzcount) = row + 2;
  jj(nzcount) = col + 2;
  ss(nzcount) = invWChol(3, 3) * Hxkm1(3, 3);
  nzcount = nzcount + 1;
  ii(nzcount) = row;
  jj(nzcount) = col + 3;
  ss(nzcount) = invWChol(1, 1) * Hxk(1, 1);
  nzcount = nzcount + 1;
  ii(nzcount) = row;
  jj(nzcount) = col + 4;
  ss(nzcount) = invWChol(1, 1) * Hxk(1, 2);
  nzcount = nzcount + 1;
  ii(nzcount) = row + 1;
  jj(nzcount) = col + 3;
  ss(nzcount) = invWChol(2, 2) * Hxk(2, 1);
  nzcount = nzcount + 1;
  ii(nzcount) = row + 1;
  jj(nzcount) = col + 4;
  ss(nzcount) = invWChol(2, 2) * Hxk(2, 2);
  nzcount = nzcount + 1;
  ii(nzcount) = row + 2;
  jj(nzcount) = col + 5;
  ss(nzcount) = invWChol(3, 3) * Hxk(3, 3);
  nzcount = nzcount + 1;

  % update row/col counter
  row = row + 3;
  col = col + 3;

  % some pre-computations
  st1 = sin(x(i, 3));
  ct1 = cos(x(i, 3));

  % loop over the observations
  for j = 1:nl
    if r(i, j) > 0
      % some pre-computations
      dxct = theta(1) * ct1;
      dxst = theta(1) * st1;
      dyct = theta(2) * ct1;
      dyst = theta(2) * st1;
      aa = l(j, 1) - x(i, 1) - dxct + dyst;
      bb = l(j, 2) - x(i, 2) - dxst - dyct;
      temp1 = aa^2 + bb^2;
      temp2 = sqrt(temp1);

      % Jacobian of observation model with respect to state variable
      Gxk(1, 1) = -aa / temp2;
      Gxk(1, 2) = -bb / temp2;
      Gxk(1, 3) = (aa * (dxst + dyct) + bb * (-dxct + dyst)) / temp2;
      Gxk(2, 1) = bb / temp1;
      Gxk(2, 2) = -aa / temp1;
      Gxk(2, 3) = (aa * (-dxct + dyst) - bb * (dxst + dyct)) / temp1 - 1;
      Gxk = -Gxk;

      % update sparse matrix filling
      ii(nzcount) = row;
      jj(nzcount) = col;
      ss(nzcount) = invNChol(1, 1) * Gxk(1, 1);
      nzcount = nzcount + 1;
      ii(nzcount) = row;
      jj(nzcount) = col + 1;
      ss(nzcount) = invNChol(1, 1) * Gxk(1, 2);
      nzcount = nzcount + 1;
      ii(nzcount) = row;
      jj(nzcount) = col + 2;
      ss(nzcount) = invNChol(1, 1) * Gxk(1, 3);
      nzcount = nzcount + 1;
      ii(nzcount) = row + 1;
      jj(nzcount) = col;
      ss(nzcount) = invNChol(2, 2) * Gxk(2, 1);
      nzcount = nzcount + 1;
      ii(nzcount) = row + 1;
      jj(nzcount) = col + 1;
      ss(nzcount) = invNChol(2, 2) * Gxk(2, 2);
      nzcount = nzcount + 1;
      ii(nzcount) = row + 1;
      jj(nzcount) = col + 2;
      ss(nzcount) = invNChol(2, 2) * Gxk(2, 3);
      nzcount = nzcount + 1;

      % Jacobian of observation model with respect to landmark variable
      Gl(1, 1) = aa / temp2;
      Gl(1, 2) = bb / temp2;
      Gl(2, 1) = -bb / temp1;
      Gl(2, 2) = aa / temp1;
      Gl = -Gl;

      % update sparse matrix filling
      temp3 = ns * 3 + 1 + (j - 1) * 2;
      ii(nzcount) = row;
      jj(nzcount) = temp3;
      ss(nzcount) = invNChol(1, 1) * Gl(1, 1);
      nzcount = nzcount + 1;
      ii(nzcount) = row;
      jj(nzcount) = temp3 + 1;
      ss(nzcount) = invNChol(1, 1) * Gl(1, 2);
      nzcount = nzcount + 1;
      ii(nzcount) = row + 1;
      jj(nzcount) = temp3;
      ss(nzcount) = invNChol(2, 2) * Gl(2, 1);
      nzcount = nzcount + 1;
      ii(nzcount) = row + 1;
      jj(nzcount) = temp3 + 1;
      ss(nzcount) = invNChol(2, 2) * Gl(2, 2);
      nzcount = nzcount + 1;

      % Jacobian of observation model with respect to calibration variable
      Gt(1, 1) = -(aa * ct1 + bb * st1) / temp2;
      Gt(1, 2) = (aa * st1 - bb * ct1) / temp2;
      Gt(2, 1) = (-aa * st1 + bb * ct1) / temp1;
      Gt(2, 2) = -(aa * ct1 + bb * st1) / temp1;
      Gt(2, 3) = -1;
      Gt = -Gt;

      % update sparse matrix filling
      ii(nzcount) = row;
      jj(nzcount) = numVar - 2;
      ss(nzcount) = invNChol(1, 1) * Gt(1, 1);
      nzcount = nzcount + 1;
      ii(nzcount) = row;
      jj(nzcount) = numVar - 1;
      ss(nzcount) = invNChol(1, 1) * Gt(1, 2);
      nzcount = nzcount + 1;
      ii(nzcount) = row + 1;
      jj(nzcount) = numVar - 2;
      ss(nzcount) = invNChol(2, 2) * Gt(2, 1);
      nzcount = nzcount + 1;
      ii(nzcount) = row + 1;
      jj(nzcount) = numVar - 1;
      ss(nzcount) = invNChol(2, 2) * Gt(2, 2);
      nzcount = nzcount + 1;
      ii(nzcount) = row + 1;
      jj(nzcount) = numVar;
      ss(nzcount) = invNChol(2, 2) * Gt(2, 3);
      nzcount = nzcount + 1;

      % update error term
      e(row) = invNChol(1, 1) * (r(i, j) - temp2);
      e(row + 1) = invNChol(2, 2) * anglemod(b(i, j) - (atan2(bb, aa) - ...
        x(i, 3) - theta(3)));
      row = row + 2;
    end
  end
  i = i + 1;
end

% create sparse Jacobian matrix
if nzmax ~= nzcount - 1 || (ns - numBatches) * 3 + numObs ~= row - 1
  disp('warning: inconsistent data');
end
J = sparse(ii, jj, ss, rows(e), numVar, nzmax);
