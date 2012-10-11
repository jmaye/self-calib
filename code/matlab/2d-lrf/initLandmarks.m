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

% This function inits the landmark positions for laser range data.

function l_hat = initLandmarks(x, Theta, r, b)

% number of landmarks
nl = cols(r);

% number of steps
steps = rows(r);

% init output
l_hat = zeros(nl, 2);

% landmarks found
l_found = zeros(nl);

% number of calibration parameters
numCalib = length(Theta);

for i = 1:steps
  for j = 1:nl
    if l_found(j) == 0 && r(i, j) > 0
      if numCalib < 3
        l_hat(j, 1) = x(i, 1) + Theta(1) * cos(x(i, 3)) + r(i, j) *...
          cos(b(i, j) + x(i, 3));
        l_hat(j, 2) = x(i, 2) + Theta(1) * sin(x(i, 3)) + r(i, j) *...
          sin(b(i, j) + x(i, 3));
      else
        l_hat(j, 1) = x(i, 1) + Theta(1) * cos(x(i, 3)) - Theta(2) *...
          sin(x(i, 3)) + r(i, j) * cos(b(i, j) + Theta(3) + x(i, 3));
        l_hat(j, 2) = x(i, 2) + Theta(1) * sin(x(i, 3)) + Theta(2) *...
          cos(x(i, 3)) + r(i, j) * sin(b(i, j) + Theta(3) + x(i, 3));
      end
      l_found(j) = 1;
    end
  end
  if nnz(l_found) == nl
    break;
  end
end
