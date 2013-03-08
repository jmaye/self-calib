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

% This script test Blake-Zisserman M-Estimator

% covariance matrix of the error
R = zeros(3, 3);
R(1, 1) = 1e-3;
R(2, 2) = 1e-3;
R(3, 3) = 1e-3;

% number of error terms to generate
numErrors = 1000;

% outliers ratio
outliersRatio = 0.1;

% errors
errors = zeros(numErrors, cols(R));

% cumulative errors
errorsCum = zeros(numErrors, 1);

% weights for each errors
errorWeights = zeros(numErrors, cols(R));

% generate errors following a normal distribution
for i = 1:numErrors
  if mod(i, 1 / outliersRatio) == 0
    errors(i, :) = mvnrnd(-100 + 200 .* rand(1, cols(R)), R); % outlier
  else
    errors(i, :) = mvnrnd(zeros(1, cols(R)), R); % inlier
  end
  errorsCum(i) = errors(i, :) * (diag(1 ./ diag(R))) * errors(i, :)';
  for j = 1:cols(R)
    errorWeights(i, j) = wbz(errors(i, j), R(j, j), 1e-3);
  end
end
