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

% This function exports a MATLAB figure to EPS

function exportEPSFig(figHandle, title, width, height)

% set the options
opts = struct('Format', 'eps', ...
              'Width', width, ...
              'Height', height, ...
              'Bounds', 'tight', ...
              'Color', 'rgb', ...
              'Renderer', 'painters', ...
              'Resolution', 1200, ...
              'FontMode', 'fixed', ...
              'FontSize', 4, ...
              'LineMode', 'fixed', ...
              'LineWidth', 0.5);

exportfig(figHandle, title, opts);
