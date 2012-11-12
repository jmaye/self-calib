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

% Fix the way that the default dotted line and dash-dot lines look terrible
% when exported as eps

function fixDottedLine(filename)

fid = fopen(filename,'r');
tempfile = tempname;
outfid = fopen(tempfile,'w');

repeat = 1;
while repeat == 1
  thisLine = fgetl(fid);
  iStart = strfind(thisLine, '/DO { [.5');
  if iStart
    thisLine(iStart + 7:iStart + 8) = '03';
  end
  iStart = strfind(thisLine, '/DD { [.5');
  if iStart
    thisLine(iStart + 7:iStart + 9) = '1.5';
    thisLine(iStart + 10:end + 1) = [' ' thisLine(iStart + 10:end)];
  end

  if ~ischar(thisLine)
    repeat = 0;
  else
    fprintf(outfid,'%s\n',thisLine);
  end
end

fclose(fid);
fclose(outfid);
copyfile(tempfile, filename);
delete(tempfile);
