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

% this script processes a morsel log file and compare with B-spline

figure;
plot(data(:,2),'g');
hold on;
plot(bsplinePoses(:,2),'r');
ylabel('x');
figure;
plot(data(:,3),'g');
hold on;
plot(bsplinePoses(:,3),'r');
ylabel('y');
figure;
plot(data(:,4),'g');
hold on;
plot(bsplinePoses(:,4),'r');
ylabel('z');
figure;
plot(deg2rad(data(:,7)),'g');
hold on;
plot(bsplinePoses(:,5),'r');
ylabel('yaw');
figure;
plot(deg2rad(data(:,6)),'g');
hold on;
plot(bsplinePoses(:,6),'r');
ylabel('pitch');
figure;
plot(deg2rad(data(:,5)),'g');
hold on;
plot(bsplinePoses(:,7),'r');
ylabel('roll');
figure;
plot(data(:,8),'g');
hold on;
plot(bsplinePoses(:,8),'r');
ylabel('v_x');
figure;
plot(data(:,9),'g');
hold on;
plot(bsplinePoses(:,9),'r');
ylabel('v_y');
figure;
plot(data(:,10),'g');
hold on;
plot(bsplinePoses(:,10),'r');
ylabel('v_z');
%figure;
%plot(log(:,11),'g');
%hold on;
%plot(bsplinePoses(:,11),'r');
%ylabel('\omega_x');
%figure;
%plot(log(:,12),'g');
%hold on;
%plot(bsplinePoses(:,12),'r');
%ylabel('\omega_y');
%figure;
%plot(log(:,13),'g');
%hold on;
%plot(bsplinePoses(:,13),'r');
%ylabel('\omega_z');
%figure;
%plot(log(:,14),'g');
%hold on;
