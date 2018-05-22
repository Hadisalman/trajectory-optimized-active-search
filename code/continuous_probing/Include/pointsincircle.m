function [points] = pointsincircle(xydata,radius,hk)
% POINTSINCIRCLE Identify points lying inside a circle

%   COMPLETELY VECTORIZED
%   Version 1.1 (Feb 29, 2012, 9:02 pm)
%   Version 1.0 (Feb 23, 2012) of this function was not uploaded

%   Find all points within a specified cut-off radius

%   Copyright, Sunil Anandatheertha Feb 23, 2012
%   http://soundcloud.com/sunilanandatheertha

%   INPUT:
%   xydata: a cell in the form xydata = {XvaluesOfCoordinates YvaluesOfCoordinates}.
%   radius: cut-off distance / cut-off radius
%   hk - h and k are the starting point of the radius vector.
%   OUTPUT:
%   points: this is a structure containing upto 2 values. 
%             points.in ---- cell of x- and y- coordinate data of all points inside the cut-off
%                                   radius
%             points.out ---- cell of x- and y- coordinate data of all points outside the cut-off
%                                   radius

dist = sqrt((hk(1) - xydata{1}(:)).^2   +   (hk(2) - xydata{2}(:)).^2); % distance calc.

% Find points inside circle
in = find(dist<radius);
points.in = {xydata{1}(in) xydata{2}(in)};

% Find points outside circle
% comment it out if its not needed
% its needed to run the testrun.m file though !!
out = find(dist>radius);  
points.out = {xydata{1}(out) xydata{2}(out)};

% ADDITIONAL --- USE TESTRUN.M TO SEE HOW TO USE THIS CODE

end