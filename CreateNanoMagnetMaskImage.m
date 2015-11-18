function [X Y Field] = CreateNanoMagnetMaskImage(X, Y, diameter, gap, MaxX, MaxY, ImageFlag)

% Takes a continuous valued layout, streches and discretizes the layout so that we
% can place a magnet of "diameter" (nanometers) dimensions to be placed at
% each location and the total layout dimension is retricted to be withing
% MaxX and MaxY nanometers. The minimum distance between magnets is
% constrained to be "gap" nanometers. 
% The discretized locations are returned as X and Y
% An image is returned as "Field" for fabrication, depending if ImageFlag  == 1
% 
%% Copyright (c) 2011,  Sudeep Sarkar
%% All rights reserved
%%     This program is free software for research and educational purposes
%%     ONLY. For these purposes, you can redistribute it and/or modify it
%%     under the terms of the GNU General Public License as published by
%%     the Free Software Foundation, either version 3 of the License, or
%%     (at your option) any later version. This program is distributed in
%%     the hope that it will be useful, but WITHOUT ANY WARRANTY; without
%%     even the implied warranty of MERCHANTABILITY or FITNESS FOR A
%%     PARTICULAR PURPOSE.  See the GNU General Public License for more
%%     details. You should have received a copy of the GNU General Public
%%     License along with this program.  If not, see
%%     <http://www.gnu.org/licenses/>.

MagnetShape = CreateMagnetShape (diameter);
mX = round(MaxX/(diameter+gap));
mY = round(MaxY/(diameter+gap));
%% Figure out the scale parameter so that the layout is constrained by the
% overall size of the layout
GridUnit = round(diameter+gap);
scale = min(MaxX-diameter, MaxY-diameter)/max(max(X)-min(X), max(Y)-min(Y));
%fprintf (1, '\n Scaling design by %f', scale)
XX = scale*X/GridUnit; YY = scale*Y/GridUnit;

%subplot (2,2,1); plot(X, Y, 'o'); hold on; plot(X(1), Y(1), 'ro'); hold off
X = floor(XX-min(XX)+1); Y = floor(YY-min(YY)+1);
%subplot(2,2,2); plot(X, Y, 'o');

% Create am empty array
Field = zeros(mX, mY);
%% Mark the entry where the magnet should be centered
for (i=1:length(X))
    % Check if no other magnet has been located in the nearby vicinity
    % initiate a spiraling search for the next empty location
    found = 0; k = 1;
    while (found == 0)
        ii = [-k:k           k*ones(1, 2*k+1) -k:k               -k*ones(1, 2*k+1)];
        jj = [k*ones(1, 2*k+1) -k:k            -k*ones(1, 2*k+1) -k:k ];
        for kk=1:length(ii)
            xi = X(i) + ii(kk); yi = Y(i) + jj (kk);
            if ((xi > 0) && (yi > 0) && (xi < mX) && (yi < mY))
                if (Field (xi, yi) == 0)
                    Field(xi, yi) = 1; X(i) = xi; Y(i) = yi; found = 1;
                    break;
                end;
            end;
        end;
        k = k+1;
    end;
end;
%fprintf(1, '\n Wanted to embed %d magnets embedded %d magnet', length(X), sum(sum(Field)));
X = X *GridUnit; Y = Y * GridUnit;
%subplot(2,2,3); plot(X, Y, 'o');

if (ImageFlag == 0), return; end;


ExpandedField = zeros(MaxX, MaxY);
for (i=1:size(Field, 1))
    for (j=1:size(Field, 2))
        if (Field(i, j) == 1)
            ExpandedField (round(i*GridUnit), round(j*GridUnit)) = 1;
        end;
    end;
end;
% Convolve to place the circular magnet at each location to generate the
% bitmap
Field = conv2(ExpandedField, MagnetShape, 'same');



%% Create mask for the nanomagnet shape (circle)
function y = CreateMagnetShape (diameter)
xc = ceil(diameter/2); yc = ceil(diameter/2); r = floor(diameter/2);
y = zeros(ceil(diameter), ceil(diameter));
for i=xc-r:xc+r
    for j=yc-r:yc+r
        if ((i>0)&&(j>0)&&((i-xc)^2+(j-yc)^2 < r^2))
            y(ceil(i), ceil(j)) = 1;
        end;
    end;
end;


       
