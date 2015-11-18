function [Lines Layout Images t_mag] = MagnetVision (Lines, AffinityEnergyMatrix, MinGap, MagnetDiameter, MagnetThickness, Units, SizeX, SizeY, SaveFlag)

%% This function produces the Magnet mask for OOMMF simulations. It produces bmp file 
%% called VisionDots.bmp can be included in the *.mif file. The output
%% Layout is also returned as structure of X and y coordinates, Layout.X
%% and Layout.Y, respectively.
%%
%% MagnetDiameter: -- Diagmeter of the magnet (nm)
%%
%% MinGap: Minimum gap between the magnets (nm)
%%
%% Units: Length of each pixel in nm. By choosing
%%        large values of this parameter (i.e. 2, 5, 10, etc)
%%        one can reduce the size of the output image (bmp)
%%
%% SizeX and SizeY: size of the layout bitmap in nm 
%%
%% SaveFlag == 1 to save bitmap output
%% 
%% Output: 
%% Files VisionDots*.bmp is created, one corresponding to each component. Each pixel is 1nm
%% Files Coords*.txt is created, listing the (x y) locations of the center of the magnets
%% for each component
%% 
% 
%% Example: MCAVision (Lines, AffinityEnergyMatrix, 0.2, 30, 100, 10)
%%  Creates a design with 100 nm diameter magnets, with minimum spacing of
%%  30nm and each pixel in the bmp file representing 10nm. It thresholds the
%%  affinity matrix by 0.2.
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

tic;
MinGap = MinGap/Units;
MagnetDiameter = MagnetDiameter/Units;

%%%---------------------------------------------------
%%Compute Layout using MDS


N = size(AffinityEnergyMatrix, 1);


Distances = ChangeToDistances(AffinityEnergyMatrix);
% Distances = RestoreMetricProperty (Distances);
% Distances = AffinityToMetricDistance (AffinityEnergyMatrix);
Y = cmdscale(Distances);
%X = [0; Y(:,1)]; Y = [0; Y(:, 2)]; % Input cell at (0,0) useful for
%ferromagnetically coupled cells. Not needed for anti-ferromagnetic
%coupling.
X = Y(:,1); Y = Y(:, 2); 

%%%---------------------------------------------------
%% Scale and discretize the layout to fit to a grid layout appropriate to
%% layout the magnets
%% Takes a continuous valued layout, streches and discretizes the layout so that we
%% can place a magnet of "diameter" (nanometers) dimensions to be placed at
%% each location and the total layout dimension is retricted to be withing
%% MaxX and MaxY nanometers. The minimum distance between magnets is
%% constrained to be "gap" nanometers. 
%% The discretized locations are returned as X and Y

mX = round(SizeX/(MagnetDiameter + MinGap));
mY = round(SizeY/(MagnetDiameter + MinGap));

% Figure out the scale parameter so that the layout is constrained by the
% overall size of the layout
GridUnit = round(MagnetDiameter + MinGap);
scale = min(SizeX - MagnetDiameter, SizeY - MagnetDiameter)/max(max(X)-min(X), max(Y)-min(Y));
%fprintf (1, '\n Scaling design by %f', scale)
XX = scale*X/GridUnit; YY = scale*Y/GridUnit;
%subplot (2,2,1); plot(X, Y, 'o'); hold on; plot(X(1), Y(1), 'ro'); hold off
X = floor(XX-min(XX)+1); Y = floor(YY-min(YY)+1);
%subplot(2,2,2); plot(X, Y, 'o');

% Create am empty array
Field = zeros(mX, mY);
% Mark the entry where the magnet should be centered
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
fprintf(1, '\n Wanted to embed %d magnets embedded %d magnet', length(X), sum(sum(Field)));
X = X *GridUnit; Y = Y * GridUnit;
Layout.X = X(1:N)*Units; Layout.Y = Y(1:N)*Units; Layout.indices = [1:N];

%pause(1);


%% Save the final layout file. This can be used for fabrication.
if SaveFlag
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
    MagnetShape = CreateMagnetShape (diameter);
    Field = conv2(ExpandedField, MagnetShape, 'same');
    
    fprintf(1, '\n Writing coordinates of magnets in text file Coords.txt');
    fp = fopen (sprintf('Coords_%s_%d_g_%d_d_%d_u_%d.txt', ImageFile, i, MinGap*Units, MagnetDiameter*Units, Units), 'w');
    fprintf(fp, 'Material = 80Co20Pt');
    fprintf(fp, '\nRadius = %f um', MagnetDiameter*Units/1000);
    fprintf(fp, '\nThickness = %f um', 40/1000);
    fprintf(fp, '\nmagnet       X (um)    Y (um)');
    for (ii=1:length(X)) fprintf(fp, '\nm%d = (%f, %f)', ii, X(ii)/1000, Y(ii)/1000); end;
    fclose (fp);
    fprintf(1, '\n Image mask in VisionDots.bmp');
    imwrite(mat2gray(MaskImage), sprintf('VisionDots_%s_%f_%d_g_%d_d_%d_u_%d.bmp', ImageFile, AffinityThreshold, i, MinGap*Units, MagnetDiameter*Units, Units), 'bmp');
    figure; imagesc(MaskImage'); colormap('gray');
end;

%% % This part simulates the protocol used for vision computing using an
% nano-magnet array. We have experimented with combination of different
% initializations and different global fields. The protocol that we have
% setlled upon consists of applying external field along two orthogonal
% directions and computing the sine of the angles of each magnet with
% respect to the two external field directions. Coupled magnets should not
% deviate when external field is applied.
t_mag = toc;

%Input = [pi/2; -10*ones(N, 1)]; 
Input = -10*ones(N, 1); % No Input Node
h = 0.5; % external field strength

fprintf(1, '\n Relaxing nanodisk array with pi/4 orientation as initial states for all disks and external field along x.');
[Angles States TimeTaken] = NanodiskArray (X, Y, MagnetDiameter/2, MagnetThickness, Input, pi*ones(size(X))/4, [h 0]);
%[Angles States TimeTaken] = NanodiskArray (X, Y, MagnetRadius, MagnetThickness, Input, zeros(size(X)), [0 0]);
fprintf(1, '\n Time Taken = %f picoseconds', TimeTaken/1e-12);
t_mag = t_mag + TimeTaken;
Images(1).CoupledStates = abs(sin(Angles/2))'; %0 if aligned with external field, 1 if anti-aligned %abs(States);
Images(1).OutputImage = DrawSelectedLines (Lines, Layout, Images(1).CoupledStates);



fprintf(1, '\n Relaxing nanodisk array with pi/4 orientation as initial states for all disks and external field along y.');
[Angles States TTimeTaken]  = NanodiskArray (X, Y, MagnetDiameter/2, MagnetThickness, Input, pi*ones(size(X))/4, [0 h]);
%[Angles States TTimeTaken]  = NanodiskArray (X, Y, MagnetRadius, MagnetThickness, Input, pi*ones(size(X))/2, [0 0]);
fprintf(1, '\n Time Taken = %f picoseconds', TimeTaken/1e-12);
t_mag = t_mag + TimeTaken;
Images(2).CoupledStates =  abs(cos(Angles/2))'; %abs(States);
Images(2).OutputImage = DrawSelectedLines (Lines, Layout, Images(2).CoupledStates);



%-----------------------------------------------------------------------



%% Restore metric property, i.e. d(x, z) <= d(X, y) + d(y, z)
function y = RestoreMetricProperty (D)
N = size(D, 1);

c = 0;
for i=1:N 
    for j=1:N
        for k=1:N
            c = max(c, max(0, D(i, j) - (D(i, k) + D(k, j))));
        end;
    end;
end;
fprintf(1, '\n Restoring metric property by adding %f', c);
y = D+c;
for i=1:N, y(i,i) = 0; end;

return;
%-----------------------------------------------------------------------

%% Restore metric property, i.e. d(x, z) <= d(X, y) + d(y, z)
function D = AffinityToMetricDistance (A)
N = size(A, 1);

cc = 0;
for i=1:N 
    for j=i+1:N
        for k=1:N
            c = cc;
            while ((A(i, j) + c)^(-1/3) > ((A(i, k) + c)^(-1/3) + (A(k, j) + c)^(-1/3)))
                c = c +1;
            end;
            cc = max (cc, c);
        end;
    end;
end;
fprintf(1, '\n Restoring metric property by adding %f', cc);
D = (A+cc).^(-1/3);

fprintf(1, '\n min and max values of the distance matrix [%f %f]', min(min(D)), max(max(D)));

for (i=1:N) D(i, i) = 0.0; end; %% Restore diagonal to zero distances
D = 0.5*(D+D');
return;

%% Change Affinities to Distances. Affinities could be larger than one.

function D = ChangeToDistances (A)


%D = (A+eps).^-(1.0/3);
%D = (-log(A)); %maps affinity in the range 0 to 1 into the range infinity and zero -- is this
% linear manifold approximation? measure of "repulsion" if 0 < A < 1.

%sum(sum(D<0))

%D = (D > 0).*D; %% to remove very small negative created values due to numerical issues

D = (A+10).^(-1/3);

%D = (log(A+1)).^(1/3);

fprintf(1, '\n min and max values of the distance matrix [%f %f]', min(min(D)), max(max(D)));
N = size(D, 1);
for (i=1:N) D(i, i) = 0.0; end; %% Restore diagonal to zero distances
D = 0.5*(D+D'); % make the matrix symmetric.
%-----------------------------------------------------------------------



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
%-----------------------------------------------------------------------
