function [Lines Layout] = MCAVision (Lines, AffinityEnergyMatrix, MinGap, MagnetDiameter, Units, SizeX, SizeY, SaveFlag)
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


MinGap = MinGap/Units;
MagnetDiameter = MagnetDiameter/Units;

%Compute Layout using MDS
Layout = AffinityMDS(AffinityEnergyMatrix);

% Scale and discretize the layout to fit to a grid layout appropriate to
% layout the magnets
[X Y MaskImage] = CreateNanoMagnetMaskImage([Layout.X'], [Layout.Y'], MagnetDiameter, MinGap, SizeX, SizeY, SaveFlag);
Layout.X = X*Units; Layout.Y = Y*Units;
%subplot (2,2,1); plot(X, Y, 'o'); hold on; plot(X(1), Y(1), 'ro'); hold off;


% Save the final layout file. This can be used for fabrication.
if SaveFlag
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



%------------------------------------------------------------------
