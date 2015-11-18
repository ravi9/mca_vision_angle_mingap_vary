function main (ImageFile, MinGap)

%% This is the main function that performs the following for images in the
%% Data directory. It assumes that edge detection has been performed and a
%% *.atts file containing all the edge segments is present for each image
%% 1. Reads the *.atts file
%% 2. Computes the affinity matrix.
%% 3. Computes the groups that would be found by traditional vision, i.e.
%%    optimization of the binary quadratic form based on stochastic search that
%%    is initialized by spectral methods.
%% 4. Computes the groups that would be created by magnet vision.
%% 5. Compares the groups in 3 and 4, using detection and false alarm rates.
%% 6. Saves the input edges, groups in 3 and 4 above as images and also all
%%    the local variables in a *.mat file.
%%
%% Function Call Graph
%% % main calls|- ReadEdgeAttsFile
%%             |- ComputeGroupingAffinities
%%             |- TraditionalVision
%%             |- MagnetVision calls |- NanodiskArray
%%             |- DrawSelectedLines
%%             |- RocAnalysis
%%
%% To run in matlab
%%   a = dir(fullfile('.','Data/*.atts')); % Creates a list of atts files
%%   for i=1:length(a), main(sprintf('Data/%s', a(i).name)); end;
%%   RocAnalysis ('Data/');
%%
%% To compile executable: mcc -m main
%% To execute the compiled version: ./run_main.sh /Applications/MATLAB_R2008a/ Data/rhino.pgm.atts
%% 
%% Last edited by Sudeep Sarkar (sarkar@cse.usf.edu)
%%
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


close all;

disp(['## Starting ' ImageFile]);

%% Reads the input edge segment file
Lines = ReadEdgeAttsFile (ImageFile);
%Lines = CreateExampleArrangement (); 

%% Computes the affinity matrix for the edges
AffinityEnergyMatrix = ComputeGroupingAffinities (Lines);


%% Computes the group of edges that optimizes the quadratic form using a
%% stochastic search stratgegy, which is initialized using the top
%% eigenvector of the afinitiy matrix.
tic
TOutput = TraditionalVision(AffinityEnergyMatrix);
t_vision = toc;
fprintf(1, '\n Time taken by traditional algorithm = %f (sec)', t_vision);
% figure; imagesc(TIOutput); title ( 'Software Vision'); colormap(gray);
% return;

%% Computes the saliencies (measured using a number between 0 and 1) of each edge based on magnet vision
%% Changes made to similar to that of fabrication

%MinGap = 10; %nm
MagnetThickness = 20; %nm
MagnetDiameter = 100; %nm, 
Units = 5; %nm, %pixel resolution in the bmp file is 5nm
SizeX = 1000; % maximum dimension of the 2D field in nm
SizeY = 1000;
SaveFlag = 0;% save a bmp file

%% 
% MinGap = 10; %nm
% MagnetThickness = 20; %nm
% MagnetDiameter = 100; %nm, 
% Units = 5; %nm, %pixel resolution in the bmp file is 5nm
% SizeX = 1000; % maximum dimension of the 2D field in nm
% SizeY = 1000;
% SaveFlag = 0;% save a bmp file

%%

% Finds the layout of magnets to use.
% Computes the states of the magnets under orthogonal global fields and uses
% it to figure out the coupled magnets.
[Lines Layout MOutput t_magnet] = MagnetVision (Lines, AffinityEnergyMatrix, MinGap, MagnetDiameter, MagnetThickness, Units, SizeX, SizeY, SaveFlag);
%subplot(2,2,3); imagesc(DrawSelectedLines (Lines, Layout,
%ones(size(Lines)))); colormap(gray);

%figure;
fprintf(1, '\n Time taken by nanomagnet computing = %f (sec)', t_magnet);
MCoupled = max (MOutput(1).CoupledStates, MOutput(2).CoupledStates)';


%% Performance comparison of the traditional vision group and the magnet
%% group using ROC. The traditional vision algorithm is considered the
%% "ground truth". Magnet vision assigns a saliency number between 0 and 1
%% to each edge segment. These saliencies are compared with a threshold to
%% decide salient/non-salient, which is then compared with the "ground
%% truth" At the end we have a curve to represent the performance on each
%% image. We save the edge map corresponding to the performance with 10%
%% false alarm

L = [];
for (i=1:length(Lines)), L = [L; Lines(i).length]; end;

Pd = []; Pf = []; ShowPerfThreshold = -1;
delt = (max(MCoupled) - min(MCoupled))/100;
for t=min(MCoupled):delt:max(MCoupled)
    p_d = sum(TOutput.*(MCoupled > t).*L)/sum(TOutput.*L);
    p_f = sum((1-TOutput).*(MCoupled > t).*L)/sum(L);
    Pd = [Pd p_d];
    Pf = [Pf p_f];
end;

% show visual performance at 10% false alarm
% Find the index for the false alarm rate of 0.1 or one that is closest
% (maximum value that is smaller to it).
i = find(Pf <= 0.1, 1);
if (i == length(Pf)), i = i -1; end; % avoid the 0 false alarm point, which is usually associated with 0 detection rate too.
ShowPerfThreshold = min(MCoupled) + delt*(i-1);
pdd = Pd(i); pff = Pf(i);
fprintf(1, '\n Detection Rate %f at %f false alarm', pdd, pff);
            
%% Save and plot stuff

%ImageFile = 'Data\zebra.pgm.atts'

ImageFileNameSplit = regexp(ImageFile , 'Data', 'split');
%resultDir = Results_mingap_6\zebra.pgm.atts
resultDir = ['Results_mingap_' num2str(MinGap) ImageFileNameSplit{2}];

% Save the input edges as an image
InputEdges = DrawSelectedLines (Lines, Layout, ones(size(Lines)));
imwrite(uint8(InputEdges), sprintf('%s_pd_%3.2f_pf_%3.2f_input.gif', resultDir, pdd, pff), 'GIF');

% Save the traditional vision edges as an image
TIOutput = DrawSelectedLines (Lines, Layout, TOutput);
imwrite(uint8(TIOutput), sprintf('%s_pd_%3.2f_pf_%3.2f_software.gif', resultDir, pdd, pff), 'GIF');

% Save the magnet vision edges as an image
MIOutput = DrawSelectedLines (Lines, Layout, MCoupled>ShowPerfThreshold );
imwrite(uint8(MIOutput), sprintf('%s_pd_%3.2f_pf_%3.2f_magnetic.gif', resultDir, pdd, pff), 'GIF');


subplot(1,4,2); imagesc(InputEdges); colormap(gray);
subplot(1,4,3); imagesc(TIOutput); title ('Software Vision'); colormap(gray);
subplot(1,4,4); imagesc(MIOutput); title ('Magnetic Vision'); colormap(gray);
% Draw the ROC on screen
subplot(1,4,1); plot (Pf, Pd, 'o-'); xlabel ('P_f'); ylabel('P_d'); hold on;
subplot(1,4,1); plot (Pf(ceil(ShowPerfThreshold *100)), Pd(ceil(ShowPerfThreshold *100)), 'ro'); 

 OFile = sprintf('%s_AllOutput.mat', resultDir);
 save(OFile, 'Lines', 'Layout', 'MOutput', 'MCoupled', 'TOutput', 'Pd', 'Pf', 't_vision', 't_magnet');

%figure; VisualizeConfiguration (Layout.X, Layout.Y, Lines)

