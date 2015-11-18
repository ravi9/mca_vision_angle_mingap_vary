function [Images t_mag] = VisionCompute (Lines, Layout, MagnetThickness, MagnetRadius)
%
% This function simulates the protocol used for vision computing using an
% nano-magnet array. We have experimented with combination of different
% initializations and different global fields. The protocol that we have
% setlled upon consists of applying external field along two orthogonal
% directions and computing the sine of the angles of each magnet with
% respect to the two external field directions. Coupled magnets should not
% deviate when external field is applied.
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

tic
X = Layout.X; Y = Layout.Y;
if (size(X, 2) ~= 1)       
    X = X'; Y = Y';
end;
N = length(X);
Input = -10*ones(N, 1); % No Input Node

% count = 1;
% fprintf(1, '\n Compute steady state. Initialize with random directions and let it settle.');
% [Angles States] = NanodiskArray (X, Y, MagnetRadius, MagnetThickness, Input, pi*ones(size(X))/2, [0 0]);
% 
% %Coupled = sin(Angles);
% Coupled = abs(States);
% 
% Images(count).CoupledStates = Coupled;
% Images(count).OutputImage = DrawSelectedLines (Lines, Layout, Coupled);
% count = count + 1;
% Images(count).CoupledStates = Coupled;
% Images(count).OutputImage = DrawSelectedLines (Lines, Layout, Coupled);
% count = count + 1;
%     
% InitialState = NanodiskArray (X, Y, MagnetRadius, MagnetThickness, Input, pi*ones(size(X))/2, [0 0]);
% Coupled = cos(InitialState);
% Images(count).CoupledStates = Coupled;
% Images(count).OutputImage = DrawSelectedLines (Lines, Layout, Coupled);
% count = count + 1;

    
count = 1;
for h = 0.5:0.5 %0.01:0.01:0.05
    %     NextState = NanodiskArray (X, Y, MagnetRadius, MagnetThickness, Input, InitialState, [0 h]);
    %     Coupled = ([cos(NextState) sin(NextState)].*[cos(InitialState) sin(InitialState)]);
    %     Coupled = (1 + sum(Coupled'))/2;
    t_mag = toc;
    fprintf(1, '\n Relaxing nanodisk array with pi/4 orientation as initial states for all disks and external field along x.');
    %[Angles States TimeTaken] = NanodiskArray (X, Y, MagnetRadius, MagnetThickness, Input, pi*ones(size(X))/4, [h 0]);
    [Angles States TimeTaken] = NanodiskArray (X, Y, MagnetRadius, MagnetThickness, Input, zeros(size(X)), [0 0]);
    fprintf(1, '\n Time Taken = %f picoseconds', TimeTaken/1e-12);
    t_mag = t_mag + TimeTaken;
    Images(count).CoupledStates = abs(sin(Angles))'; abs(States);
    Images(count).OutputImage = DrawSelectedLines (Lines, Layout, Images(count).CoupledStates);
    count = count + 1;

    pause (1);
    %     NextState = NanodiskArray (X, Y, MagnetRadius, MagnetThickness, Input, InitialState, [h 0]);
    %     Coupled = ([cos(NextState) sin(NextState)].*[cos(InitialState) sin(InitialState)]);
    %     Coupled = (1 + sum(Coupled'))/2;
    fprintf(1, '\n Relaxing nanodisk array with pi/4 orientation as initial states for all disks and external field along y.');
    %[Angles States TTimeTaken]  = NanodiskArray (X, Y, MagnetRadius, MagnetThickness, Input, pi*ones(size(X))/4, [0 h]);
    [Angles States TTimeTaken]  = NanodiskArray (X, Y, MagnetRadius, MagnetThickness, Input, pi*ones(size(X))/2, [0 0]);
    fprintf(1, '\n Time Taken = %f picoseconds', TimeTaken/1e-12);
    t_mag = t_mag + TimeTaken;
    Images(count).CoupledStates =  abs(cos(Angles))'; %abs(States);
    Images(count).OutputImage = DrawSelectedLines (Lines, Layout, Images(count).CoupledStates);
    count = count + 1;
end
%save VisionCompute.mat InitialState Images Lines Layout

% % Slowly change the state of the cell at (0, 0) -- the first one in the
% % list by 180 degree from the starting state. 
% This strategy just picks up the cluster at (0,0). Method deprecated
% 
% StartState = InitialState; 
% for t=1:10
%     Input(1) = InitialState(1) + t*pi/10; 
%     NextState = NanodiskArray (X, Y, 50, 20, Input, StartState, 1000); 
%     StartState = NextState;
% end;





% Simulate to compute the states of the disk. Initial state of all disks at
% -pi/4 angle. Figure out the disks whose states changed from pi/4, due to
% coupling and report it in State0 as a vector of 1 (changed) and 0 (no
% change)



% % Simulate to compute the states of the disk. Initial state of all disks at
% % pi/4 angle. Figure out the disks whose states changed from pi/4, due to
% % coupling and report it in State0 as a vector of 1 (changed) and 0 (no
% % change)
% fprintf(1, '\n With all disk vectors to be at pi/4');
% FinalState0 = NanodiskArray (X, Y, 50, 20, Input, 0.25*pi*ones(N,1), 3000);
% State0 = (abs(FinalState0 - pi/4) > 0.3); 
% figure; colormap(gray); 
% fprintf(1, '\n Coupled edges are shown with intensity 0 (black) and uncoupled with intensity 120 (gray).');
% imagesc(DrawSelectedLines (Lines, Layout, State0));
% 
% % Simulate to compute the states of the disk. Initial state of all disks at
% % -pi/4 angle. Figure out the disks whose states changed from pi/4, due to
% % coupling and report it in State0 as a vector of 1 (changed) and 0 (no
% % change)
% fprintf(1, '\n With all disk vectors to be at pi/4');
% Input(1) = FinalState0(1) + pi/2;
% FinalState1 = NanodiskArray (X, Y, 50, 20, Input, FinalState0);
% State1 = (abs(FinalState1 + pi/4)> 0.3);
% figure; colormap(gray); 
% fprintf(1, '\n Coupled edges are shown with intensity 0 (black) and uncoupled with intensity 120 (gray).');
% imagesc(DrawSelectedLines (Lines,  Layout, State0)); 

% % Simulate to compute the states of the disk. Initial state of all disks at
% % -3pi/4 angle. Figure out the disks whose states changed from pi/4, due to
% % coupling and report it in State0 as a vector of 1 (changed) and 0 (no
% % change)
% fprintf(1, '\n With all disk vectors to be at pi/4');
% State2 = (abs(NanodiskArray (X, Y, 50, 20, Input, -0.75*pi*ones(N,1)) + 0.75*pi) > 0.3);
% figure; colormap(gray); 
% fprintf(1, '\n Coupled edges are shown with intensity 0 (black) and uncoupled with intensity 120 (gray).');
% imagesc(DrawSelectedLines (Lines,  Layout, State0)); 
% 
% % Simulate to compute the states of the disk. Initial state of all disks at
% % 3pi/4 angle. Figure out the disks whose states changed from pi/4, due to
% % coupling and report it in State0 as a vector of 1 (changed) and 0 (no
% % change)
% fprintf(1, '\n With all disk vectors to be at pi/4');
% State3 = (abs(NanodiskArray (X, Y, 50, 20, Input, 0.75*pi*ones(N,1))  - 0.75*pi) > 0.3);
% figure; colormap(gray); 
% fprintf(1, '\n Coupled edges are shown with intensity 0 (black) and uncoupled with intensity 120 (gray).');
% imagesc(DrawSelectedLines (Lines,  Layout, State0));


%%%%--------------------------------------------

