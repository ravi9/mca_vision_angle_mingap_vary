function [Phi State TotalTime] = NanodiskArray (X, Y, r, h, Input, InitialState, H_ext)
%% 
%% This function simulates the magnetic states of a collection of nanodisks
%% based on LLG equation adapted for circular disks. It is assumed that each
%% disk is single domain which an inplane magnetization vector. It does not
%% allow for vortex states. This is based on the equations derived in 
%% 
%% Simulating collective magnetic dynamics in nanodisk arrays, A. J.
%% Bennetta and J. M. Xu, APPLIED PHYSICS LETTERS VOLUME 82, NUMBER 15, 2003
%% 
%% X and Y: location of disks in nanometers
%% 
%% r: radius of the disks
%% 
%% h: thickness of the disks
%% 
%% Input: input vector for each dot. length is same as the X and Y vector.
%% disks with input has the value of the phi(angle) and disks without input
%% has -10 as the corresponding value
%% 
%% InitialState: Initial angles of each of disks. The state in the Input
%% vector, if indicated, supercedes the state indicated in this vector
%% all spatial units are in r nanometers, i.e in terms of the dimensions of
%% the radius of the disks.
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

%% The cells are anti-ferromagnetically coupled. Did experiment with
%% ferromagnetic coupling but not successful -- hard to change coupled
%% arrays of cells using just changing one magnet.

gamma = 17.6e+9/(4*pi); %m/(A.s) gyromagnetic ratio
alpha = 0.5; %dimensionless sampling constant
M_o = 1e+5; %A/m magnetization
% r = 25; %nm radius of the disks
% h = 20; %nm thickness of the disks
% gap = 10; %nm
del_t = 1e-15; %seconds simulation time step
C = gamma * M_o * h*del_t/(4 * r * alpha);
C_e = 4*pi*(1e-7)*pi*M_o*M_o*r*h*h*(1e-8)/2 ;%units: 1e-19 J, distances normalized w.r.t radius, r
%note these above constants are the actual constants divided by r^3. This is
%done as we are working with normalized spatial distances, i.e. distances
%are in terms of number of radii, r.

global x22_y2 y22_x2 xy x2y2;
N = size(X, 1);
X = X/r; Y = Y/r; %normalize distances.
if (size(X, 1) == 1) X = X'; end; % turn into column vectors
if (size(Y, 1) == 1) Y = Y'; end;

x = repmat(X, 1, N); y = repmat(Y, 1, N); 
x = x - x'; y = y - y'; % distance of all cells from the i-th cell
x2 = x.^2; y2 = y.^2; % pre-compute the squares and the coefficients
x22_y2 = (2.*x2 - y2); 
y22_x2 = (2.*y2 - x2);
xy = x .* y;
x2y2 = (x2 + y2).^2.5+eps;


%DiskToDiskDistance = (2*r + gap)/r; 
% %% hexagonal lattice
% Y = [1:N]'*ones(1, N); X = ones(N, 1)*[1:N];
% X = [X X+0.5]; Y = [Y Y + 0.5];
% X = DiskToDiskDistance.*reshape (X, 2*N*N, 1); Y = DiskToDiskDistance.*reshape (Y, 2*N*N, 1);
% %% Rectagular lattice
% Y = [1:N]'*ones(1, N); X = ones(N, 1)*[1:N];
% X = DiskToDiskDistance.*reshape (X, N*N, 1); Y = DiskToDiskDistance.*reshape (Y, N*N, 1);
%% Line of cells
%  Y = DiskToDiskDistance.*ones(N, 1); X = DiskToDiskDistance.*[1:N]';
%% -------------------------------------


Phi = InitialState; %pi * (rand(N, 1)-0.5); % array of angles
Phi(find(Input ~= -10)) = Input(find(Input ~= -10));
FreeCells = find(Input == -10); %FreeCells are the indices of the non-input cells
%% First tried matlab ode code did not work. Need to study the details
%% about the paramters needs more careful.
% Tspan = [0 1e4*del_t]; % Solve from t=1 to t=5
% options = odeset('InitialStep', del_t, 'MaxStep', 5*del_t);
% [T fPhi] = ode45(@(t,y) dPhi (t, Phi, X, Y, N, C, C_e),Tspan, Phi, options); % Solve ODE
% 
% for i=1:size(T,1)
%      subplot(1,2,1) ;plot(X, Y, 'o'); hold on; quiver(X, Y, sin(fPhi(i,:))', cos(fPhi(i,:))' ,0); hold off;
%      Energy = Hamiltonian (fPhi(i,:)', X, Y, N, C_e);
%      subplot(1,2,2) ; plot(T(i), Energy, 'o'); hold on;
%      pause(0.01);
% end;

%Normalized distance between neighboring disks
Steps = 1; LastChange = 0; 
[DelPhi Energy] = dPhi (Steps, Phi, N, C, C_e, H_ext);
NewChange = max(abs(DelPhi));
%while (Steps < NumSteps)
while (abs(LastChange - NewChange) > 1e-11)
    Steps = Steps + 1;
    % Update only the non input cells
    Phi (FreeCells)  = Phi (FreeCells) + DelPhi (FreeCells);
    LastChange = NewChange;
    [DelPhi Energy] = dPhi (Steps, Phi, N, C, C_e, H_ext);
    NewChange = max(abs(DelPhi)); 
    %
    %if (rem(Steps, 100) == 0), Phi = Phi + (rand(N, 1)- 0.5); end; %
    %random perturbation at regular intervals does not make much
    %difference.
%     if (rem(Steps, 200) == 1)
%        % State = InferState (N, Phi);
%         subplot(2,1,1); plot(X, Y, 'ro'); hold on; plot(X(FreeCells), Y(FreeCells), 'o'); 
% %         plot(X(find(abs(State)>0.5)), Y(find(abs(State)>0.5)), 'o');  hold on; 
% %         plot(X(find(abs(State)<=0.5)), Y(find(abs(State)<=0.5)), 'ro'); hold on;
% %         axis equal; 
% %         for (i=1:N)
% %             rectangle('Position',[X(i)-1,Y(i)-1,2,2], 'Curvature',[1,1]); 
% %             hold on; 
% %         end;
%         quiver(X, Y, cos(Phi), sin(Phi) ,0.2); hold off;
%         
%         subplot(2,1,2) ; plot(Steps, Energy, '.'); hold on;
%         fprintf(1, '\n %d %f %f', Steps, LastChange, NewChange);
%         pause(0.0000001);
%     end;
end;


State = InferState (N, Phi);
TotalTime = del_t *Steps;

%% Compute the change in phi based on LLG evoluation. Need to monitor the
%% values of H_x, H_y and dPhi -- they should not be too large.
   
function [dPhi_dt Energy]= dPhi (t, Phi, N, C, C_e, H_external)

Energy = 0;

global x22_y2 y22_x2 xy x2y2; %precomputed coefficients that are dependent on location of cells.
if (size(Phi, 1) == 1) Phi = Phi'; end;
cosPhi = repmat(cos(Phi), 1, N); % precompute the sin and cos
sinPhi = repmat(sin(Phi), 1, N);

% H_x = (cosPhi.* (2.*x2 - y2) + 3 * sinPhi.*(x .* y));
% H_x = sum(H_x./((x2 + y2).^2.5+eps)) + H_external(1);
% 
% H_y = (sinPhi .* (2.*y2 - x2) + 3 * cosPhi.*(x .* y));
% H_y = sum(H_y./((x2 + y2).^2.5+eps)) + H_external(2);

H_x = (cosPhi.* (x22_y2) + 3 * sinPhi.*(xy));
H_x = sum(H_x./x2y2) + H_external(1);

H_y = (sinPhi .* (y22_x2) + 3 * cosPhi.*(xy));
H_y = sum(H_y./x2y2) + H_external(2);

dPhi_dt = C * (sinPhi(:,1)'.*H_x - cosPhi(:,1)'.*H_y);
Energy = sum(C_e*(cosPhi(:,1)'.*H_x + sinPhi(:,1)'.*H_y));

% for i=1:N
%     x = X - X(i); y = Y - Y(i); % distance of all cells from the i-th cell
%     
%     H_x_i = (cos(Phi) .* (2.*x.^2 - y.^2) + 3 * sin(Phi).*(x .* y));
% 
%     H_x = sum(H_x_i./((x.^2 + y.^2).^2.5+eps)) + H_external(1);
%     
%     H_y_i = (sin(Phi) .* (2.*y.^2 - x.^2) + 3 * cos(Phi).*(x .* y));
%     H_y = sum(H_y_i./((x.^2 + y.^2).^2.5+eps)) + H_external(2);
%     
%     dPhi_dt (i) = C * (sin(Phi(i)) * H_x - cos (Phi(i))*H_y);
%     
%     Energy = Energy + C_e*(cos(Phi(i)) * H_x + sin (Phi(i))*H_y);
%     %subplot(2,2,3); plot3(X(i), Y(i), H_x); hold on; subplot(2,2,4); plot3(X(i), Y(i), H_y); hold on;
%     %fprintf(1, '\n H_x = %f, H_y = %f, C = %f, DelPhi(%d) = %f', H_x, H_y, C, i, dPhi_dt(i));
% end;
%fprintf(1, '\n H_x = %f, H_y = %f', max(H_x), max(H_y));
dPhi_dt = dPhi_dt';

%%
function  State = InferState (N, Phi)
%% Figuring out the states as the rank 1 decomposition of the coupling
%% constants (v_i v_j) -- see paper for notation.
%% The idea does not work in practice.

global x22_y2 y22_x2 xy x2y2; %precomputed coefficients that are dependent on location of cells.

if (size(Phi, 1) == 1) Phi = Phi'; end;

cosPhi = repmat(cos(Phi), 1, N); % precompute the sin and cos
sinPhi = repmat(sin(Phi), 1, N);

d = x2y2.^(2/5);
H_x = (cosPhi.* (x22_y2) + 3 * sinPhi.*(xy));

H_y = (sinPhi .* (y22_x2) + 3 * cosPhi.*(xy));


xi_xj = (cosPhi'.*H_x + sinPhi'.*H_y)./d;

xi_xj = xi_xj + diag(2*ones(1, N));
[State Lambda] = eigs(xi_xj, 1);

State = State'./max(abs(State));



%%  ------------------------------another options is to optimize the
%%  following Hamiltonian by some other means
function Energy = Hamiltonian (Phi, X, Y, N, C_e)
Energy = 0;
for i=1:N
    x = X - X(i); y = Y - Y(i); %distance from i-th cell.
    H_x_i = (cos(Phi) .* (2.*x.^2 - y.^2) + 3 * sin(Phi).*(x .* y));
    H_x = sum(H_x_i./((x.^2 + y.^2).^2.5+eps));
    H_y_i = (sin(Phi) .* (2.*y.^2 - x.^2) + 3 * cos(Phi).*(x .* y));
    H_y = sum(H_y_i./((x.^2 + y.^2).^2.5+eps));
    Energy = Energy + C_e*(cos(Phi(i)) * H_x + sin (Phi(i))*H_y);
    %fprintf(1, '\n H_x = %f, H_y = %f, C = %f, DelPhi(%d) = %f', H_x, H_y, C, i, DelPhi(i));
end;
