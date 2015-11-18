function layout = AffinityMDS(Affinity)
% Assumes that it is passed an affinity matrix
% The output "layout" is an array of structures with the layout of each
% component,i.e. layout.X and layout.Y, the indices of the original features
% that participate in the component,

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

layout = [];

%% Find Connected Components
N = size(Affinity, 1);

Distances = ChangeToDistances(Affinity);
Distances = RestoreMetricProperty (Distances);


% options.dims = [1:4]; %(row) vector of embedding dimensionalities to use
% %                        (1:10 = default)
% %options.comp = which connected component to embed, if more than one. 
% %                        (1 = largest (default), 2 = second largest, ...)
% options.display = 1; %plot residual variance and 2-D embedding?
% %                        (1 = yes (default), 0 = no)
% options.overlay = 1; %overlay graph on 2-D embedding?  
% %                        (1 = yes (default), 0 = no)
% options.verbose = 1; %display progress reports? 
% %                        (1 = yes (default), 0 = no)
% [C, R, E] = Isomap(Distances, 'k', 10, options);
% layout.X = C.coords{2}(1,:); layout.Y = C.coords{2}(2,:); layout.indices
% = C.index;
% fprintf(1, '\n Embedded %d points out of %d points', length(C.index), N);

[X Y] = MDS (Distances); layout.X = X; layout.Y = Y; layout.indices = [1:N];
    

%% Perform MDS on and return first two coordinates
function [X Y] = MDS (D)
% opts = statset('display', 'iter');
% [Y,stress] = mdscale(D,2, 'criterion','stress', 'Options', opts );
Y = cmdscale(D);
X = Y(:,1); Y  = Y(:, 2); 

% subplot(2,2,1); plot(X, Y, 'o'); 

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
fprintf(1, '\n Restoring metrix property by adding %f', c);
y = D+c;
for i=1:N, y(i,i) = 0; end;

return;
%% Change Affinities to Distances. Affinities could be larger than one.

function D = ChangeToDistances (A)

maxA = max(max(A));
minA = min(min(A));

if ((maxA > 1.0) || (minA < 0)) 
    fprintf(1, '\n Values of the affinity matrix are beyond  0 to 1 range. [%f %f]', minA, maxA);
    maxA = maxA + eps; minA = eps;
    A = (A-minA)/(maxA - minA);
end;
%D = (A+eps).^-(1.0/3);
D = (-log(A)); %maps affinity in the range 0 to 1 into the range infinity and zero -- is this
% linear manifold approximation? measure of "repulsion" if 0 < A < 1.

sum(sum(D<0))

%D = (D > 0).*D; %% to remove very small negative created values due to numerical issues
D = (D).^(1/3);

% D = (log(A+1)).^(1/3);

fprintf(1, '\n min and max values of the distance matrix [%f %f]', min(min(D)), max(max(D)));
N = size(D, 1);
for (i=1:N) D(i, i) = 0.0; end; %% Restore diagonal to zero distances
D = 0.5*(D+D'); % make the matrix symmetric.

