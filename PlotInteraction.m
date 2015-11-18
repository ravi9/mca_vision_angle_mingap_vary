function y = PlotInteraction (r, h)

% A service standalone function to plot the interaction
% function between two circular nanometers of radius r
% and height h. This is used to visualize the interaction
% form.

gamma = 17.6e+9/(4*pi); %m/A/s gyromagnetic ratio
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

x = [-5:0.1:5]; x = x'*ones(size(x));
y = [-5:0.1:5]; y = ones(size(y))'*y;


for Phi_j = 0:0.1:2*pi

    Phi_i = 0;

    H_x = (cos(Phi_j).* (2.*x.^2 - y.^2) + 3 * sin(Phi_j).*(x .* y));
    H_x = (H_x./((x.^2 + y.^2).^2.5+eps));

    H_y = (sin(Phi_j) .* (2.*y.^2 - x.^2) + 3 * cos(Phi_j).*(x .* y));
    H_y = (H_y./((x.^2 + y.^2).^2.5+eps));


    Energy = C_e*(cos(Phi_i).* H_x + sin (Phi_i).*H_y);

    Energy = (x.^2 + y.^2 > 1).*Energy; %mask out the disk area.
    mesh (x, y, Energy);
    pause(1);
end