% Calculation of the bridge by the Finite Element Method
% Example of calculation using Finite Element Method with Euler Bernoulli 
% Beams implemented in single GNU Octave file.
% Author: Milos D. Petrasinovic <mpetrasinovic@mas.bg.ac.rs>
% Structural Analysis of Flying Vehicles
% Faculty of Mechanical Engineering, University of Belgrade
% Department of Aerospace Engineering, Flying structures
% https://vazmfb.com
% Belgrade, 2022
%
% ---------------
%
% Copyright (C) 2022 Milos Petrasinovic <info@vazmfb.com>
%  
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as 
% published by the Free Software Foundation, either version 3 of the 
% License, or (at your option) any later version.
%   
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%   
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
% ---------------
clear('all'), clc, close('all'), tic
disp([' --- ' mfilename ' --- ']);

% - Input parameters
% Bridge dimensions
H = 3000; % [mm] bridge height
B = 2500; % [mm] bridge width
L = 12000; % [mm] bridge length

% Characteristics of cross section
Ds = 60; % [mm] outer diameter of the pipe
t = 5; % [mm] Pipe wall thickness
% Material characteristics
E = 210e9; % [Pa] Elastic modulus
ni = 0.3; % [-] Poisson's ratio

% Coordinates of nodes
% nc(node, 1:3) = [x, y, z];
% Left side of the bridge
nc = [0, 0, 0; L/4, 0, 0; L/2, 0, 0; L*3/4, 0, 0; L, 0, 0; 
  L/4, 0, H; L/2, 0, H;L*3/4, 0, H;];
% Right side of the bridge
nc = [nc; nc+[0, B, 0];]; 

% Nodes of elements
% en(element, :) = [node1, node2];
% Left side of the bridge
en = [1, 2; 1, 6; 6, 2; 2, 7; 6, 7; 2, 3; 
  3, 7; 7, 8; 3, 4; 7, 4; 4, 8; 4, 5; 8, 5;]; 
% Right side of the bridge
en = [en; en+[8, 8]]; 
% Spreads
en = [en; (1:8).', (9:16).']; 

% Nodes connections
% ec(i, 1:8) = [node1, node2, Tx, Ty, Tz, Rx, Ry, Rz];
% 0 - free
% 1 - constrained
ec = [];

% Boundary conditions
% bc(node, 1:6) = [Tx, Ty, Tz, Rx, Ry, Rz];
% 0 - free
% 1 - constrained
bc([1, 9], :) = repmat([1, 1, 1, 0, 0, 0], 2, 1);
bc([5, 13], :) = repmat([0, 1, 1, 0, 0, 0], 2, 1);

% External loads
% F(node, 1:6) = [Fx[N], Fy[N], Fz[N], Mx[Nm], My[Nm], Mz[Nm]];
F(3, :) = [0, 0, -20000, 0, 0, 0];
F(11, :) = [2000, 0, -20000, 0, 0, 0];

% Additional variables
s = [1000, 200, 300, 3*10^-2, 1, 1000, 0.2, 100]; % Scaling coefficients

% - Defining a stiffness matrix
d2r = 180/pi; % degrees to radians
ne = size(en, 1); % number of elements
nn = size(nc, 1); % number of nodes
nec = size(ec, 1); % number of connected nodes
dof = 6*nn; % number of degrees of freedom
D = zeros(dof, 1); % displacements
bc(end+1:nn, :) = zeros(length(size(bc, 1)+1:nn), 6); % boundary conditions
F(end+1:nn, :) = zeros(length(size(F, 1)+1:nn), 6); % external loads
F = reshape(F.', [], 1); 
K = zeros(dof, dof); % stiffness matrix

% Degrees of freedom of nodes
rdof = find(bc')'; % constrained degrees of freedom
fdof = find(~bc')'; % free degrees of freedom

% Connecting nodes
Cp = eye(6);
Np = 0; % number of additional degrees of freedom
if(nec)
    Np = size(ec(ec(:, 3:8)>0), 1);
end
lmn = (1:Np)+dof; % additional degrees of freedom

Q = zeros(Np, 1); % vector for nodal load expansion
C = zeros(Np, dof); % matrix for stiffness matrix extension 
j = 0;
for i = 1:nec
  Npi = size(ec(ec(i, 3:8)>0), 2);
  C(j+1:j+Npi, (ec(i, 1)-1)*6+(1:6)) = Cp(ec(i, 3:8)>0, :);
  C(j+1:j+Npi, (ec(i, 2)-1)*6+(1:6)) = -Cp(ec(i, 3:8)>0, :);
  j = j+Npi;
end

% Material and cross-sectional characteristics
G = E/(2*(1+ni)); % Shear modulus

% Cross-sectional characteristics (circular tube)
Ds = Ds/s(1); % [m]
t = t/s(1); % [m]
A = (Ds^2-(Ds-2*t)^2)*pi/4; % [m^2] Cross-sectional area
Iy = (Ds^4-(Ds-2*t)^4)*pi/64; % [m^4] Axial moment of inertia for the y-axis
Iz = Iy; % [m^4] Axial moment of inertia for the z-axis
J = (Ds^4-(Ds-2*t)^4)*pi/32; % [m^4] Polar moment of inertia

% Determination of stiffness matrix
nc = nc/s(1); % [m]
for i=1:ne 
    j = en(i, :);       
    edof = [6*j(1)-5 6*j(1)-4 6*j(1)-3 ... % degrees of freedom of the element
           6*j(1)-2 6*j(1)-1 6*j(1)...
           6*j(2)-5 6*j(2)-4 6*j(2)-3 ...
           6*j(2)-2 6*j(2)-1 6*j(2)]; 
    L_e = sqrt((nc(j(2), 1)-nc(j(1), 1))*(nc(j(2), 1)-... % element length
        nc(j(1), 1))+(nc(j(2), 2)-nc(j(1), 2))*(nc(j(2), 2)-...
        nc(j(1), 2))+(nc(j(2), 3)-nc(j(1), 3))*(nc(j(2), 3)-nc(j(1), 3)));
    ki = [E*A/L_e, 12*E*Iz/(L_e^3), 6*E*Iz/(L_e^2), 4*E*Iz/L_e, ...
        2*E*Iz/L_e, 12*E*Iy/(L_e^3), 6*E*Iy/(L_e^2), ...
        4*E*Iy/L_e, 2*E*Iy/L_e, G*J/L_e];
    a = diag([ki(1), ki(2), ki(6)]);
    b(2, 3) = ki(3); b(3, 2) = -ki(7);
    c = diag([ki(10), ki(8), ki(4)]);
    d = diag([-ki(10), ki(9), ki(5)]);
    
    % Element stiffness matrix in the local coordinate system of the element
    k_e = [a, b, -a, b; b.', c, b, d; -a.', b.', a, -b; b.', d.', -b.', c];
  
    if nc(j(1), 1) == nc(j(2), 1) && nc(j(1), 2) == nc(j(2), 2)
        if nc(j(2), 3) > nc(j(1), 3)
            r = [0, 0, 1; 0, 1, 0; -1, 0, 0];
        else
            r = [0, 0, -1; 0, 1, 0; 1, 0, 0];
        end
    else
        CXx = (nc(j(2), 1)-nc(j(1), 1))/L_e;
        CYx = (nc(j(2), 2)-nc(j(1), 2))/L_e;
        CZx = (nc(j(2), 3)-nc(j(1), 3))/L_e;
        i_eXY = sqrt(CXx^2+CYx^2);
        CXy = -CYx/i_eXY;
        CYy = CXx/i_eXY;
        CZy = 0;
        CXz = -CXx*CZx/i_eXY;
        CYz = -CYx*CZx/i_eXY;
        CZz = i_eXY;
        r = [CXx, CYx, CZx; CXy, CYy, CZy; CXz, CYz, CZz];
    end
    
    T = [r, zeros(3, 9);  % Transformation matrix
          zeros(3), r, zeros(3, 6);
          zeros(3, 6), r, zeros(3); 
          zeros(3, 9), r];
        
    K(edof, edof) = K(edof, edof)+T.'*k_e*T; % Element stiffness matrix
end  

% - Solving FEM
Kp = [K, C.'; C, zeros(Np, Np)]; % extended stiffness matrix 
Fp = [F; Q]; % extended nodal load vector

% Solving the reduced finite element equation
D1 = Kp([fdof, lmn], [fdof, lmn])\Fp([fdof, lmn]); 

D([fdof, lmn]) = D1;
Rp = Kp*D; % extended reaction vector
lam = D((end+1-Np):end); % Lagrange multipliers
R = Rp(rdof); % reactions
D = reshape(D(1:(end-Np)), 3, []); % node displacement
DT = D(:, 1:2:end); % translation 
DR = D(:, 2:2:end); % angular displacement

% - Display of results
% Display of displacement values and reactions
disp(' -------------------- ');
disp(' Displacements')
i = 1:dof/2;
DTv = [reshape(repmat(1:nn, 3, 1), 1, []); ...
  reshape(repmat(1:3, nn, 1).', [], 1).'; DT(i)*s(1)];
disp(' Node | Component | Displacement [mm]');
fprintf(' %3d | %3d | %14.10f\n', DTv);

DRv = [reshape(repmat(1:nn, 3, 1), 1, []); ...
  reshape(repmat(4:6, nn, 1).', [], 1).'; DR(i)/d2r];
disp(' Node | Component | Displacement [deg]');
fprintf(' %3d | %3d | %14.10f\n', DRv);

% Display of reactions
disp(' -------------------- ');
disp(' Reactions')
nrdof = mod(rdof, 6); % constrained degrees of freedom
nrdof(nrdof == 0) = ones(length(find(nrdof == 0)), 1)*6;
nr = (rdof-nrdof)/6+1; % node with constraints
FTv = [nr(nrdof < 4).', nrdof(nrdof < 4).', R(nrdof < 4)].';
disp(' Node | Component | Reaction [N]');
fprintf(' %3d | %3d | %14.10f\n', FTv);

FRv = [nr(nrdof > 3).', nrdof(nrdof > 3).', R(nrdof > 3)].';
disp(' Node | Component | Reaction [Nm]');
fprintf(' %3d | %3d | %14.10f\n', FRv);

% - Display of model
% Display of initial model with loads
disp(' -------------------- ');
disp(' Display of initial model with loads... ');
drawArrow = @(x, y, z, varargin) quiver3(x(1), y(1), z(1), ...
    x(2)+10^-5, y(2)+10^-5, z(2)+10^-5, 0, varargin{:});   
drawArrowMarker = @(x, y, z, m, varargin) [plot3([x(1); x(1)+x(2)], ...
  [y(1); y(1)+y(2)], [z(1); z(1)+z(2)], '-', varargin{:}), ...
  plot3(x(1)+x(2), y(1)+y(2), z(1)+z(2), m, varargin{:})];    

nc = nc*s(1);
bms = [[nc(en(:, 1), 1), nc(en(:, 2), 1)], ...
    [nc(en(:, 1), 2), nc(en(:, 2), 2)], ...
    [nc(en(:, 1), 3), nc(en(:, 2), 3)]]; % beams
  
figure(1);
box on, grid on, hold on
c = get(gca, 'colororder');
plot3(bms(:, 1:2).', bms(:, 3:4).', bms(:, 5:6).', ...
    'LineWidth', 2, 'Color', c(1, :));
plot3(nc(:, 1), nc(:, 2), nc(:, 3), 'o', 'LineWidth', 2, ...
    'MarkerFaceColor', c(2, :), 'Color', c(2, :));
text(nc(:, 1)+s(8), nc(:, 2)+s(8), nc(:, 3), ...
   num2str((1:nn).'), 'FontSize', 18, 'Color', c(2, :));
rotate3d, view(45, 15), axis equal, set(gca, 'FontSize', 18);
xlim([-s(6), L+s(6)]), ylim([-s(6), B+s(6)]), zlim([-s(6), H+s(6)])

% Display of boundary conditions
for i=1:length(rdof)
  p = zeros(1, 3);
  if(nrdof(i) < 4) 
    p(nrdof(i)) = 1*s(3);
    drawArrowMarker([nc(nr(i), 1), p(1)], [nc(nr(i), 2), p(2)], ...
      [nc(nr(i), 3), p(3)], '+', 'Color', 'b', 'LineWidth', 2, ...
      'MarkerFaceColor', 'b');
  else
    p(nrdof(i)-3) = -1*s(3);
    drawArrowMarker([nc(nr(i), 1), p(1)], [nc(nr(i), 2), p(2)], ...
      [nc(nr(i), 3), p(3)], 'o', 'Color', 'r', 'LineWidth', 2, ...
      'MarkerFaceColor', 'r');
  end
end

% Display of external forces
for i=1:nn
  p = zeros(1, 3);
  if(norm(F((i-1)*6+1:(i-1)*6+3)) > 0)
    p = (F((i-1)*6+1:(i-1)*6+3))*s(4);
    drawArrow([nc(i, 1), p(1)], [nc(i, 2), p(2)], ...
      [nc(i, 3), p(3)], 'Color', 'b', 'LineWidth', 2, ...
      'MarkerFaceColor', 'b', 'MaxHeadSize', s(7));
  end
  if(norm(F((i-1)*6+4:(i-1)*6+6)) > 0)
    p = (F((i-1)*6+4:(i-1)*6+6))*s(5);
    drawArrow([nc(i, 1), p(1)], [nc(i, 2), p(2)], ...
      [nc(i, 3), p(3)], 'Color', 'r', 'LineWidth', 2, ...
      'MarkerFaceColor', 'r', 'MaxHeadSize', s(7));
  end
end

% Display of deformed model with loads
disp(' Display of deformed model with loads... ');
ncd = nc+(DT.'.*s(1)).*s(2);
bmsd = [[ncd(en(:, 1), 1), ncd(en(:, 2), 1)], ...
    [ncd(en(:, 1), 2), ncd(en(:, 2), 2)], ...
    [ncd(en(:, 1), 3), ncd(en(:, 2), 3)]]; % beams

figure(2);
box on, grid on, hold on
plot3(bms(:, 1:2).', bms(:, 3:4).', bms(:, 5:6).', '--', ...
    'LineWidth', 1, 'Color', c(1, :));
plot3(bmsd(:, 1:2).', bmsd(:, 3:4).', bmsd(:, 5:6).', ...
    'LineWidth', 2, 'Color', c(2, :));
plot3(ncd(:, 1), ncd(:, 2), ncd(:, 3), 'o', 'LineWidth', 2, ...
    'MarkerFaceColor', c(2, :), 'Color', c(2, :));
rotate3d, view(45, 15), axis equal, set(gca, 'FontSize', 18);
xlim([-s(6), L+s(6)]), ylim([-s(6), B+s(6)]), zlim([-s(6), H+s(6)])

% - End of program
disp(' -------------------- ');
disp(' They belong to all and treat all alike; they are useful,');
disp(' always built for a purpose, at a spot where most human');
disp(' needs entwine; they are more durable than other buildings');
disp(' and serve no secret or evil purpose.');
disp('  - Bridges, Ivo Andric');
disp(' -------------------- ');
disp(' The program was successfully executed... ');
disp([' Execution time: ' num2str(toc, '%.2f') ' seconds']);
disp(' -------------------- ');