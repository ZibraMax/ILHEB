% Geometry and Material Data
E = 200e6; % Young's modulus (Pa)
G = 80e6; % Shear modulus (Pa)

A = 0.001; % Cross-sectional area (m²)
J = 2e-3; % Torsional constant (m⁴)
Iz = 1e-3; % Second moment of area about z-axis (m⁴)
Iy = 1e-3; % Second moment of area about y-axis (m⁴)

% Joint coordinates (x, y, z)
jointCoords = [0, 0, 5;
               4, 0, 5;
               4, 0, 0;
               0, -3, 0;
               0, 3, 0];

% Support data (node, x-restrained, y-restrained, z-restrained, 
%               rotation-x, rotation-y, rotation-z)
supportData = [3, 1, 1, 1, 1, 1, 1; % Node 2 in Python → Node 3 in MATLAB
               4, 1, 1, 1, 1, 1, 1; % Node 3 in Python → Node 4 in MATLAB
               5, 1, 1, 1, 0, 0, 0]; % Node 4 in Python → Node 5 in MATLAB

% Material properties (E, G)
materialProperties = [E, G];

% Section properties (A, J, Iz, Iy)
sectionProperties = [A, J, Iz, Iy];

% Member information (start node, end node, section type)
MembersInfo = [1, 2, 1 1;
               2, 3, 1 1;
               4, 1, 1 1;
               1, 5, 1 1];

% Joint load data (node, Fx, Fy, Fz,Mx, My, Mz)
jointLoadData = [1, 0, 0, -120 0 0 0; % Node 0 in Python → Node 1 in MATLAB
                 2, 60, 0, 0 0 0 0];  % Node 1 in Python → Node 2 in MATLAB
distLoadData = [];
jointDispData = [];
disp('Data loaded!');
