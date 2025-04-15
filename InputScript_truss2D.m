% Geometry and Material Data
L1 = 12.0;
L2 = 8.0;
settlement = -0.5;

E = 29000; % Young's modulus
v = 0.3; % Poisson's ratio

A1 = 3.5; % Section area 1
A2 = 2.5; % Section area 2

% Joint coordinates
jointCoords = [0, 0;
               L1, L2;
               2*L1, L2;
               3*L1, L2;
               4*L1, 0;
               L1, 0;
               2*L1, 0;
               3*L1, 0];

% Support data (node, x-restrained, y-restrained)
supportData = [1, 1, 1; % Node 0 in Python → Node 1 in MATLAB
               5, 0, 1]; % Node 4 in Python → Node 5 in MATLAB

% Prescribed displacements (node, Ux, Uy)
jointDispData = [5, 0, settlement];

% Young's moduli
youngsModuli = [E;
                E];

% Section areas
sectionAreas = [A1;
                A2];

% Number of cross-section types
NCT = 2;

% Member information (start node, end node, section type)
MembersInfo = [1, 2, 1 1; 
               2, 3, 1 1;
               3, 4, 1 1;
               4, 5, 1 1;
               1, 6, 2 1;
               6, 7, 2 1;
               7, 8, 2 1;
               8, 5, 2 1;
               6, 2, 1 1;
               7, 2, 1 1;
               7, 3, 1 1; 
               7, 4, 1 1;
               8, 4, 1 1];

% Joint load data (node, Fx, Fy)
jointLoadData = [6, 0, -25; % Node 5 in Python → Node 6 in MATLAB
                 7, 0, -20;
                 8, 0, -20];

disp('Data loaded!');
