% Geometry and Material Data
L1 = 15.0;
L2 = 12.0;
L3 = 5.0;

E = 29000; % Young's modulus
v = 0.3; % Poisson's ratio

A1 = 14; % Section area
Iz = 800; % Moment of inertia

% Joint coordinates
jointCoords = [0, 0;
               0, L2;
               L1, L2;
               L1 + L3, L2;
               L1, 0];

% Support data (node, x-restrained, y-restrained, rotation-restrained)
supportData = [1, 1, 1, 1; % Node 0 in Python → Node 1 in MATLAB
               5, 0, 1, 0]; % Node 4 in Python → Node 5 in MATLAB

% Young's moduli
youngsModuli = [E];

% Section properties (area, moment of inertia)
sectionProperties = [A1, Iz];

% Member information (start node, end node, section type)
MembersInfo = [1, 2, 1 1; 
               2, 3, 1 1;
               3, 4, 1 1;
               3, 5, 1 1];

% Joint load data (node, Fx, Fy, Mz)
jointLoadData = [1, 5, 0 0;
                 4, 0, -15 0]; 

jointDispData = [];

% Distributed load data (element, wx, wy)
distLoadData = [2, 0, -2]; % Applied on element 2 (in Python indexing)

disp('Data loaded!');
