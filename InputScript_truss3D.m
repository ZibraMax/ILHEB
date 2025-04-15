% Geometry and Material Data
E = 70; % Young's modulus

A = 10000; % Section area

% Joint coordinates (x, y, z)
jointCoords = [0, 0, 8;
               -4, -2, 0;
               4, -2, 0;
               4, 2, 0;
               -4, 2, 0];

% Support data (node, x-restrained, y-restrained, z-restrained, 
%               rotation-x, rotation-y, rotation-z)
supportData = [2, 0, 0, 1; % Node 1 in Python → Node 2 in MATLAB
               3, 1, 1, 1; % Node 2 in Python → Node 3 in MATLAB
               4, 1, 0, 0; % Node 3 in Python → Node 4 in MATLAB
               5, 0, 0, 1]; % Node 4 in Python → Node 5 in MATLAB

% Young's modulus and shear modulus
youngsModuli = [E];

% Section properties (area)
sectionAreas = [A];


% Member information (start node, end node, section type)
MembersInfo = [2, 3, 1 1;
               3, 4, 1 1;
               4, 5, 1 1;
               5, 2, 1 1;
               5, 3, 1 1;
               1, 2, 1 1;
               1, 3, 1 1;
               1, 4, 1 1;
               1, 5, 1 1];

% Joint load data (node, Fx, Fy, Fz)
jointLoadData = [1, 30, -40, -60]; % Node 0 in Python → Node 1 in MATLAB
jointDispData = [];
disp('Data loaded!');
