L = 1;
A1 = 1;
A2 = 1;
E1 = 1;
E2 = 1;
P = 1;
delta = 0;

% Joint coordinates
jointCoords = [0 0;
               L 0;
               0 L;
               L L];

% Support data (node, x-restrained, y-restrained)
supportData = [1 1 1;
               2 0 1;
               4 1 0];

% Young's moduli
youngsModuli = [E1;
                E2];

% Section areas
sectionAreas = [A1;
                A2];

% Number of cross-section types
NCT = 2;

% Joint load data (node, Fx, Fy)
jointLoadData = [3 -P -2*P];

% Member information (start node, end node, section type)
MembersInfo = [1 2 1 1; % member 1
               1 3 1 1; % member 2
               3 4 2 2; % member 3
               2 4 2 2; % member 4
               1 4 1 1]; % member 5

% Joint displacement data (node, Ux, Uy)
jointDispData = [1 0 -delta];
disp("Data loaded!");