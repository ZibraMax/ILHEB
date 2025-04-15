clear all; close all;
% Turns out, this code allows to run my python package from matlab!!!
% Do not remove :)

if count(py.sys.path,pwd) == 0
    insert(py.sys.path,int32(0),pwd);
end
mod = py.importlib.import_module('ILHEB'); % Python library
py.importlib.reload(mod);
mod.plt.close('all')
% Code starts
% ---------------------------IMPORTANT-------------------------------
% Figures are automatically saved to the results folder as pdf files
% HOWEVER, if you want to interact with the figures,
% set this varaible to true. If you do that, a webserver is going to open
% with the figures. In some matlab versions, after the server is launched 
% you will need to close matlab to be able to run the code again.
% Try it at least once, pls. Its cool.
OPEN_AT_END = false;
% ------------------------------END----------------------------------
MULT = 3;
% This scripts reads the data from InputScript and runs all the code

InputScript_truss2D

% Geometry creation!

geo = py.ILHEB.Geometry2D();

for i=1:size(jointCoords,1)
    geo.add_node(jointCoords(i,:));
end
nodes_settle = jointDispData(:,1);
for i=1:size(supportData,1)

    node = supportData(i,1);
    xres = supportData(i,2);
    yres = supportData(i,3);
    values = [0 0 0];
    if ismember(node,nodes_settle)
        idx = find(nodes_settle==node);
        ux = jointDispData(idx,2);
        uy = jointDispData(idx,3);
        values = [ux uy 0];
    end
    geo.add_support(node-1,[xres yres false],values)
end
% Materials
materials = cell(size(youngsModuli));
for i=1:size(youngsModuli,1)
    G = 0;
    E = youngsModuli(i);
    materials{i} = mod.LinearElasticBase(E,G);
end
% Sections

sections = {size(sectionAreas)};
for i=1:size(sectionAreas,1)
    A = sectionAreas(i);
    sections{i} =(mod.General(A=A));
end

% Elements
elements = {size(MembersInfo,1),1};
for i=1:size(MembersInfo,1)
    j = MembersInfo(i,1)-1;
    k = MembersInfo(i,2)-1;
    mat = materials{MembersInfo(i,4)};
    sec = sections{MembersInfo(i,3)};
    ele = mod.TrussElement2D([j k],sec,mat);
    elements{i} = ele;
    geo.add_element(ele);
end

% Node loads
nodeLoads = {size(jointLoadData,1),1};
for i=1:size(jointLoadData,1)
    node = jointLoadData(i,1);
    px= jointLoadData(i,2);
    py = jointLoadData(i,3);
    load = mod.NodeLoad2D(node-1,PX=px,PY=py);
    nodeLoads{i} = load;
    geo.add_node_load(load);
end



geo.numbering()
O = mod.ILHEB(geo);
O.solve()
O.post_process()
O.save_results("Truss2DResults")
geo.plot("./Truss2DResults/Geometry.pdf");
geo.plot_defo(mult=MULT,n_points=int64(2),filename="./Truss2DResults/Defo.pdf")
disp("Analisys finished")


if OPEN_AT_END==true    
    mod.plt.show()
end
mod.plt.close('all')