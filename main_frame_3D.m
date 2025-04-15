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

MULT = 500;

% This scripts reads the data from InputScript and runs all the code

InputScript_frame3D

% Geometry creation!

geo = py.ILHEB.Geometry3D();

for i=1:size(jointCoords,1)
    geo.add_node(jointCoords(i,:));
end
nodes_settle  = [];
if size(jointDispData) >0
    nodes_settle = jointDispData(:,1);
end
for i=1:size(supportData,1)

    node = supportData(i,1);
    xres = supportData(i,2);
    yres = supportData(i,3);
    zres = supportData(i,4);
    rxres = supportData(i,5);
    ryres = supportData(i,6);
    rzres = supportData(i,7);
    values = [0 0 0 0 0 0];
    if ismember(node,nodes_settle)
        idx = find(nodes_settle==node);
        ux = jointDispData(idx,2);
        uy = jointDispData(idx,3);
        uz = jointDispData(idx,4);
        rx = jointDispData(idx,5);
        ry = jointDispData(idx,6);
        rz = jointDispData(idx,7);
        values = [ux uy uz rx ry rz];
    end
    geo.add_support(node-1,[xres yres zres rxres ryres rzres],values)
end
% Materials
materials = cell(size(materialProperties));
for i=1:size(materialProperties,1)
    
    E = materialProperties(i,1);
    G = materialProperties(i,2);
    materials{i} = mod.LinearElasticBase(E,G);
end
% Sections

sections = {size(sectionProperties,1),1};
for i=1:size(sectionProperties,1)
    A = sectionProperties(i,1);
    J = sectionProperties(i,2);
    Iz = sectionProperties(i,3);
    Iy = sectionProperties(i,4);
    sections{i} =(mod.General(A=A,J=J,Iy=Iy,Iz=Iz));
end

% Elements
elements = {size(MembersInfo,1),1};
for i=1:size(MembersInfo,1)
    j = MembersInfo(i,1)-1;
    k = MembersInfo(i,2)-1;
    mat = materials{MembersInfo(i,4)};
    sec = sections{MembersInfo(i,3)};
    ele = mod.FrameElement3D([j k],sec,mat);
    elements{i} = ele;
    geo.add_element(ele);
end

% Node loads
if size(jointLoadData,1) > 0
    nodeLoads = {size(jointLoadData,1),1};
    for i=1:size(jointLoadData,1)
        node = jointLoadData(i,1);
        px= jointLoadData(i,2);
        py = jointLoadData(i,3);
        pz = jointLoadData(i,4);
        mx = jointLoadData(i,5);
        my = jointLoadData(i,6);
        mz = jointLoadData(i,7);
        load = mod.NodeLoad(node-1,PX=px,PY=py,PZ=pz,MX=mx,MY=my,MZ=mz);
        nodeLoads{i} = load;
        geo.add_node_load(load);
    end
end
geo.numbering()

% Distributed loads in elements

if size(distLoadData) >0
    frame_loads = {size(distLoadData,1)};
    for i=1:size(distLoadData,1)
        e = distLoadData(i,1);
        wx= distLoadData(i,2);
        wy = distLoadData(i,3);
        wz = distLoadData(i,4);
        load = mod.DistributedLoad(elements{e},WX=wx,WY=wy,WZ=wz);
        frame_loads{i} = load;
    end
end



O = mod.ILHEB(geo);
O.solve()
O.post_process()
O.save_results("Frame3DResults")
geo.plot("./Frame3DResults/Geometry.pdf");
geo.plot_defo(mult=MULT,n_points=int64(20),filename="./Frame3DResults/Defo.pdf")
disp("Analisys finished")


if OPEN_AT_END==true    
    mod.plt.show()
end
mod.plt.close('all')