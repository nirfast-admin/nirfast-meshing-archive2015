function mesh = checkerboard3d_mm(filename, type)
% checkerbaord3d_mm(fnprefix, type)
% Reads input surfaces which are either in .inp format (Mimics exported in Abaqus
% file format) or .ele format (tetgen format) and then generates a 3D 
% multiple material tetrahedral mesh.
% 
% The routine returns a 'mesh' structure which is based on NIRFAST method, and
% its fields are:
% .node, .elements, .bndvtx, .region, .name, .type, .dimension
% Written By:
%           Hamid R Ghadyani, May 2010

%% Check the proper input filename
if nargin==0
    [fname, pname, filterindex] = uigetfile({'*.inp;*.ele','Select one of exported files from Mimics (*.inp)'},...
                             'Please select a 3D surface mesh file',...
                             'MultiSelect','off');
    filename =[pname fname];
end
[path fnprefix num_flag myext] = GetFilenameNumbering(filename);
fnprefix=fullfile(path,fnprefix);

if isempty(myext) || (~strcmpi(myext,'.inp') && ~strcmpi(myext,'.ele'))
    errordlg('Surface filenaem should have either .inp or .ele as its extension','Meshing Error');
    error('Surface filenaem should have either .inp or .ele as its extension');
end
if nargin~=2
    type='generic';
end

fprintf('\n\n--> Beginning mesh generation process...\n\n');
%% Read in the mesh
if strcmpi(myext,'.inp') 
    % Each INP file represents a surface in mesh with disjoint sub regions
    % The filename with smallest counter in its name should always be the most
    % exterior region enclosing other inp surfaces.
    no_regions = length(dir([fnprefix '*.inp']));
    if no_regions==0
        errordlg(['Cannot find file .inp files whose prefix is ' fnprefix],'Mesh Error');
        error(['Cannot find file .inp files whose prefix is ' fnprefix]);
    end

    fprintf('\n\tConverting inp files and re-orienting\n');

    telem = [];
    tnode = [];
    fcounter = num_flag;
    if num_flag==0
        fn = [fnprefix '.inp'];
    else
        fn = [fnprefix num2str(fcounter) '.inp'];
    end
    
    newmatc = 1;
    tags={};
    while true
        fid = fopen(fn,'rt');
        if fid < 0, break; end
        fclose(fid);
        [celem,cnode] = abaqus2nodele_surface(fn);
        if fcounter == num_flag
            extelem = celem;
        end
        ind_regions = GetIndRegions(celem,cnode);
        for i=1:length(ind_regions)
            pin = GetOneInteriorNode(celem,cnode);
            if isempty(pin)
                errordlg('Could not find an interior point in surface.','Mesh Error');
                error('Could not find an interior point in surface.' );
            end
            tags{newmatc,1} = pin;
            tags{newmatc,2} = newmatc;
            
            foo_ele = celem(ind_regions{i},1:3);
            foo_nodes = unique(foo_ele(:));
            foo_coords = cnode(foo_nodes,:);
            [tf foo_ele] = ismember(foo_ele,foo_nodes);
        
            foo_ele = FixPatchOrientation(foo_coords,foo_ele,[],1);
            telem = [telem; foo_ele+size(tnode,1)];
            tnode = [tnode; foo_coords];
            newmatc = newmatc + 1;
        end
        fcounter = fcounter + 1;
        fn = [fnprefix num2str(fcounter) '.inp'];
    end
   
elseif strcmpi(myext,'.ele')
%     One single .node/.ele file is used to represent the mesh. This is for
%     meshes whose interior sub-surfaces share some triangles. This format
%     is the default output of the MMC.m marching cube algorithm.
%     This routine assumes that the smallest region id (excluding 0) is the
%     id for the most exterior region which encloses all other sub regions
    if num_flag~=0
        fnprefix = [fnprefix num2str(num_flag)];
    end
    [telem tnode] = read_nod_elm(fnprefix,1);
    output = SeparateSubVolumes(telem, tnode);
    tags = output.tags;
    extelem = output.extelem;
end



extelem = FixPatchOrientation(tnode,extelem,[],1);


[foo ix] = unique(sort(telem(:,1:3),2),'rows');
clear foo
telem=telem(ix,:);
myargs.silentflag=1;
myargs.bdyfn = fnprefix;
myargs.regions = tags;
% Remove the following line to create a mesh based on average size of the
% input surfaces
% myargs.edgesize = 4.1;

myargs.examineinpmesh=0; % Do not run inspection checks on input surface
myargs.extelem=extelem;

clear output
%% Call the main checkerboard3d routine
[mesh.elements, mesh.nodes] = checkerboard3d(telem(:,1:3),tnode,myargs);
delete('input4delaunay.*','junk.txt');

%% Write NIRFAST-format mesh files
mesh.dimension = 3;
mesh.type = type;
mesh.name = fnprefix;

% First figure out the exterior nodes
faces=[mesh.elements(:,[1,2,3]);
       mesh.elements(:,[1,2,4]);
       mesh.elements(:,[1,3,4]);
       mesh.elements(:,[2,3,4])];
faces=sort(faces,2);
[foo,ix,jx]=unique(faces,'rows');
vec=histc(jx,1:max(jx));
qx = vec==1;
bdy_faces=faces(ix(qx),:);
exterior_nodes_id = unique(bdy_faces(:));
mesh.bndvtx = zeros(size(mesh.nodes,1),1);
mesh.bndvtx(exterior_nodes_id) = 1;

% Get region info
% Since NIRFAST is a node based FEM package we need to assign
% material/region IDs to nodes rather than elements!
region_ids = unique(mesh.elements(:,5));
if min(region_ids)==0, region_ids=region_ids+1; end

mesh.region = zeros(size(mesh.nodes,1),1);
for i=1:length(region_ids)
    relem = mesh.elements(mesh.elements(:,5)==region_ids(i),:);
    rnodes = unique([relem(:,1);relem(:,2);relem(:,3);relem(:,4)]);
    mesh.region(rnodes) = region_ids(i);
end
mesh.elements=mesh.elements(:,1:4);

mesh.elements=mesh.elements(:,1:4);
fprintf('\n\n--> Finished mesh generation.\n');

