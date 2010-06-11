function mesh = checkerboard3d_mm(fnprefix, type)
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
    [fname, pname, filterindex] = uigetfile({'*.inp','Select one of exported files from Mimics (*.inp)'},...
                             'Please select a 3D surface mesh file',...
                             'MultiSelect','off');
    [fnprefix myext] = remove_extension(fname);
    fnprefix = [pname fnprefix];
else
    [fnprefix myext] = remove_extension(fnprefix);
end
if isempty(myext) || (~strcmpi(myext,'.inp') && ~strcmpi(myext,'.ele'))
    errordlg('Surface filenaem should have either .inp or .ele as its extension','Meshing Error');
    error('Surface filenaem should have either .inp or .ele as its extension');
end
if nargin~=2
    type='generic';
end

fprintf('\n\n--> Beginning mesh generation process...\n\n');
%% Read in the mesh
if strcmpi(myext,'.inp') % Each INP file represents a surface in mesh with disjoint sub regions
    no_regions = length(dir([fnprefix '*.inp']));
    if no_regions==0
        errordlg(['Cannot find file .inp files whose prefix is ' fnprefix],'Mesh Error');
        error(['Cannot find file .inp files whose prefix is ' fnprefix]);
    end

    interior_nodes=zeros(no_regions,3);

    fprintf('\n\tConverting inp files and re-orienting\n');
    tnn = 0;
    telem = [];
    tnode = [];
    flag= true; fcounter = 1;
    while flag
        fn = [fnprefix num2str(fcounter) '.inp'];
        [fid msg]=fopen(fn,'rt');
        if ~isempty(msg), flag=false; continue; end
        fclose(fid);
        fprintf('  Reading mesh file: %s\n', fn);
        [elem,node] = read_abaqus_inp(fn);
        % Fix the orientation of the triangles
        [elem] = FixPatchOrientation(node,elem,[],1);
        tmp = GetOneInteriorNode(elem,node);
        if isempty(tmp)
            errordlg(['Could not find an interior point in surface: ' fn],'Mesh Error');
            error(['Could not find an interior point in surface: ' fn]);
        end
        interior_nodes(fcounter,:) = tmp;
        % Save the exteriro surface
        if fcounter == 1;
            extelem = elem; extnode = node;
        end
        cnn = size(node,1);
        elem = elem + tnn;
        tnn = tnn + cnn;
        telem = [telem;elem];
        tnode = [tnode;node];
        fcounter = fcounter + 1;
    end
    if fcounter==1 % couldn't read the first file
        errordlg({[msg ':'], fn},'Meshing Error');
        error(msg,'Meshing Error');
    end
    regions=cell(size(interior_nodes,1),2);
    for i=1:size(interior_nodes,1)
        regions{i,1} = interior_nodes(i,:);
        regions{i,2} = i;
    end
elseif strcmpi(myext,'.ele')
%     One single .node/.ele file is used to represent the mesh. This is for
%     meshes whose interior sub-surfaces share some triangles. This format
%     is the default output of the MMC.m marching cube algorithm.
%     This routine assumes that the smallest region id (excluding 0) is the
%     id for the most exterior region which encloses all other sub regions

    [tele tnode] = read_nod_elm(fnprefix,1);
    if size(tele,2)~=3 && size(tele,2)<5
        error(['Input surface mesh does not seem to be a triangular mesh: ' fnprefix]);
    end
    if size(tele,2)==3
        nregions=1;
        regions=1;
        tele = [tele repmat([1 0], size(tele,1), 1)];
    else
        regions = tele(:,4:5);
        regions = unique(regions(:));
        nregions = length(regions);
    end
    [tf idx]=ismember(0,regions);
    if tf, regions(idx)=[]; nregions=nregions-1; end % Remove 0 (which represents outside space)
    
    interior_nodes=zeros(nregions,3);
    telem=[];
    for i=1:nregions
        bf=tele(:,4)==regions(i) | tele(:,5)==regions(i);
        reg_ele=tele(bf,1:3);
        tmp = GetOneInteriorNode(reg_ele,tnode);
        if isempty(tmp)
            errordlg(['Could not find an interior point in surface: ' fn],'Mesh Error');
            error(['Could not find an interior point in surface: ' fn]);
        end
        interior_nodes(i,:) = tmp;
        if i == 1;
            extelem = reg_ele;
        end
        telem=[telem;reg_ele];
    end
    regions=cell(size(interior_nodes,1),2);
    for i=1:size(interior_nodes,1)
        regions{i,1} = interior_nodes(i,:);
        regions{i,2} = i;
    end
end

[foo ix] = unique(sort(telem(:,1:3),2),'rows');
clear foo
telem=telem(ix,:);
myargs.silentflag=1;
myargs.bdyfn = [fnprefix];
myargs.regions = regions;
% Remove the following line to create a mesh based on average size of the
% input surfaces
% myargs.edgesize = 4.1;

myargs.examineinpmesh=0; % Do not run inspection checks on input surface
myargs.extelem=extelem;

% writenodelm_nod_elm(myargs.bdyfn,extelem,extnode)
[mesh.elements, mesh.nodes] = checkerboard3d(telem,tnode,myargs);
delete('input4delaunay.*','junk.txt');

%% Write NIRFAST-format mesh files
fn = [fnprefix '-nirfast'];
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

fprintf('\n\n--> Finished mesh generation.\n');












