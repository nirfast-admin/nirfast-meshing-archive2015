function mesh = checkerboard3d_mm(fnprefix, type)
% checkerbaord3d_mm(fnprefix)
% Reads input surfaces which are in .inp format (Mimics exported in Abaqus
% file format) and generates a 3D multiple material tetrahedral mesh.
% It also writes 3 files that are compatible with NIRFAST mesh format:
% .node, .elem and .region
% 
% The other mesh format that is written by this routine is 'tetgen' style
% mesh formats: .node and .ele files. Their names are in following format:
%  fnprefix-cb3d.node and fnprefix-cb3d.ele
% This routines a 'mesh' structure which is based on NIRFAST method, and
% its fields are:
% .node, .elements, .bndvtx, .region, .name, .type, .dimension

if nargin==0
    [fname, pname, filterindex] = uigetfile({'*.inp','Select one of exported files from Mimics (*.inp)'},...
                             'Please select a 3D surface mesh file',...
                             'MultiSelect','off');
    idx = regexpi(fname,'[0-9]*\.');
    fnprefix = fname(1:idx(end)-1);
    fnprefix = [pname fnprefix];
end
if nargin~=2
    type='generic';
end

fprintf('\n\n--> Beginning mesh generation process...\n\n');
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
myargs.silentflag=1;
myargs.bdyfn = [fnprefix];
myargs.regions = regions;
% Remove the following line to create a mesh based on average size of the
% input surfaces
% myargs.edgesize = 4.1;

myargs.examineinpmesh=0; % Do not run inspection checks on input surface

writenodelm_nod_elm(myargs.bdyfn,extelem,extnode)
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

fprintf('\n\n--> Finished writing NIRFAST-format mesh files: %s {.node, .elem, .region}\n\n',fn);












