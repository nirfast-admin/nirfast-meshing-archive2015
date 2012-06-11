function [e p] = RunCGALMeshGenerator(mask,param)
% Runs CGAL mesh generator on 3D matrix 'mask'
% For info on 'param' and 'info' refer to 'image2mesh_cgal.m'
% 
% Written by:
%   Hamid Ghadyani May 2011

if ~isfield(param,'tmppath')
    tmppath=tempdir;
else
    tmppath=param.tmppath;
end

if ~isfield(param,'delmedit')
    delmedit = 1;
elseif param.delmedit==0
    delmedit = 0;
else
    delmedit = 1;
end

tmpmeshfn = [tmppath filesep '._out.mesh'];
tmpinrfn  = [tmppath filesep '._cgal_mesher.inr'];
cgalparam_fn = [tmppath filesep '._criteria.txt'];

savefn = add_extension(tmpinrfn,'.inr');
saveinr(mask,savefn,param);

% Set up the necessary parameters for meshing
facet_angle = 25; facet_size = 3; facet_distance = 2;
cell_radius_edge = 3; cell_size = 3; % general tet size of all regions
special_subdomain_label = 0; % label of region to be refined
special_size = 0; % tet size of the region 'special_subdomain_label'
if isfield(param,'facet_angle'), facet_angle = param.facet_angle; end
if isfield(param,'facet_size'),  facet_size  = param.facet_size; end
if isfield(param,'facet_distance'), facet_distance = param.facet_distance; end
if isfield(param,'cell_radius_edge'), cell_radius_edge = param.cell_radius_edge; end
if isfield(param,'cell_size'), cell_size = param.cell_size; end
if isfield(param,'special_subdomain_label'), special_subdomain_label = param.special_subdomain_label; end
if isfield(param,'special_subdomain_size'), special_size = param.special_subdomain_size; end

% Write up the parameter files
fid = fopen(cgalparam_fn,'wt');
fprintf(fid,'%f\n',facet_angle);
fprintf(fid,'%f\n',facet_size);
fprintf(fid,'%f\n',facet_distance);
fprintf(fid,'%f\n',cell_radius_edge);
fprintf(fid,'%f\n',cell_size);
fprintf(fid,'%d\n',special_subdomain_label);
fprintf(fid,'%f\n',special_size);
fclose(fid);

% Run the executable
syscommand = GetSystemCommand('image2mesh_cgal');
if ~ispc
    eval(['! chmod u+x "' syscommand '"']);
end
makemeshcommand = ['! "' syscommand '" "' savefn '" "' cgalparam_fn '" "' tmpmeshfn '"'];
eval(makemeshcommand);

% Read the resulting mesh
[e p] = readMEDIT(tmpmeshfn);
if isfield(param,'Offset')
    p = p + repmat(param.Offset,size(p,1),1);
end

% Remove possible extra nodes that might be left out in 'p' list
% CGAL tends to do this.

nodes = unique([e(:,1);e(:,2);e(:,3);e(:,4)]);
p = p(nodes,:);
[tf ee] = ismember(e(:,1:4),nodes);
if size(e,2) > 4
	e = [ee e(:,5:end)];
else
    e = ee;
end
warning('off','MATLAB:DELETE:FileNotFound');
delete(cgalparam_fn,tmpinrfn);
if delmedit, delete(tmpmeshfn); end
warning('on','MATLAB:DELETE:FileNotFound');

