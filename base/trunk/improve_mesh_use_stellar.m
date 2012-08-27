function [e p st] = improve_mesh_use_stellar(e, p, opt_params)
if nargin < 3 || isempty(opt_params)
    qualmeasure = 0;
    facetsmooth = 0;
    usequadrics = 1;
    opt_params = [];
end
eorig = e;

mycf = pwd;
foodir = fileparts(which('CheckMesh3D.m'));
stellar_config_fn = fullfile(foodir,'.stellar_config');

h = waitbar(0,'Initializing optimization.');

% Make sure element orientations are OK
e = check_element_orientation_3d(e,p);

% Write temp element files for stellar
tmpfolder = tempdir;
cd(tmpfolder)

fnprefix = fullfile(tmpfolder,'stellar_input');
writenodelm_nod_elm(fnprefix, e, p)

if ~isempty(opt_params) && isfield(opt_params,'qualmeasure')
    qualmeasure = opt_params.qualmeasure;
end
if ~isempty(opt_params) && isfield(opt_params,'facetsmooth')
    facetsmooth = opt_params.facetsmooth;
end
if ~isempty(opt_params) && isfield(opt_params,'usequadrics')
    usequadrics = opt_params.usequadrics;
end
% Create a config file for Stellar based on the template
fid = fopen([stellar_config_fn '_temp'],'rt');
foos = fscanf(fid,'%c',Inf);
fclose(fid);
foos = strrep(foos,'__QUALMEASURE__', num2str(qualmeasure));
foos = strrep(foos,'__FACETSMOOTH__', num2str(facetsmooth));
foos = strrep(foos,'__USEQUADRICS__', num2str(usequadrics));
fid = fopen(stellar_config_fn,'wt');
fprintf(fid,'%c',foos);
fclose(fid);

execname = GetSystemCommand('improve_mesh_stellar');

waitbar(0.1,h,'Running optimization.');

command = ['"' execname '" -s ' '"' stellar_config_fn ...
    '" -j -C -V "' fnprefix '"'];

fprintf('\n ** Running Optimizier, please wait');
tic
[st result] = system(command);
t2=toc;
fprintf(' **\n');

if st~=0
    warning('nirfast:improvemesh', ' --> Could not improve the mesh:');
    fprintf('--------------------------------\n%s',result);
    fprintf('--------------------------------\n');
    e = eorig;
    close(h);
else
    [ms me foo mstr] = regexp(result,'worstqual:\s*[\d.+-]+');
    bqual = mstr{1};
    aqual = mstr{2};
    fprintf('\n -- Quality Optimziation --\n\tBefore: %s\n\tAfter : %s\n%s'...
        ,bqual,aqual);
    fprintf('\tTime: %.2f secs\n',t2);
    waitbar(0.96,h,'Reading optimized mesh.');
    [e p] = read_nod_elm([fnprefix '.1.'],1);
    % Remove unused nodes
    ee = e(:,1:4);
    nodes = unique(ee(:));
    p = p(nodes,:);
    [tf ee] = ismember(ee,nodes);
    e(:,1:4) = ee;
    close(h);
end

delete([fnprefix '.*'])
cd(mycf)


