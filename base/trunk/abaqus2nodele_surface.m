function [elem,node] = abaqus2nodele_surface(files)
% abaqus2nodele_surface(files)
% Reads list of files from cell array 'files' or from GUI and converts them
% to .node/.ele format (tetgen) format

% Get multiple abaqus inp files as input
if nargin==0
    [files, pname, filterindex] = uigetfile({'*.inp','Select Surface Files (*.inp)'},...
                         'Please select exported surfaces from Mimics.)',...
                         'MultiSelect','on');
else
    pname = '';
end
% change the type to cell if only one file is selected
if length(files)==1 || ~iscell(files)
    files = {files};
end

% read them one by one and convert them to .ele/.node format
for i=1:length(files)
    % remove extention of file name
    [fname ext] = remove_extension(files{i});
    if isempty(strfind(ext,'.inp'))
        error('Input files have to have .inp extensions!')
    end
    foo = sprintf('%s%s%s',pname,fname,ext);
    % read the actual file
    [elem,node] = read_abaqus_inp(foo);
    % Fix the orientation of the triangles
    [elem] = FixPatchOrientation(node,elem,[],1);
    % Output in 'tetgen' format
    writenodelm_nod_elm(fname,elem,node,[],[],1);
end