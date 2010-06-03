function [elem, points] = read_abaqus_inp_3D(mesh_fn)
% Reads abaqus mesh files exported by Mimics V.13
% 
% Written by Hamid Ghadyani, Apr 2010
% The old version which is very slow is now called read_abaqus_inp_line
% 

if(nargin~=1)
    [fname,pname]=uigetfile('*.inp','Please pick the inp file');
else
    pname = [];
	fname=mesh_fn;
end
[fid st]=OpenFile([pname fname],'r');
if st~=0, error(' '); end

% read header junk
for i=1:4, fgetl(fid); end
% read node coordinates
data=textscan(fid,' %u64, %.54f, %.54f, %.54f%*[^\n]');
points = [data{2} data{3} data{4}];

% read header
s = fgetl(fid);
if ~isempty(strfind(s,'S3R'))
    pattern = '%u64, %u64, %u64, %u64%*[^\n]'; % surface mesh
elseif ~isempty(strfind(s,'C3D4'))
    pattern = ' %u64, %u64, %u64, %u64, %u64%*[^\n]'; % solid/tetrahedral mesh
end

% read element list
data=textscan(fid,pattern);
if ~isempty(strfind(s,'S3R'))
    elem  = double([data{2} data{3} data{4}]); % surface mesh
elseif ~isempty(strfind(s,'C3D4'))
    elem  = double([data{2} data{3} data{4} data{5}]); % solid/tetrahedral mesh
end


fclose(fid);
