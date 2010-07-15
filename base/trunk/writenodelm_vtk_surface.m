function writenodelm_vtk_surface(fn,t,p)
% writenodelm_vtk_solid(fn,tet,p)
% Writes the mesh defined in 'tet' and 'p' to a file called 'fn'
% It will renumber the nodes from 1 to N

fprintf('%s','Writing data to file... ')

fn = add_extension(fn,'.vtk');
fid=OpenFile(fn,'wt');

t=t(:,1:3);
np = size(p,1); ne = size(t,1);

[mypath myfn myext]=fileparts(fn);

% Write the header info
fprintf(fid,'# vtk DataFile Version 2.0\n');
fprintf(fid,['vtk ' myfn '\n']);
fprintf(fid,'ASCII\n');
fprintf(fid,'DATASET POLYDATA\n');

fprintf(fid,'POINTS %d double\n',np);
fprintf(fid,'%.18f %.18f %.18f\n',p');
fprintf(fid,'\n');

fprintf(fid,'POLYGONS %d %d\n',ne,4*ne);
fprintf(fid,'3 %d %d %d\n',(t-1)');
fprintf(fid,'\n');

fclose(fid);

fprintf('\b%s\n','Done writing mesh to:')
if ~isempty(mypath)
    fprintf('\t%s\n',['Path: ' mypath])
end
fprintf('\t%s\n',['Filename: ' myfn myext])
