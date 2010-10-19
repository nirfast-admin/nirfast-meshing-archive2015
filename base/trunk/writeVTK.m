function writeVTK(fn,e,p)

if size(e,2)==3
    writevtk_polydata(fn,e,p)
else
    error('Can not handle this type of mesh')
end

function writevtk_polydata(fn,e,p)

fid = OpenFile(fn,'wt');

fprintf(fid,'# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET POLYDATA\n');

fprintf(fid,'POINTS %d double\n',size(p,1));

fprintf(fid,'%f %f %f\n',(p(:,1:3))');

fprintf(fid,'POLYGONS %d %d\n',size(e,1), size(e,1)*4);
fprintf(fid,'3 %d %d %d\n',(e(:,1:3)-1)');

fclose(fid);
