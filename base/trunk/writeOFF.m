function writeOFF(fn,t,p)
% Writes surface mesh defined in 't' and 'p' to a filed called fn
% Note: At the momemnt it assumes mesh is a triangulated surface mesh only
% http://shape.cs.princeton.edu/benchmark/documentation/off_format.html
%
% Written by:
%           Hamid Ghadyani Oct 2010

fid = OpenFile(fn,'wt');

fprintf(fid,'OFF\n%d %d 0\n',size(p,1),size(t,1));

fprintf(fid,'%.16f %.16f %.16f\n',(p(:,1:3))');

fprintf(fid,'3 %d %d %d\n',(t(:,1:3)-1)');


fclose(fid);
