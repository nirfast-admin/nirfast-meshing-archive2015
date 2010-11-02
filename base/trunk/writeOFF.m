function writeOFF(fn,t,p)
% Writes surface mesh defined in 't' and 'p' to a filed called fn
% Note: At the momemnt it assumes mesh is a triangulated surface mesh only

fid = OpenFile(fn,'wt');

fprintf(fid,'%d %d 0\n',size(t,1),size(p,1));

fprintf(fid,'%d %d %d\n',(t(:,1:3)-1)');

fprintf(fid,'%.16f %.16f %.16f\n',(p(:,1:3))');

fclose(fid);
