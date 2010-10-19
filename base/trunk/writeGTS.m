function writeGTS(fn,t,p)
% Writes the triangulated surface in 't' and 'p' into a .gts file which is
% based on GNU Triangulated Surface format:
% http://gts.sourceforge.net
% 
% Written by:
%   Hamid Ghadyani Sep 2010

% if size(t,2)~=3 || size(p,2)~=3
%     error('Can only handle triangulated surfaces.')
% end

fprintf('Writing data to file... ')

nt = size(t,1);
edges = [t(:,[1 2]);t(:,[1 3]);t(:,[2 3])];
[e ix jx] = unique(sort(edges,2),'rows');

fid = OpenFile(fn,'wt');

fprintf(fid,'%d %d %d\n',size(p,1),size(e,1),nt);
fprintf(fid,'%.18f %.18f %.18f\n',p');
fprintf(fid,'%d %d\n',e');

eidx = [(1:nt)' (1:nt)'+nt (1:nt)'+2*nt];
te = jx(eidx);
fprintf(fid,'%d %d %d\n',te');
fclose(fid);



[mypath myfn myext]=fileparts(fn);
fprintf('\b\b%s\n','Done writing surface mesh to:')
if ~isempty(mypath)
    fprintf('\t%s\n',['Path: ' mypath])
end
fprintf('\t%s\n',['Filename: ' myfn myext])

