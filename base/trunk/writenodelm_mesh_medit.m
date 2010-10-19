function writenodelm_mesh_medit(fn,tet,p,nodemap)
% writenodelm_mesh_medit(fn,tet,p)
% Writes a tetrahedral mesh defined in 'tet' and 'p' to a file called 'fn'
% in medit format (.mesh) so it can be viewed using medit.
% Note that this routine renumbers the node numbers from 1 to N and accordingly
% adjust the node numbers in mesh list.

fprintf('%s','Writing data to file... ')
ref=0;
tet_mat=10;
np=size(p,1); dim=size(p,2);
ntet=size(tet,1);
if size(tet,2)<3 && size(tet,2)~=3
    disp('The provided mesh is not a tetrahedral or triangular mesh!');
    error('writenodelm_mesh_medit() can only write tetrahedral/triangular mesh to medit format');
end
if dim~=3
    error('Input points are not 3D');
end
if isempty(fn), error('Filename should be specified'), end

if size(tet,2)==5 % material list for a multiple material mesh
    mmflag = true;
else
    mmflag = false;
end

os=computer;
if ~isempty(strfind(os,'PCWIN')) % Windows
    newlinech ='pc';
elseif ~isempty(strfind(os,'MAC')) ||  ~isempty(strfind(os,'GLNX86')) % Mac OS or Linux
    newlinech ='unix';
end
if nargin==4
    if size(nodemap,1)==1
        nn=nodemap';
    elseif size(nodemap,2)==1
        nn=nodemap;
    elseif size(nodemap,2)~=1
        error('The nodemap list provided should be an n by 1 vector');
    end
    [mytflag newtets]=ismember(tet(:,1:4),nodemap);
    if mmflag
        newtets = [newtets tet(:,5)];
    end
    sum1=sum(mytflag,2);
    sumbflag=sum1==4;
    if sum(sumbflag)~=ntet % not all vertices of 'faces' are given in nodemap list!
        disp(' ');
        disp('writenodelm_surface_medit: Not all vertices of the input tet could be found in given nodemap list!')
        error('writenodelm_surface_medit: tets and nodemap lists are not compatible!')
    end
%     newfaces=nodemap(faces);
%     newfaces=faces;
else
    nn=(1:np)';
    newtets=tet;
end

fid=fopen(fn,'wt');
fprintf(fid,'MeshVersionFormatted 1\n');
fprintf(fid,'Dimension\n3\nVertices\n%u\n',np);

% fclose(fid);
% dlmwrite(fn,[nn p],'precision',10,'delimiter',' ','-append','newline',newlinech);
% for i=1:np
%     fprintf(fid,'%.12f %.12f %.12f %d\n',p(i,1),p(i,2),p(i,3),ref);
% end

fprintf(fid,'%.12f %.12f %.12f %d\n',[p(:,1:3) ones(np,1)*ref]');

% fid=fopen(fn,'at+');
if size(tet,2)>=4
    fprintf(fid,'Tetrahedra\n%u\n',ntet);
elseif size(tet,2)==3
    fprintf(fid,'Triangles\n%u\n',ntet);
end
% fclose(fid);
if mmflag
%     dlmwrite(fn,newtets,'precision','%u','delimiter',' ','-append','newline',newlinech);
    fprintf(fid,'%d %d %d %d %d\n',(newtets(:,1:5))');
else
%     dlmwrite(fn,[newtets ones(ntet,1)*tet_mat],'precision','%u','delimiter',' ','-append','newline',newlinech);
    fprintf(fid,'%d %d %d %d %d\n',[newtets(:,1:4) ones(ntet,1)*tet_mat]');
end

% fid=fopen(fn,'at+');
fprintf(fid,'End\n');
fclose(fid);

[mypath myfn myext]=fileparts(fn);
fprintf('\b\b%s\n','Done writing mesh to:')
if ~isempty(mypath)
    fprintf('\t%s\n',['Path: ' mypath])
end
fprintf('\t%s\n',['Filename: ' myfn myext])
