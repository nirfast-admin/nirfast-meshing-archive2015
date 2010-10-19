function writenodelm_surface_medit(fn,faces,p,nodemap)
% writenodelm_surface_medit(fn,faces,p,nodemap)
% Writes a shell surface defined in 'faces' and 'p' to a file called 'fn'
% in medit format (.mesh) so it can be viewed using medit.
% Note that this routine renumbers the node numbers from 1 to N and accordingly
% adjust the node numbers in face list.
fprintf('%s','Writing data to file... ')

ref=0;
face_att=1;
np=size(p,1); dim=size(p,2);
nf=size(faces,1);
if size(faces,2)~=3
    disp('The provided face list is not a triangular facet list!');
    error('writenodelm_surface_medit() can only write shell surfaces to medit format');
end
if dim~=3
    error('Input points are not 3D');
end
if isempty(fn), error('Filename should be specified'), end
os=computer;
if ~isempty(strfind(os,'PCWIN')) % Windows
    newlinech ='pc';
elseif ~isempty(strfind(os,'MAC')) ||  ~isempty(strfind(os,'GLNX86')) % Mac OS or Linux
    newlinech ='unix';
end
fn = add_extension(fn,'.mesh');
    
fid=fopen(fn,'wt');
fprintf(fid,'MeshVersionFormatted 1\r');
fprintf(fid,'Dimension\n3\nVertices\n%u\r',np);

if nargin==4
    if size(nodemap,1)==1
        nn=nodemap';
    elseif size(nodemap,2)==1
        nn=nodemap;
    elseif size(nodemap,2)~=1
        error('The nodemap list provided should be an n by 1 vector');
    end
    [mytflag newfaces]=ismember(faces,nodemap);
    sum1=sum(mytflag,2);
    sumbflag=sum1==3;
    if sum(sumbflag)~=nf % not all vertices of 'faces' are given in nodemap list!
        disp(' ');
        disp('writenodelm_surface_medit: Not all vertices of the input face could be found in given nodemap list!')
        error('writenodelm_surface_medit: faces and nodemap lists are not compatible!')
    end
%     newfaces=nodemap(faces);
%     newfaces=faces;
else
    nn=(1:np)';
    newfaces=faces;
end
% dlmwrite(fn,[nn p],'precision',8,'delimiter',' ','-append','newline',newlinech);
fprintf(fid,'%.10f %.10f %.10f %d\n',[p(:,1:3) ones(size(p,1),1)*ref]');
fprintf(fid,'Triangles\n%u\n',nf);

fprintf(fid,'%d %d %d %d\n',[newfaces ones(nf,1)*face_att]');
fprintf(fid,'End\n');
fclose(fid);


[mypath myfn myext]=fileparts(fn);
fprintf('\b%s\n','Done writing surface mesh to:')
if ~isempty(mypath)
    fprintf('\t%s\n',['Path: ' mypath])
end
fprintf('\t%s\n',['Filename: ' myfn myext])
