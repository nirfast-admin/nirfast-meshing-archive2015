function [vol,vol_ratio,zeroflag]=checkmesh3d_solid(e,p,type,nodemap)
% If type=0 then it will only check the volume of each tet, otherwise it
% will also check the quality and output some general ifno
global TetrahedronFailQuality
if isempty(TetrahedronFailQuality)
    TetrahedronFailQuality=0.03;
end
ntet=size(e,1); np=size(p,1);
os=computer;
if ~isempty(strfind(os,'PCWIN')) % Windows
    newlinech ='pc';
elseif ~isempty(strfind(os,'MAC')) ||  ~isempty(strfind(os,'GLNX86')) % Mac OS or Linux
    newlinech ='unix';
end


if nargin~=4
    nodenumbers=(1:size(p,1))';
else
    nodenumbers=nodemap;
end

% check to see if all nodes belong to at least one tet.
tf =ismember(nodenumbers,e);
s=sum(~tf);
if s~=0
    fprintf(...
        'checkmesh3d_solid:\n The provided mesh has extra nodes that are not used in the element list\n');
    warning('Meshing:check', ' Not all nodes are used in the element connectivity list');
end

nodes = unique([e(:,1);e(:,2);e(:,3);e(:,4)]);
if nargin==3
    nodemap=nodes;
end
p=p(nodes,:);
[tf1 ren_tet]=ismember(e(:,1:4),nodemap);
e=ren_tet;
sumtf=sum(tf1,2);
bf=sumtf~=4;
if sum(bf)~=0
    tempe=1:ntet;
    disp('checkmesh3d_solid: Some of the tets are using nodes that are not defined in node list!');
    dlmwrite('tets_with_extra_nodes.txt',tempe(bf),'newline',newlinech);
    error('Some of the tets are using nodes that are not defined in node list!');
end

global tiny
if isempty(tiny)
    bbx=[min(p(:,1)) min(p(:,2)) min(p(:,3)) ...
        max(p(:,1)) max(p(:,2)) max(p(:,3))];
    for i=1:3
        temp(i)=bbx(i+3)-bbx(i);
    end
    tiny=max(temp)*1e-6;
%     tiny=1e-6;
end

% Getting the min volume and finding all tets that have volumes close to
% the minimum one
vol=signed_tetrahedron_vol(e(:,1:4),p(:,1),p(:,2),p(:,3));
fprintf('Avg Min Max volume: %f %f %f\n',mean(abs(vol)),min(abs(vol)),max(abs(vol)));

% vol_ratio=quality_vol_ratio(ren_tet,p);
vol_ratio=simpqual(p,ren_tet,2);

zeroflag=vol_ratio<=TetrahedronFailQuality;

nvoids=sum(zeroflag);

if nvoids~=0
    disp(['There are ' num2str(nvoids) ' elements with undesirable quality.']);
    disp('Check voidelements.txt');
    dlmwrite('voidelements.txt',e(zeroflag,:),'delimiter',' ','newline',newlinech);
end

if type~=0
    fprintf('Avg Min Max volume ratio quality: %f %f %f\n',...
        mean(vol_ratio), min(vol_ratio), max(vol_ratio));
end

% check faces to make sure every face is only used by 1 or 2 tetrahedrons
fprintf('Checking faces..... ')
faces=[e(:,[1 2 3]);e(:,[1 2 4]);e(:,[2 3 4]);e(:,[1 3 4])];
faces=sort(faces,2);
[foo ix jx]=unique(faces,'rows');
range=1:max(jx);
vec=histc(jx,range);
bf=vec>2;
nbadfaces=sum(bf);
jx2=range(bf);
badfaces=foo(jx2,:);
badtets=[];
if nbadfaces~=0 % Some of faces are shared by more than tetrahedron: a definite problem
    fprintf('\t\n------------ Invalid solid mesh! ------------\n')
    fprintf('\tA total %d faces of the mesh are shared by more than two tetrahedrons!\n',nbadfaces)
    fprintf('\tThose faces can be found in bad_faces_extra_shared_solid.txt\n')
    fid = OpenFile('bad_faces_extra_shared_solid.txt','wt');
    for i=1:nbadfaces
        fprintf(fid,'Face: %d %d %d\t',badfaces(i,:));
        [tf idx]=ismember(badfaces(i,:),foo,'rows');
        MyAssert(vec(idx)>2);
        tmpjx = jx;
        fprintf(fid,'Tets: ');
        for j=1:vec(idx)
            [tf2 idx2]=ismember(idx,tmpjx);
            MyAssert(tf2);
            idx3 = mod(idx2,ntet); if idx3==0, idx3=ntet; end
            badtets=[badtets idx3];
            fprintf(fid,'%d ',idx3);
            tmpjx(idx2)=0;
        end
        fprintf(fid,'\n');
    end
    fclose(fid);
    badtets=unique(badtets);
%     writenodelm_3dm('bad_faces_extra_shared_solid.3dm',e(badtets,:),p);
end
fprintf('\bDone\n');


[foo p]=boundfaces(p,e(:,1:4));
fprintf('\n\n----> Checking integrity of the surface of the solid mesh...\n')
CheckMesh3D(foo,p);
fprintf('----> Done.\n\n');

