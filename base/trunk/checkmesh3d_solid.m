function [vol,vol_ratio,zeroflag]=checkmesh3d_solid(e,p,type,nodemap,eb,pb)
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
nodes = unique([e(:,1);e(:,2);e(:,3);e(:,4)]);
if nargin==3
    nodemap=nodes;
end
[tf1 ren_tet]=ismember(e(:,1:4),nodemap);
sumtf=sum(tf1,2);
bf=sumtf~=4;
if sum(bf)~=0
    tempe=1:ntet;
    disp('checkmesh3d_solid: Some of the tets are using nodes that are not defined in node list!');
    dlmwrite('tets_with_extra_nodes.txt',tempe(bf),'newline',newlinech);
    error('Some of the tets are using nodes that are not defined in node list!');
end
q=[];
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

tetra=zeros(4,3,ntet);
for i=1:ntet
    tetra(:,:,i)=[p(ren_tet(i,1),:); p(ren_tet(i,2),:); p(ren_tet(i,3),:); p(ren_tet(i,4),:)];
end
% Getting the min volume and finding all tets that have volumes close to
% the minimum one
vol=signed_tetrahedron_vol(tetra);
disp(' ');
disp(['Avg., Min, Max volume: ' num2str(mean(abs(vol))) ', ' ...
    num2str(min(abs(vol))) ', ' num2str(max(abs(vol)))]);
% vol_ratio=quality_vol_ratio(ren_tet,p);
vol_ratio=simpqual(p,ren_tet,2);

zeroflag=vol_ratio<=TetrahedronFailQuality;
% if type~=0
%     q=simpqual(p,ren_tet);
% end
nvoids=sum(zeroflag);

if nvoids~=0
    disp(['There are ' num2str(nvoids) ' elements with undesirable quality.']);
    disp('Check voidelements.txt');
    dlmwrite('voidelements.txt',e(zeroflag,:),'delimiter',' ','newline',newlinech);
end

if type~=0
    disp(['Avg., Min, Max volume ratio quality: ' num2str(mean(vol_ratio)) ', ' num2str(min(vol_ratio)) ', ' num2str(max(vol_ratio))]);
end

% check to see if all nodes belong to at least one tet and also that no tet
% is using a node which is not in node list
[tf idx]=ismember(nodemap,e);
if sum(~tf)~=0
    disp(' ==== checkmesh3d_solid: There are some nodes that have not been used in any tet. ===== ');
    disp('Check unused_nodes_in_tet.txt file!');
    badnodes=nodemap(~tf,:);
    dlmwrite('unused_nodes_in_tet.txt', badnodes, 'delimiter', '\t','newline',newlinech);
end

% check faces to make sure every face is only used by 1 or 2 tetrahedrons
disp(sprintf('Checking faces..... '))
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
    disp(sprintf('\t\n@@@@@@@@@@@@@@@@@ Invalid solid mesh! @@@@@@@@@@@@@@@@@\n'))
    disp(sprintf('\tA total %d faces of the mesh are shared by more than two tetrahedrons!\n',nbadfaces))
    disp(sprintf('\tThose faces can be found in bad_faces_extra_shared_solid.txt\n'))
    [fid st]=OpenFile('bad_faces_extra_shared_solid.txt','wt');
    if st==1, error('Error in writing the file bad_faces_extra_shared_solid.txt'); end
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
disp(sprintf('\bDone\n'))


% if nargin>3 && ~isempty(eb) && ~isempty(pb)
%     CalcEnclosedVolBySurfaceMesh(eb,pb);
% end


