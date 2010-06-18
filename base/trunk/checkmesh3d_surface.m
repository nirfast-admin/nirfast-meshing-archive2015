function [q_radius_ratio,q_area_ratio,area,zeroflag,...
    edges,edgecnt,edgeIndexing,edgeflag]=checkmesh3d_surface(e,p,type,nodemap)
% Checks the integrity of a surface mesh by running following diagnosis:
% 'type' = 
% 0- all
% 1- area of each patch
% 2- face connectivity
% 3- edge connectivity
% 4- quality of each surface patch
e=e(:,1:3);
if nargin~=4
    nodenumbers=(1:size(p,1))';
else
    nodenumbers=nodemap;
end
% [tf idx]=ismember([nodenumbers(e(:,1)) nodenumbers(e(:,2)) nodenumbers(e(:,3))],nodenumbers);
[tf ee]=ismember(e,nodenumbers);
if sum(sum(~tf))~=0
    error('checkmesh3d_surface:BadSurfaceMesh',...
        'The provided mesh uses node numbers that are not part of the node list!');
end
[tf idx]=ismember(nodenumbers,e);
s=sum(~tf);
if s~=0
    Message={'The provided 3D surface mesh has extra nodes that';...
            'are not used in the patch element list!';...
            'Please select a filename so the fixed mesh can be written to!';...
            '';...
            };
    Title='Method of Fixing ?';
    Icon='help';
    h=msgbox(Message,Title,Icon);
    uiwait(h);
    disp(' ');
    disp('The provided mesh has extra nodes that are not used in the patch element list!');
    [FileName,PathName,FilterIndex] = uiputfile({'*.node;*.ele','Mesh Files';...
          '*.*','All Files' },'Save fixed mesh');
    FileName=remove_extension(FileName);
    FileName=sprintf('%s%s',PathName,FileName);
    pp=p(tf,:);
    writenodes_tetgen([FileName '.node'],pp,nodenumbers(tf,:));
    writeelms_tetgen([FileName '.ele'],e);
    disp(['A new file without extra nodes has been written to: ' FileName]);
    disp('Note that the new files are NOT renumbered.');
    disp(' ');
    error('checkmesh3d_surface:BadSurfaceMesh',...
        'The provided mesh has extra nodes that are not used in the patch element list');
end

% s=input('Do you want a full list of edges and their adjacent patches? (y/n)','s');
s='y';

[edgeflag,edges,edgecnt,edgeIndexing]=checkedges(e,p,nodenumbers);

% renumber based on nodemap
% [tf ee]=ismember(e,nodenumbers);

ee=uint32(ee);
[area,zeroflag]=checkarea(ee,p);
q_radius_ratio=quality_triangle_radius(p,ee);
q_area_ratio=quality_triangle_area(p,ee,area);

if sum(zeroflag)>0;
    disp(' ');
    disp('At least one of the patches has a very small area!');
    disp(' ');
end
if s~='y'
    edges=[];edgecnt=[];edgeIndexing=[];
end



function [area,zeroflag] = checkarea(e,p)
global tiny
zeroflag=false;
if isempty(tiny)
    tiny=1e-6;
end
ne=size(e,1);
area=zeros(ne,1);
area = triangle_area_3d(p(e(:,1),:),p(e(:,2),:),p(e(:,3),:));
zeroflag = IsEqual(area,0,tiny);
if sum(zeroflag)~=0
    zeroflag=true;
else
    zeroflag=false;
end



function [retflag,edges,edgecnt,indexing]=checkedges(e,p,nodemap)
% retflag=0 no problem
% retflag=1 3d mesh might be multiple material or invalid
% retfalg=2 3d mesh is invalid (or at least open)

retflag=0;
edges = [e(:,[1 2]);e(:,[1 3]);e(:,[2 3])];
edges = unique(sortrows(sort(edges,2)),'rows');
nedges=size(edges,1);
edgecnt=zeros(nedges,3);
ne=size(e,1);
temp=zeros(3,2);
i=1;
flag=0;
indexing=zeros(nedges,2);
while i<=nedges
    c=edges(i,1);
    indexing(c,1)=i;
    if i==nedges
        indexing(c,2)=nedges;
        break;
    end
    for j=i+1:nedges
        if edges(j,1)~=c
            indexing(c,2)=j-1;
            flag=1;
            break;
        end
    end
    if flag==0
        indexing(c,2)=j;
        break;
    end
    i=j;
    flag=0;
end

handle=waitbar(0,'Please wait.... Checking the edges!');
for i=1:size(e,1)
    temp(1,:)=sort([e(i,1) e(i,2)]);
    temp(2,:)=sort([e(i,1) e(i,3)]);
    temp(3,:)=sort([e(i,2) e(i,3)]);
    temp=sortrows(temp);
    for j=1:3
        starti=indexing(temp(j,1),1);
        endi  =indexing(temp(j,1),2);
        [tf idx]=ismember(temp(j,:),edges(starti:endi,:),'rows');
        if ~tf,error('Impossible!!!'),end
        ij=starti+idx-1;
        if edgecnt(ij,1)==0
            edgecnt(ij,1)=i;
        elseif edgecnt(ij,2)==0
            edgecnt(ij,2)=i;
        elseif edgecnt(ij,3)==0
            retflag=1;
            edgecnt(ij,3)=i;
        end
    end
    waitbar(i/ne,handle);
end
close(handle);
b = edgecnt(:,1)>0;
c = edgecnt(:,2)>0;
b=~b; c=~c;
sumb=sum(b); sumc=sum(c);
os=computer;
if ~isempty(strfind(os,'PCWIN')) % Windows
    newlinech ='pc';
elseif ~isempty(strfind(os,'MAC')) ||  ~isempty(strfind(os,'GLNX86')) % Mac OS or Linux
    newlinech ='unix';
end
if sumc~=0 || sumb~=0
    disp(' ');
    disp('Connectivity issue in surface mesh.');
    disp('It is possible that the mesh is not a closed surface!');
    disp(' ');
    disp('List of problematic edges were written to SurfaceMesh-InvalidEdges.txt');
    dlmwrite('SurfaceMesh-InvalidEdges.txt',edges(c,:),'delimiter',' ','newline',newlinech);
    disp(' ');
    retflag=2;
end
d=edgecnt(:,3)>0;
if sum(d)>0
    disp(' ');
    disp('*** Is this mesh a multiple material surface ? ***');
    disp('List of edges being shared by more than two triangular patches');
    disp('were written to SurfaceMesh-SharedEdges-Problem.txt');
    dlmwrite('SurfaceMesh-SharedEdges-Problem.txt',edges(d,:),'delimiter',' ','newline',newlinech);
    disp(' ');
    retflag=3;
end



