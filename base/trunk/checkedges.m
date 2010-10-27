function [retflag, edges] = checkedges(e)
% retflag=0 no problem
% retflag=1 3d mesh might be multiple material or invalid
% retfalg=2 3d mesh is invalid (or at least open)
% retflag=3 both open surface and multi-material

retflag=0;
edges = [e(:,[1 2]);e(:,[1 3]);e(:,[2 3])];
[edges ix jx] = unique(sort(edges,2),'rows');
vec=histc(jx,1:max(jx));

b = vec==1; 
c = vec>2;
sumb=sum(b); sumc=sum(c);


os=computer;
if ~isempty(strfind(os,'PCWIN')) % Windows
    newlinech ='pc';
elseif ~isempty(strfind(os,'MAC')) ||  ~isempty(strfind(os,'GLNX86')) % Mac OS or Linux
    newlinech ='unix';
end


if sumb~=0 % Some edges are used only once.
    dlmwrite('SurfaceMesh-InvalidEdges.txt',edges(b,:),'delimiter',' ','newline',newlinech);
    retflag=2;
end
if sumc~=0 % Some esges are used more than twice
    if sumb~=0
        dlmwrite('SurfaceMesh-InvalidEdges.txt',edges(c,:),'delimiter',' ','-append','newline',newlinech);
    else
        dlmwrite('SurfaceMesh-InvalidEdges.txt',edges(c,:),'delimiter',' ','newline',newlinech);
    end
    if retflag==2
        retflag=3;
    else
        retflag=1;
    end
end

