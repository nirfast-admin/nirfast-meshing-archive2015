function [totvol tet_vol]=shell_volume(t,p,p0)
% calculates the volume of 3D shell defined in t and p.
% if p0 is provided then it will use it as its fixed point.
% it will return the volume of the shell in totvol and individual vol of
% each tetrahedron formed by p0 and triangular patches of the shell in
% tet_vol
% Note: this routine assumes that all triangular patches defined in t are
% following the same right hand rule.

if size(t,2)~= 3 || size(p,2)~=3
    error('Can only handle 3D shell surfaces');
end

% Make sure nodes in 't' are numbered from 1 to N
nodes=unique(t(:));
pp=p(nodes,:);
[tf ee]=ismember(t,nodes);

ne = size(ee,1);
% calculate centroid of shell as p0 if it's not provided.
if nargin==2
    p0 = (pp(ee(:,1),:)+pp(ee(:,2),:)+pp(ee(:,3),:));
    p0 = sum(p0)/ne;
end

tets=zeros(4,3,ne);
for i=1:3
    tets(i,:,:)=reshape(pp(ee(:,i),:)',[1 3 ne]);
end
tets(4,:,:)=reshape(repmat(p0,ne,1)',[1 3 ne]);
tet_vol = signed_tetrahedron_vol(tets);
totvol = sum(tet_vol);
