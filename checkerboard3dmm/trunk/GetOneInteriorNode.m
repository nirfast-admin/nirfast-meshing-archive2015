function pin = GetOneInteriorNode(t,p)
% Asssuming a closed manifold surface is defined in 't' and 'p'
% this routine returns coordinates of a point within the surface.
% It assumes the surface is oriented and the normals of the triangles are
% pointing outward.

ne=size(t,1);
foundflag=false;

bbx1=max(p);
bbx2=min(p);
ld=min(bbx1-bbx2);

tiny = eps;

offset = tiny*100 * ld;

if offset < eps
    pin=[];
    return
end

% Make sure nodes in 't' are numbered from 1 to N
nodes=unique(t(:));
pp=p(nodes,:);
[tf ee]=ismember(t,nodes);

% Check the subvolume's integrity
input_args.verbose=0;
input_args.type=1;
[~,~,~,myst] = CheckMesh3D(ee,pp,[],input_args);

if isfield(myst,'b') && myst.b~=0 && myst.b~=4
%     writenodelm_surface_medit('foo.mesh',ee,pp)
%     system('/usr/local/bin/medit foo.mesh')
%     disp(' ')
    cprintf([1 0.5 0.5],'Warning (GetOneInteriorNode.m):\n   The given surface is not closed or single material!\n');
end

for i=1:ne
    v1=pp(ee(i,2),1:3)-pp(ee(i,1),1:3);
    v2=pp(ee(i,3),1:3)-pp(ee(i,2),1:3);
    normal = cross(v1,v2);
    normal = -normal / norm(normal);
    
    if any(isnan(normal)), continue; end
    offsetp = mean(pp(ee(i,:),1:3)) + offset * normal;
%     st = PointInPolyhedron_mex(offsetp,double(ee),pp,tiny*100);
    clear mex
    st = involume_mex(offsetp, double(ee), pp, 200, min(pp(:,1)), max(pp(:,1)), tiny);
    if st == 1
        foundflag = true;
        break
    end
end

if ~foundflag
    pin=[];
    return
else
    pin=offsetp;
end
