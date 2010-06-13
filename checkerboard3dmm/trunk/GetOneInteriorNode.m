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

for i=1:ne
    v1=p(t(i,2),1:3)-p(t(i,1),1:3);
    v2=p(t(i,3),1:3)-p(t(i,2),1:3);
    normal = cross(v1,v2);
    normal = -normal / norm(normal);
    
    if any(isnan(normal)), continue; end
    offsetp = mean(p(t(i,:),1:3)) + offset * normal;
%     st = PointInPolyhedron_mex(offsetp,double(t),p,tiny*100);
    st = involume_mex(offsetp, double(t), p, 200, min(p(:,1)), max(p(:,1)), tiny);
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
