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
    v1=p(t(i,2),:)-p(t(i,1),:);
    v2=p(t(i,3),:)-p(t(i,2),:);
    normal = cross(v1,v2);
    normal = -normal / norm(normal);
    
    offsetp = mean(p(t(i,:),:)) + offset * normal;
    if PointInPolyhedron_mex(offsetp,double(t),p,tiny*100) == 1
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