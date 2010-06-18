function vol = signed_tetrahedron_vol(tetra)
% computes the volume of a tetrahedron in 3D.
%  Author:
%
%    Hamid Ghadyani
%
%  Parameters:
%
%    Input, real TETRA(4,3,:), the vertices of the tetrahedron.
%
%    Output, real VOLUME, the signed volume of the tetrahedron.
%  It returns a positive value if one follows the right hand rule for first
%  three nodes and his thumb points away from the forth one
% tetra=tetra';
s=size(tetra);
if size(s,2)==2
    ntet=1;
else
    ntet=s(3);
end
vol=zeros(ntet,1);
u=ones(4,ntet);
tetra(:,4,:)=u;
for i=1:ntet
%     vol(i) = 1/6 * det([tetra(2,:,i); tetra(1,:,i); tetra(3,:,i); tetra(4,:,i)]);
    vol(i) = 1/6 * det(tetra(:,:,i));
end

