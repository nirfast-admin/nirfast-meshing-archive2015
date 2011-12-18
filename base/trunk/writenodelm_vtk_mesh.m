function writenodelm_vtk_mesh(fn,e,p,soldata)
% Writes a tetrahedra mesh defined in e and p to 'fn'
% 'soldata' is the nodal solution for every point in p:
%   soldata{1} = Name of fields/solutions
%   soldata{2} = a numnodes by d matrix where d is number of fields
%   defined in soldata{1,:}

nodes = p;
numnodes = size(nodes,1);
elems = e;
numelems = size(elems,1);

fid = fopen(fn,'wt');

%define an VTK header for FEM mesh representation
line0 = '# vtk DataFile Version 2.0';
line1 = 'NIRFAST mesh with solutions';
line2 = 'ASCII';
line3 = 'DATASET UNSTRUCTURED_GRID';
fprintf(fid,'%s\n%s\n%s\n',line0,line1,line2,line3);

line4 = ['POINTS ', num2str(numnodes), ' double']; %node defs
fprintf(fid,'%s\n',line4);
fprintf(fid, '%.18f %.18f %.18f\n', nodes');

line5 = ['CELLS ',num2str(numelems),' ',num2str(numelems*4+numelems)]; %connectivity maps
fprintf(fid,'%s\n',line5);
fprintf(fid,'%d %d %d %d %d\n',[4*ones(numelems,1) elems(:,1:4)-1]');
line6 = ['CELL_TYPES ', num2str(numelems)]; %specify the mesh basis 10-tetrahedral for all connectivity maps 
fprintf(fid,'%s\n',line6);
fprintf(fid,'%d\n', ones(numelems,1)*10);

if nargin>3 && ~isempty(soldata)
    fprintf(fid,'POINT_DATA %d\n',numnodes); %specify the data that follows is defined on the nodes

    for i = 1:size(soldata,2)
        fprintf(fid,'%s\n',['SCALARS ', soldata{1}{i}, ' float 1']);
        fprintf(fid,'%s\n','LOOKUP_TABLE default');
        fprintf(fid,'%f\n', soldata{2}(:,i));
    end;
end

if size(elems,2)>4
    fprintf(fid,'CELL_DATA %d\n', numelems);
    fprintf(fid,'SCALARS materials float 1\n');
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n', elems(:,5));

end
fclose(fid);