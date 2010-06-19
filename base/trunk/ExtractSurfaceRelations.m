function mesh = ExtractSurfaceRelations(fnprefix, no_regions)
% Reads disjoint surface meshes from each file (total of 'no_region' files)
% and establishes their relational hierarchy.
% The output is a structure compatible with nirfast format. The relation
% ship information is stored in .region field, which provides the
% information about which two regions each triagnle is separating.
% This routine assumes the surface filenames are in following format:
% fnprefix1.node fnprefix1.ele
% fnprefix2.node fnprefix2.ele
% fnprefix3.node fnprefix3.ele
% ...
% fnprefixno_regions.node fnprefixno_regions.ele

fprintf('\tSurface detection.\n\t\t');
clear mex
surface_relations_mex(fnprefix,no_regions);
relations = textread('surface_relations.txt');

region_id=zeros(no_regions,1);
idcounter=1;
for i=1:size(relations,1)
    if i==1
        fn = [fnprefix num2str(relations(i,1))];
        [mesh.elements,mesh.nodes] = read_nod_elm(fn,1);
        region1 = repmat(1,[size(mesh.elements,1) 1]);
        region2 = repmat(0,[size(mesh.elements,1) 1]);
        no_ext_nodes = size(mesh.nodes,1);
        NoBdyNodes(i) = no_ext_nodes;
        NoBdyElems(i) = size(mesh.elements,1);
    end
    if region_id(relations(i,1)) == 0
        region_id(relations(i,1)) = idcounter;
        idcounter = idcounter + 1;
    end
    momid = region_id(relations(i,1));
    for j=2:size(relations,2)
        if relations(i,j)==0, continue; end
            fn = [fnprefix num2str(relations(i,j))];
            kidid = idcounter;
            region_id(relations(i,j)) = idcounter;
            
            [fooele,foonode] = read_nod_elm(fn,1);
            mesh.elements=[mesh.elements;fooele+size(mesh.nodes,1)]; 
            mesh.nodes=[mesh.nodes; foonode]; 
            region1 = cat(1,region1,repmat(momid,[size(fooele,1) 1]));
            region2 = cat(1,region2,repmat(kidid,[size(fooele,1) 1]));
            
            NoBdyNodes(idcounter) = size(foonode,1);
            NoBdyElems(idcounter) = size(fooele,1);
            idcounter = idcounter + 1;
    end
end
foo=zeros(size(mesh.nodes,1),1);
foo(1:no_ext_nodes,:) = 1;

mesh.bndvtx = foo;
mesh.region = [region1 region2];
