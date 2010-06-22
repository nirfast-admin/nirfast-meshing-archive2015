function [tags extelem] = SeparateSubVolumes(telem,tnode)
% Assuming multiple material regions are delineated in 'telem', this routines
% returns:
% tags : list of nodes inside _every_ subvolume (one node per subvolume)
%        tags{:,1} = coordiantes
%        tags{:,2} = region ID
% 
% extelem : indices to the elements that compose the largest volume
% 
% 
% telem(:,1:3) defines the triagnels in the surface
% telem(:,4:5) defines the two materials on each side of a triangle

if size(telem,2)~=3 && size(telem,2)<5
        error(['Input surface mesh does not seem to be a triangular mesh: ' fnprefix]);
end
if size(telem,2)==3
    nregions=1;
    regions=1;
    telem = [telem repmat([1 0], size(telem,1), 1)];
else
    regions = telem(:,4:5);
    regions = unique(regions(:));
    nregions = length(regions);
end
[tf idx]=ismember(0,regions);
if tf, regions(idx)=[]; nregions=nregions-1; end % Remove 0 (which represents outside space)


% Find a point within each region (note that regions might consist of
% multiple sub-surfaces)
tag_counter = 1;
maxvol = -realmax;
for i=1:nregions
    maxvol_local = -realmax;
    bf=telem(:,4)==regions(i) | telem(:,5)==regions(i);
    idx = find(bf);
    subvol=telem(bf,1:3);
    % Get all the sub-surfaces of current region
    ind_regions = GetIndRegions(subvol,tnode);
    % Figure out which sub-surface should be used for point coordinate
    % calculations.
    if size(ind_regions,1)>1
        for j=1:size(ind_regions,1)
            vol = shell_volume(subvol(ind_regions{j},1:3),tnode);
            if vol > maxvol_local
                maxvol_local = vol;
                myidx = j;
            end
        end
    else
        myidx = 1;
    end
    % Orient the surface
    fooelem = FixPatchOrientation(tnode,subvol(ind_regions{myidx},1:3),[],1);
    % Calculate a node within the surface
    interior_nodes = GetOneInteriorNode(fooelem,tnode);
    if isempty(interior_nodes)
            errordlg('Could not find an interior point in surface.','Mesh Error');
            error('Could not find an interior point in surface.' );
    end
    % Check if this surface is the exterior one, we need it for
    % actual checkerboard3d to work. We do this by calculating the
    % volume of each sub-volume.
    vol = shell_volume(fooelem,tnode);
    if vol > maxvol
        maxvol = vol;
        extelem = telem(idx(ind_regions{myidx}),1:3);
    end
    for j=1:size(interior_nodes,1);
        tags{tag_counter,1} = interior_nodes(j,:);
        tags{tag_counter,2} = regions(i);
        tag_counter = tag_counter + 1;
    end
end
