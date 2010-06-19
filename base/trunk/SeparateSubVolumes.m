function [tags extelem regionvols] = SeparateSubVolumes(telem,tnode)
% Assuming multiple material regions are delineated in 'telem', this routines
% returns:
% tags : list of nodes inside _every_ subvolume (one node per subvolume)
%        tags{:,1} = coordiantes
%        tags{:,2} = region ID
% 
% extelem : indices to the elements that compose the largest volume
% 
% regionvols: a cell containing subvolumes that share a material ID
%             regionvols{:,1} = material ID
%             regionvols{:,2} = ind_regions
%                               ind_regions{i} = indices of elements for
%                               every subvolume belonging to region ID
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
regionvols = cell(nregions,2);
tag_counter = 1;
maxvol = -realmax;
invc=0;
for i=1:nregions
    interior_nodes = [];
    bf=telem(:,4)==regions(i) | telem(:,5)==regions(i);
    subvol=telem(bf,1:3);
    % Get all the sub-surfaces of current region
    ind_regions = GetIndRegions(subvol,tnode);
    regionvols{i,1} = regions(i); % Storing the material ID of each region
    regionvols{i,2} = ind_regions;
    % Calculate an interior node for each sub-surface
    for j=1:size(ind_regions,1)
        % Orient each subvolume and get a point within its boundary
        fooelem = FixPatchOrientation(tnode,subvol(ind_regions{j},:),[],1);
        tmp = GetOneInteriorNode(fooelem,tnode);
        if isempty(tmp)
            invc=invc+1;
            fprintf('No of invalid subvolumes so far: %d\n',invc);
%             errordlg('Could not find an interior point in surface.','Mesh Error');
%             error('Could not find an interior point in surface.' );
        end
        interior_nodes = [interior_nodes; tmp];
        % Check if this surface is the exterior one, we need it for
        % actual checkerboard3d to work. We do this by calculating the
        % volume of each sub-volume.
        vol = shell_volume(subvol,tnode);
        if vol > maxvol
            maxvol = vol;
            extelem = subvol;
        end
    end
    for j=1:size(interior_nodes,1);
        tags{tag_counter,1} = interior_nodes(i,:);
        tags{tag_counter,2} = regions(i);
        tag_counter = tag_counter + 1;
    end
end
