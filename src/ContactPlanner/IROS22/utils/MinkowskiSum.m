function free_cslice = FreeSpaceAtSlice(object,free_workspace_verts)
	% inputs
	% 	- free_workspace_verts: set of N vertices descring a non-convex compact region of the environment
	% 	- object: set descring a non-convex object (cell of M convex polygons)
	% outputs
	%	- free_cslice: vertices of compact regions from the free-cspace
	
	% checks if env is given
	free_workspace = {}
	for i = 1:size(free_workspace_verts,2)
		idx_1 = i+1;
		if i == size(free_workspace_verts,2); idx_1 = 1; end;;
		free_workspace{end+1} = [free_workspace_verts(:,i),free_workspace_verts(:,idx_1)];
	end
	
	% performs sum for each pair of convex polygons
	pol_sum = {};
	for i = 1:length(free_workspace)
		for j = 1:length(object)
			pol_sum{end+1} = MinkowskiCVX2CVX(free_workspace{i},object{j});
		end
	end

	free_cslice = polyshape(free_workspace_verts(1,:),free_workspace_verts(2,:));

	% returns the vertex of the sum
	for i = 1:length(pol_sum)
		free_cslice = substract(free_cslice,polyshape(pol_sum{i}(1,:),pol_sum{i}(2,:)));
	end

	% decomposes in compact regions
	if free_cslice.NumRegions > 1
		regions = free_cslice.regions()
		free_cslice = {};
		for i = 1:free_cslice.NumRegions
			free_cslice{end+1} = regions(i).Vertices;
		end
	else
		free_cslice = {free_cslice.Vertices};
	end
end