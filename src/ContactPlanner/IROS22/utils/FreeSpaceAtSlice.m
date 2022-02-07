function free_cslice = FreeSpaceAtSlice(free_workspace,object)
	% inputs
	% 	- free_workspace: cell descring a non-convex compact region of the environment (cell of N line segments)
	% 	- object: cell descring a non-convex object (cell of M convex polygons)
	% outputs
	%	- free_cslice: vertices of compact regions from the free-cspace
	
	% performs sum for each pair of convex polygons
	pol_sum = {};
	for i = 1:length(free_workspace)
		for j = 1:length(object)
			pol_sum{end+1} = MinkowskiCVX2CVX(free_workspace{i},-object{j}.v);
		end
	end

	free_workspace_verts = [];
	for i = 1:length(free_workspace)
		free_workspace_verts = [free_workspace_verts, free_workspace{i}(:,2)];
	end

	free_cslice = polyshape(free_workspace_verts');


	% returns the vertex of the sum
	for i = 1:length(pol_sum)
		free_cslice = subtract(free_cslice,polyshape(pol_sum{i}(1,:),pol_sum{i}(2,:)));
	end

	% decomposes in compact regions
	if free_cslice.NumRegions > 1
		regions = free_cslice.regions();
		free_cslice = {};
		for i = 1:length(regions)
			free_cslice{end+1} = regions(i).Vertices';
		end
	else
		free_cslice = {free_cslice.Vertices'};
	end
end