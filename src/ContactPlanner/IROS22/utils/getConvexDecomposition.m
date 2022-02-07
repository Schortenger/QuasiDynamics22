function cvx_decomp = getConvexDecomposition(polygon)
	% inputs
	%	- object: set descring a non-convex polygon (2 x N matrix)
	% outputs
	%	- cvx_decomp: cell of convex polygons

	% performs sum for each pair of convex polygons
	tri = triangulation(polyshape(polygon'));
	con = tri.ConnectivityList;
	verts = tri.Points;
	
	cvx_decomp = struct();
	
	% generates polygons
	cvx_decomp.polys = {};
	for i = 1:size(con,1)
		% gets each polygon from the triangulation
		cvx_decomp.polys{end+1} = verts(con(i,:),:)';
	end

	% connectivity between polygons
	cvx_decomp.neig = zeros(size(con,1),size(con,1));
	% checks if segments of polygon intersect
	for i = 1:length(cvx_decomp.polys)
		for j = 1:length(cvx_decomp.polys)
			if i == j
				% do nothing
			else
				for k = 1:size(cvx_decomp.polys{j},2)
					pol = polyshape(cvx_decomp.polys{i}(1,:)',cvx_decomp.polys{i}(2,:)');
					if isinterior(pol,cvx_decomp.polys{j}(:,k)')
						cvx_decomp.neig(i,j) = 1;
						cvx_decomp.neig(j,i) = 1;    
					end
				end
			end
		end
	end
end