function pol_sum = MinkowskiCVX2CVX(pol1,pol2)
	% inputs
	% 	- pol1: list of vertices of convex polygon 
	% 	- pol2: list of vertices of convex polygon
	% outputs
	%	- sum: minkowski sum of pol1 and pol2

	% sums all vertices
	pol_sum = [];
	for i = 1:size(pol1,2)
		for j = 1:size(pol2,2)
			pol_sum = [pol_sum, pol1(:,i)+pol2(:,j)];
		end
	end

	pol_sum = pol_sum(:,convhull(pol_sum'));
	pol_sum = pol_sum(:,1:end-1);
end