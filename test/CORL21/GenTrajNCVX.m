classdef GenTrajNCVX
	properties
		m = 0.01
		I = 1.2
		g = 0
		nv
		v
		regions
		polygons
		regions_init
		convexity
		env_regions
		planar = 1
		lines
		traj
		randv = 0
		ext_f
		env
	end

	methods
		function obj = GenTrajNCVX(init_,goal,nv)
			% number of vertices
			if nargin < 3; obj.nv = 4; end;

			if nv == 0
				obj.nv = 4;
				obj.v = [-1,-1;1,-1;1,1;-1,1]'/20;
			else
				if nv < 3; nv = 3; end; 
				obj.nv = nv
				% a = sort((rand([1,nv])-0.5)*2*pi);
				% x = (1/20) .* cos(a);
				% y = (1/20) .* sin(a);

				[x,y,dt] = simple_polygon(nv);
				x = x(1:end-1);
				y = y(1:end-1);

				x0 = mean(x);
				y0 = mean(y);
				x = x - x0;
				y = y - y0;

				obj.convexity = zeros(1,obj.nv)
			    obj.v = [x';y']*(1/(20*max([max(abs(x)), max(abs(y))])));
			end

			obj.convexity(1,convhull(obj.v')) = 1

			px = obj.v(1,:);
			py = obj.v(2,:);
			dt = delaunayTriangulation(px',py');   % delaunay triangulation 
	        points = dt.Points;    % points 
	        tri = dt.ConnectivityList;    % nodal connectivity
	        x = points(:,1);  
	        pX = x(tri);
	        y = points(:,2);
	        pY = y(tri);

	        X = [];
	        Y = [];
	        for i = 1:length(pX)
	            p1 = [pX(i,1);pY(i,1)];
	            p2 = [pX(i,2);pY(i,2)];
	            p3 = [pX(i,3);pY(i,3)];

	            p = (p1+p2+p3)/3;
	            % disp('fine')
	            if inpolygon(p(1),p(2),px,py)
	                X = [X;pX(i,:)];
	                Y = [Y;pY(i,:)];
	            end
	        end
	        x = px;
	        y = py;
	        % segments the object into triangles
	        obj.polygons = {};
	        M = size(X,1);
	        numvert = M*3;

	        for i = 1:M
	            obj.polygons{i} = struct('v', [], 'center', []);
	            v = [];
	            for j = 1:3
	                vert = [X(i,j);Y(i,j)];
	                v = [v, vert];
	            end
	            obj.polygons{i}.v = v;
	            obj.polygons{i}.center = [mean(v(1,:)),mean(v(2,:))]';
	        end
				
			dt = 0.1;
			obj.I = obj.m*0.01^2/3;

			% trajectory
			NT = 5;

			obj.traj.r = [linspace(init_(1),goal(1),5);linspace(init_(2),goal(2),5);linspace(init_(3),goal(3),5)];
			% obj.traj.r(1,:) = obj.traj.r(1,:) + rand(1, 5)/20-1/40;
			% obj.traj.r(2,:) = obj.traj.r(2,:) + rand(1, 5)/20-1/40;
			obj.traj.dr = zeros(3,NT);
			obj.traj.ddr = zeros(3,NT);

			% derivatives
			for t = 2:NT-1
				obj.traj.dr(:,t) = (obj.traj.r(:,t+1)-obj.traj.r(:,t-1))/(2*dt);
			end

			for t = 1:NT
				if t > 1; v1 = obj.traj.r(:,t-1); else; v1 = obj.traj.r(:,t); end;
				if t < NT; v2 = obj.traj.r(:,t+1); else; v2 = obj.traj.r(:,t); end;

				obj.traj.ddr(:,t) = (-2*obj.traj.r(:,t)+v1+v2)/(dt)^2;
			end

			% obj.traj.ddr
			% pause()

			% constraints for the entire plan
			obj.ext_f = {};
			for v = 1:obj.nv
				obj.ext_f{v} = struct();
				obj.ext_f{v}.fc1 = zeros(2,NT);
				obj.ext_f{v}.fc2 = zeros(2,NT);
				obj.ext_f{v}.jac = zeros(2,3,NT);

				for t = 2:NT-1
					% computes the jacobian and the velocity vector
					vs = [-sin(obj.traj.r(3,t)),-cos(obj.traj.r(3,t));cos(obj.traj.r(3,t)),-sin(obj.traj.r(3,t))]*obj.v(:,v);
					vel = [eye(2),vs]*obj.traj.dr(:,t);

					if v == 1
						vel(1) = vel(1);
					end

					obj.ext_f{v}.fc1(:,t) = -obj.m*9.8*normalize(vel)*0.3;
					obj.ext_f{v}.fc2(:,t) = -obj.m*9.8*normalize(vel)*0.3;
					obj.ext_f{v}.jac(:,:,t) = [eye(2),vs];
				end
			end

			% generates lines
			lines_init = {};

			% lines for reference
			for i = 1:obj.nv
				idx_1 = i;
				idx_2 = i+1;
				if idx_2 > obj.nv; idx_2 = 1; end;
				dv = obj.v(:,idx_2) - obj.v(:,idx_1);
				rot = [cos(angle(dv(1)+sqrt(-1)*dv(2))) , -sin(angle(dv(1)+sqrt(-1)*dv(2)));...
					   sin(angle(dv(1)+sqrt(-1)*dv(2))) ,  cos(angle(dv(1)+sqrt(-1)*dv(2)))];
				lines_init{i} = struct(); 
				lines_init{i}.v1 = obj.v(:,idx_1); 
				lines_init{i}.v2 = obj.v(:,idx_2);
				lines_init{i}.fc1 = rot*[-1;1]; 
				lines_init{i}.fc2 = rot*[1;1];
			end

			% lines for the entire plan
			obj.lines = {};

			% applies the transformation to each line segment
			for l = 1:obj.nv
				obj.lines{l} = struct();
				obj.lines{l}.v1 = zeros(2,NT);
				obj.lines{l}.v2 = zeros(2,NT);
				obj.lines{l}.fc1 = zeros(2,NT);
				obj.lines{l}.fc2 = zeros(2,NT);
				for t = 1:NT
					th = obj.traj.r(3,t); trans = obj.traj.r(1:2,t);
					rotmat = [cos(th),-sin(th);sin(th),cos(th)];
					obj.lines{l}.v1(:,t) = trans + rotmat*lines_init{l}.v1;
					obj.lines{l}.v2(:,t) = trans + rotmat*lines_init{l}.v2;

					obj.lines{l}.fc1(:,t) = rotmat*lines_init{l}.fc1;
					obj.lines{l}.fc2(:,t) = rotmat*lines_init{l}.fc2;
				end
			end

			% regions for reference
			% regions_init = {};
			% for i = 1:obj.nv
			% 	idx_1 = i;
			% 	idx_2 = i+1;
			% 	if idx_2 > obj.nv; idx_2 = 1; end;
			% 	regions_init{i} = struct(); 
			% 	res = inv([obj.v(:,idx_1)';obj.v(:,idx_2)'])*[1;1];
			% 	regions_init{i}.A = -res'; regions_init{i}.b = -1;
			% end 

			% convex hull regions
			k = convhull(obj.v');
			cvx_v = obj.v(:,k);
			cvx_v = cvx_v(:,1:end-1);

			obj.regions_init = {};
			for i = 1:size(cvx_v,2)
				idx_1 = i;
				idx_2 = i+1;
				if idx_2 > size(cvx_v,2); idx_2 = 1; end;
				obj.regions_init{i} = struct(); 
				res = inv([cvx_v(:,idx_1)';cvx_v(:,idx_2)'])*[1;1];
				obj.regions_init{i}.A = -res'; obj.regions_init{i}.b = -1;
			end

			% non-convex regions
			for i = 1:obj.nv
				if i == obj.nv; break; end
				indexes = []
				cvx = false
				k = i
				while cvx == false
					% indexes
					idx_0 = k-1; idx_1 = k; idx_2 = k+1;
					% 
					if idx_2 > obj.nv; idx_2 = 1; end;
					if idx_0 < 1; idx_0 = obj.nv; end;

					if obj.convexity(1,k) == 0
						if k == i
							indexes = [indexes, idx_0]
						end
						indexes = [indexes, k]
					else
						cvx = true
					end
					k = k + 1
					if k > obj.nv; k = 1; end
				end

				% creates region within non-convexity
				if length(indexes) > 0
					indexes = [indexes, idx_1]
					indexes
					obj.v(:,indexes)
					[A, b] = vert2con(obj.v(:,indexes)');
					obj.regions_init{end+1} = struct(); 
					obj.regions_init{end}.A = A; obj.regions_init{end}.b = b;
				end
			end

			% regions during the motion
			obj.regions = {};
			for r = 1:length(obj.regions_init)
				r
				obj.regions{r} = struct();
				obj.regions{r}.A = zeros(size(obj.regions_init{r}.A,1),size(obj.regions_init{r}.A,2),NT);
				obj.regions{r}.b = zeros(size(obj.regions_init{r}.b,1),NT);
				% reorients the planes through time
				for t = 1:NT
					th = obj.traj.r(3,t); trans = obj.traj.r(1:2,t);
					rotmat = [cos(th),-sin(th);sin(th),cos(th)];
					obj.regions{r}.A(:,:,t) = obj.regions_init{r}.A*inv(rotmat);
					obj.regions{r}.b(:,t) = obj.regions_init{r}.b + obj.regions_init{r}.A*inv(rotmat)*trans;
				end
			end	
		end
	end

end