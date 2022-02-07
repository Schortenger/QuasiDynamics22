classdef GenSagMotionRect
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
		lines
		traj
		randv = 0
		y_floor = 0.0
		th_floor = 0.0
		ext_f
		env
		ext_contact = zeros(1,5)
	end

	methods
		function obj = GenSagMotionRect(init_,goal,nv)
			% number of vertices
			if nargin < 3; obj.nv = 4; end;

			% polygon generation
			if nv == 0
				obj.nv = 4;
				obj.v = [-1,-1;1,-1;1,1;-1,1]'/20;
			else
				if nv < 3; nv = 3; end; 
				obj.nv = nv;
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

				obj.convexity = zeros(1,obj.nv);
			    obj.v = [x';y']*(1/(20*max([max(abs(x)), max(abs(y))])));
			end

			% polygon segmentation
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
			% obj.traj.r(1,:) = obj.traj.r(1,:) + 0*rand(1, 5)/20-1/40;
			% obj.traj.r(2,:) = obj.traj.r(2,:) + 0*rand(1, 5)/20-1/40;

			% floor angle and elevation
			th_floor = 0*rand();
			rot_floor = [cos(th_floor), sin(th_floor); -sin(th_floor), cos(th_floor)]
			y_floor  = 0.0*rand();

			obj.th_floor = th_floor;
			obj.y_floor = y_floor;


			% th_wall = pi/6*rand() + pi/2;
			% rot_wall = [cos(th_wall), sin(th_wall); -sin(th_wall), cos(th_wall)]
			% x_wall  = 0.05*rand() + 0.1;

			% floors for y <= 0
			obj.env_regions = {};
			obj.env_regions{1}.A = [sin(th_floor), -cos(th_floor)]; obj.env_regions{1}.b = -y_floor;	
			% obj.env_regions{2}.A = [sin(th_wall), -cos(th_wall)]; obj.env_regions{1}.b = x_floor;	
			
			obj.env = {}; 
			obj.env{1}.x = [-1,1,cos(th_floor),-cos(th_floor)]; obj.env{1}.y = [-1,-1, y_floor + sin(th_floor), y_floor - sin(th_floor)];
			obj.env{1}.n = [-sin(th_floor);cos(th_floor)];
			% obj.env{2}.x = [-1,1, x_wall + cos(th_wall),x_wall -cos(th_wall)]; obj.env{1}.y = [-1,-1, sin(th_wall), y_floor - sin(th_wall)];
			% obj.env{2}.n = [-sin(th_wall);cos(th_wall)];

			% corrects motion to account for envionment
			ext_mode = zeros(3,obj.nv,NT);

			for t = 1:NT
				% corrects interpenetration
				delta = 0;
				for v = 1:obj.nv
					pos_v = obj.traj.r(1:2,t) + [cos(obj.traj.r(3,t)),-sin(obj.traj.r(3,t));sin(obj.traj.r(3,t)),cos(obj.traj.r(3,t))]*obj.v(:,v);
					dist_v = obj.env{1}.n' * (pos_v - [0; y_floor]);

					if dist_v < 0 && (abs(dist_v) > delta)
						delta =  abs(dist_v);
					end
				end

				obj.traj.r(1:2,t) = obj.traj.r(1:2,t) + obj.env{1}.n*delta

				% checks if vertex is in contact with environment
				for v = 1:obj.nv
					pos_v = obj.traj.r(1:2,t) + [cos(obj.traj.r(3,t)),-sin(obj.traj.r(3,t));sin(obj.traj.r(3,t)),cos(obj.traj.r(3,t))]*obj.v(:,v);
					dist_v = obj.env{1}.n'*(pos_v - [0; y_floor])
					if abs(dist_v) < 1e-6
						ext_mode(2,v,t) = 1;
					end
				end

				% corrects interpenetration
				% delta = 0;
				% for v = 1:obj.nv
				% 	pos_v = obj.traj.r(1:2,t) + [cos(obj.traj.r(3,t)),-sin(obj.traj.r(3,t));sin(obj.traj.r(3,t)),cos(obj.traj.r(3,t))]*obj.v(:,v);
				% 	dist_v = obj.env{2}.n' * (pos_v - [x_wall; 0]);

				% 	if dist_v < 0 && (abs(dist_v) > delta)
				% 		delta =  abs(dist_v);
				% 	end
				% end

				% obj.traj.r(1:2,t) = obj.traj.r(1:2,t) + obj.env{2}.n*delta

				% % checks if vertex is in contact with environment
				% for v = 1:obj.nv
				% 	pos_v = obj.traj.r(1:2,t) + [cos(obj.traj.r(3,t)),-sin(obj.traj.r(3,t));sin(obj.traj.r(3,t)),cos(obj.traj.r(3,t))]*obj.v(:,v);
				% 	dist_v = obj.env{2}.n'*(pos_v - [x_wall; 0])
				% 	if abs(dist_v) < 1e-6
				% 		ext_mode(2,v,t) = 1;
				% 	end
				% end
			end

			obj.traj.dr = zeros(3,NT);
			obj.traj.ddr = zeros(3,NT);

			% derivatives
			obj.traj.dr(:,1) = (obj.traj.r(:,2)-obj.traj.r(:,1))/(2*dt);
			obj.traj.dr(:,NT) = (obj.traj.r(:,NT)-obj.traj.r(:,NT-1))/(2*dt);
			for t = 2:NT-1
				obj.traj.dr(:,t) = (obj.traj.r(:,t+1)-obj.traj.r(:,t-1))/(2*dt);
			end

			for t = 1:NT
				if t > 1; v1 = obj.traj.r(:,t-1); else; v1 = obj.traj.r(:,t); end;
				if t < NT; v2 = obj.traj.r(:,t+1); else; v2 = obj.traj.r(:,t); end;

				obj.traj.ddr(:,t) = (-2*obj.traj.r(:,t)+v1+v2)/(dt)^2;
			end

			% determines contact mode
			for t = 1:NT
				for v = 1:obj.nv
					% determines if it is sliding
					if t < NT
						% check if the vertex remain in contact between t and t+1
						if sum(ext_mode(:,v,t)) > 0 && sum(ext_mode(:,v,t+1)) > 0

							vs = [-sin(obj.traj.r(3,t)),-cos(obj.traj.r(3,t));cos(obj.traj.r(3,t)),-sin(obj.traj.r(3,t))]*obj.v(:,v);
							% orthogonalization to the floor normal
							vel = [eye(2),vs]*obj.traj.dr(:,t) - (obj.env{1}.n' * ([eye(2),vs]*obj.traj.dr(:,t)))*obj.env{1}.n;

							% only slides if the vertex is moving
							if vel(1) > 0 or vel(2) > 0
								ext_mode(1,v,t) = 0;
								ext_mode(2,v,t) = 0;
								ext_mode(3,v,t) = 1;
							end
							if vel(1) > 0 or vel(2) > 0
								ext_mode(1,v,t) = 1;
								ext_mode(2,v,t) = 0;
								ext_mode(3,v,t) = 0;
							end

							% vel = [eye(2),vs]*obj.traj.dr(:,t) - (obj.env{2}.n' * ([eye(2),vs]*obj.traj.dr(:,t)))*obj.env{2}.n;

							% % only slides if the vertex is moving
							% if vel(1) > 0 or vel(2) > 0
							% 	ext_mode(1,v,t) = 0;
							% 	ext_mode(2,v,t) = 0;
							% 	ext_mode(3,v,t) = 1;
							% end
							% if vel(1) > 0 or vel(2) > 0
							% 	ext_mode(1,v,t) = 1;
							% 	ext_mode(2,v,t) = 0;
							% 	ext_mode(3,v,t) = 0;
							% end
						end
					end
				end
			end

			% constraints for the entire plan
			obj.ext_f = {};
			for v = 1:obj.nv
				obj.ext_f{v} = struct();
				obj.ext_f{v}.fc1 = zeros(2,NT);
				obj.ext_f{v}.fc2 = zeros(2,NT);
				obj.ext_f{v}.jac = zeros(2,3,NT);

				for t = 1:NT
					% computes the jacobian and the velocity vector
					vs = [-sin(obj.traj.r(3,t)),-cos(obj.traj.r(3,t));cos(obj.traj.r(3,t)),-sin(obj.traj.r(3,t))]*obj.v(:,v);
					vel = [eye(2),vs]*obj.traj.dr(:,t);

					fc1 = (rot_floor*[-0.1;1])';
					fc2 = (rot_floor*[0.1;1])';

					if ext_mode(1,v,t) == 1; fc2 = fc1; end;
					if ext_mode(3,v,t) == 1; fc1 = fc2; end;

					if sum(ext_mode(:,v,t)) == 0
						fc1 = 0*[-0.1,1];
						fc2 = 0*[0.1,1];
					end

					obj.ext_f{v}.fc1(:,t) = fc1;
					obj.ext_f{v}.fc2(:,t) = fc2;
					obj.ext_f{v}.jac(:,:,t) = [eye(2),vs];

					if sum(ext_mode(:,v,t)) > 0
						obj.ext_contact(1,t) = v
					end

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
				indexes = [];
				cvx = false;
				k = i;
				while cvx == false
					% indexes
					idx_0 = k-1; idx_1 = k; idx_2 = k+1;
					% 
					if idx_2 > obj.nv; idx_2 = 1; end;
					if idx_0 < 1; idx_0 = obj.nv; end;

					if obj.convexity(1,k) == 0
						if k == i
							indexes = [indexes, idx_0];
						end
						indexes = [indexes, k];
					else
						cvx = true;
					end
					k = k + 1;
					if k > obj.nv; k = 1; end
				end

				% creates region within non-convexity
				if length(indexes) > 0
					indexes = [indexes, idx_1];
					obj.v(:,indexes);
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