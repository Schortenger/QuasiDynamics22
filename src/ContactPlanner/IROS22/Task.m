classdef Task
	properties
	    m = 0.01
		I = 0.1
		g = 9.8
		nv
		n_env_v
		v
		regions
		polygons
		convexity
		environment
		regions_init
		env_regions
		lines
		traj
		randv = 0
		y_floor = 0.0
		th_floor = 0.0
		ext_f
		v_f
		env
		ext_contact
		v_contact
		mu_f = 0.1;
		mu = 0.1;
		planar = false;
	end

	methods
		function obj = Task(object, verts, env, r, dr_prev, ddr_prev, planar)

			% if nargin < 6
			% 	planar = false;
			% end

			% obj.planar = planar;

			if nargin < 5
				dr_prev = zeros(3,5);
				ddr_prev = zeros(3,5);
			end

			obj.v = -verts;
			obj.convexity(1,convhull(obj.v')) = 1

			obj.nv = size(verts,2);

			environment = zeros(2,length(env)+1);
			obj.n_env_v = length(env);
			obj.environment = environment

			for i = 1:length(env)
				environment(:,i) = env{i}(:,1);
			end
			environment(:,end) = env{end}(:,2);

			obj.polygons = {};

			for i = 1:length(object)
	            obj.polygons{i} = struct('v', [], 'center', []);
	            v = -object{i}.v;
	            obj.polygons{i}.v = v;
	            obj.polygons{i}.center = [mean(v(1,:)),mean(v(2,:))]';
	        end

	        dt = 0.1;
			obj.I = obj.m*0.25^2/3;

			NT = size(r,2);
			obj.ext_contact = zeros(1,NT);

			obj.traj.r = r;
			obj.traj.dr = zeros(3,NT);
			obj.traj.ddr = zeros(3,NT);

			% derivatives
			obj.traj.dr(:,1) = (-3*obj.traj.r(:,1)+4*obj.traj.r(:,2)-obj.traj.r(:,3))/(2*dt);
			for t = 2:NT-1
				obj.traj.dr(:,t) = (obj.traj.r(:,t+1)-obj.traj.r(:,t-1))/(2*dt);
			end
			obj.traj.dr(:,end) = (3*obj.traj.r(:,end)-4*obj.traj.r(:,end-1)+obj.traj.r(:,end-2))/(2*dt);

			obj.traj.ddr(:,1) = (-3*obj.traj.dr(:,1)+4*obj.traj.dr(:,2)-obj.traj.dr(:,3))/(2*dt)
			for t = 2:NT-1
				obj.traj.ddr(:,t) = (obj.traj.dr(:,t+1) - obj.traj.dr(:,t-1))/(2*dt);
			end
			obj.traj.ddr(:,end) = (3*obj.traj.dr(:,end)-4*obj.traj.dr(:,end-1)+obj.traj.dr(:,end-2))/(2*dt);

			obj.traj.ddr(3,:) = -obj.traj.ddr(3,:);

			% environment regions
			env_decomp = getConvexDecomposition(environment);
			obj.env_regions = {};

			for i = 1:length(env_decomp.polys)
				poly1 = env_decomp.polys{i};
				poly_center = mean(poly1,2);
				poly_fixed = poly_center + (poly1 - poly_center)*0.95;
				[A,b] = vert2con(poly_fixed');
				obj.env_regions{i} = struct('A',A,'b',b);
			end

			% determines contact mode
			ext_mode = zeros(3,obj.nv,NT);
			
			% determines external friction cones
			obj.ext_f = {};
			for v = 1:obj.nv
				obj.ext_f{v} = struct();
				obj.ext_f{v}.fc1 = zeros(2,NT);
				obj.ext_f{v}.fc2 = zeros(2,NT);
				obj.ext_f{v}.jac = zeros(2,3,NT);
				obj.ext_f{v}.vel = zeros(2,NT);
				obj.ext_f{v}.external = zeros(1,NT);

				for t = 1:NT
					% computes the jacobian and the velocity vector
					vs = [-sin(obj.traj.r(3,t)),-cos(obj.traj.r(3,t));cos(obj.traj.r(3,t)),-sin(obj.traj.r(3,t))]*obj.v(:,v);
					vel = [eye(2),vs]*obj.traj.dr(:,t); 

					v_w = obj.traj.r(1:2,t) + [cos(obj.traj.r(3,t)),-sin(obj.traj.r(3,t));sin(obj.traj.r(3,t)),cos(obj.traj.r(3,t))]*obj.v(:,v);
					[flag, n_env] = obj.inContact(v_w, environment);
					% flag
					% n_env
					% pause()
					if obj.planar
						obj.ext_f{v}.fc1(:,t) = -obj.mu_f*vel/(norm(vel) + 1e-6);
						obj.ext_f{v}.fc2(:,t) = -obj.mu_f*vel/(norm(vel) + 1e-6);
						obj.ext_f{v}.jac(:,:,t) = [eye(2),vs];
					end

					if flag
						v_floor = environment(:,n_env + 1) - environment(:,n_env);
						v_wrt_floor = vel'*v_floor;

						if v_wrt_floor > 0
							ext_mode(1,v,t) = 1;
						end
						if v_wrt_floor < 0
							ext_mode(3,v,t) = 1;
						end
						if v_wrt_floor == 0
							ext_mode(2,v,t) = 1;
						end

						% ext_mode(2,v,t) = 1;

						th_floor = angle(v_floor(1)+1i*v_floor(2));
						rot_floor = [cos(th_floor),-sin(th_floor);sin(th_floor),cos(th_floor)];

						fc1 = (rot_floor*[-obj.mu_f;1])';
						fc2 = (rot_floor*[obj.mu_f;1])';


						if ext_mode(1,v,t) == 1; fc2 = fc1; end;
						if ext_mode(3,v,t) == 1; fc1 = fc2; end;

						obj.ext_f{v}.fc1(:,t) = fc1;
						obj.ext_f{v}.fc2(:,t) = fc2;
						obj.ext_f{v}.jac(:,:,t) = [eye(2),vs];
						obj.ext_f{v}.external(1,t) = 1;
					end

					if sum(ext_mode(:,v,t)) > 0
						obj.ext_contact(1,t) = v;
					end
				end
			end

			% determines environment vertex friction cones
			v_mode = zeros(3,size(env,2),NT);
			obj.v_f = {};
			for v = 1:obj.n_env_v
				obj.v_f{v} = struct();
				obj.v_f{v}.fc1 = zeros(2,NT);
				obj.v_f{v}.fc2 = zeros(2,NT);
				obj.v_f{v}.jac = zeros(2,3,NT);

				for t = 1:NT
					% computes the jacobian and the velocity vector
					% vs = [-sin(obj.traj.r(3,t)),-cos(obj.traj.r(3,t));cos(obj.traj.r(3,t)),-sin(obj.traj.r(3,t))]*obj.v(:,v);
					% vel = [eye(2),vs]*obj.traj.dr(:,t); 

					v_w = obj.traj.r(1:2,t) + [cos(obj.traj.r(3,t)),-sin(obj.traj.r(3,t));sin(obj.traj.r(3,t)),cos(obj.traj.r(3,t))]*obj.v;
					[flag, n_v_0] = obj.inContact(environment(:,v), v_w);
					% flag
					% pause()
					if flag

						n_v_1 = n_v_0 + 1;
						if n_v_0 == size(v_w,2); n_v_1 = 1; end;

						v_obj = v_w(:,n_v_1) - v_w(:,n_v_0);

						th_f = angle(v_obj(1)+sqrt(-1)*v_obj(2));
						rot = [cos(th_f) , -sin(th_f);...
							   sin(th_f) ,  cos(th_f)];

						fc1 = (rot*[-obj.mu;1])';
						fc2 = (rot*[obj.mu;1])';

						if obj.traj.dr(1:2,t)'*	v_obj > 0; v_mode(1,v,t) = 1; end;
						if obj.traj.dr(1:2,t)'*	v_obj < 0; v_mode(3,v,t) = 1; end;

						if v_mode(1,v,t) == 1; fc2 = fc1; end;
						if v_mode(3,v,t) == 1; fc1 = fc2; end;

						% v
						% v_obj
						% obj.traj.dr(1:2,t)
						% fc1
						% fc2

						% figure(420)
						% pgon1 = polyshape(environment(1,:)',environment(2,:)');
						% plot(pgon1,'FaceAlpha',0.2,'FaceColor',[100 100 100]/255,'EdgeColor','black')
						% hold on;

						% pgon2 = polyshape(v_w(1,:)',v_w(2,:)');
						% plot(pgon2,'FaceAlpha',0.9,'FaceColor',[221 217 195]/255,'EdgeColor','black')
						% hold on;

						% plot(environment(1,v),environment(2,v),'r*')

						% xlim([-1,1]);
						% ylim([-0.5,1]);

						% pause()

						% close all;

						obj.v_f{v}.fc1(:,t) = fc1;
						obj.v_f{v}.fc2(:,t) = fc2;
						% obj.v_f{v}.jac(:,:,t) = [eye(2),vs];
					end

					if sum(v_mode(:,v,t)) > 0
						obj.v_contact(1,t) = v;
					end
				end
			end

			% determines object facets and friction cones
			lines_init = {};

			% lines for reference
			for i = 1:obj.nv
				idx_1 = i;
				idx_2 = i+1;
				if idx_2 > obj.nv; idx_2 = 1; end;
				dv = obj.v(:,idx_2) - obj.v(:,idx_1);

				th_f = angle(dv(1)+sqrt(-1)*dv(2));
				rot = [cos(th_f) , -sin(th_f);...
					   sin(th_f) ,  cos(th_f)];
				
				lines_init{i} = struct(); 
				lines_init{i}.v1 = obj.v(:,idx_1); 
				lines_init{i}.v2 = obj.v(:,idx_2);
				lines_init{i}.fc1 = rot*[-obj.mu;1]; 
				lines_init{i}.fc2 = rot*[ obj.mu;1];
			end

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

			% determines free-space regions
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

		function [flag, facet] = inContact(obj, vertex, polygon)
			flag = 0;
			facet = -1;

			for i = 1:size(polygon,2)-1
				
				x = vertex'; %some point
				a_ = polygon(:,i)'; %segment points a,b
				b_ = polygon(:,i+1)';

				a = 0.9*a_ + 0.1*b_;
				b = 0.1*a_ + 0.9*b_;

				d_ab = norm(a-b);
				d_ax = norm(a-x);
				d_bx = norm(b-x);

				if dot(a-b,x-b)*dot(b-a,x-a)>=0
				    A = [a,1;b,1;x,1];
				    dista = abs(det(A))/d_ab;        
				else
				    dista = min(d_ax, d_bx);
				end

				if dista < 1e-2
					flag = 1;
					facet = i;
				end
			end
		end
	end
end