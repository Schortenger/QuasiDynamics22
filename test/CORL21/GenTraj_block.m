classdef GenTraj_block
	properties
		m = 0.01
		I = 1.2
		g = 0
		nv
		v
		regions
		polygons
		env_regions
		planar = 1
		lines
		traj
		randv = 0
		ext_f
		env
	end

	methods
		function obj = GenTraj_block(init_,goal,v,noise)
			% number of vertices
			if nargin < 4; noise = true; end

			obj.nv = 3;
			obj.v = [-1,-1;1,-1;1,1]'/20;
			dt = 0.1;
			obj.I = obj.m*0.01^2/3;

			obj.v = obj.v - [mean(obj.v(1,:)),mean(obj.v(2,:))]';

			% trajectory
			NT = 5;
			obj.polygons{1} = struct('v', obj.v, 'center', [mean(obj.v(1,:)),mean(obj.v(2,:))]');
			dt = 0.1;
			obj.I = obj.m*0.01^2/3;

			% trajectory
			NT = 5;
			obj.traj.r = [linspace(init_(1),goal(1),5);linspace(init_(2),goal(2),5);linspace(init_(3),goal(3),5)];

			if noise
				obj.traj.r(1,:) = obj.traj.r(1,:) + rand(1, 5)/20-1/40;
				obj.traj.r(2,:) = obj.traj.r(2,:) + rand(1, 5)/20-1/40;	
			end

			obj.traj.dr = zeros(3,NT);
			obj.traj.ddr = zeros(3,NT);

			% derivatives
			for t = 2:NT-1
				obj.traj.dr(:,t) = (obj.traj.r(:,t+1)-obj.traj.r(:,t-1))/(2*dt);
			end

			for t = 2:NT-1
				if t > 1; v1 = obj.traj.r(:,t-1); else; v1 = zeros(3,1); end;
				if t < NT; v2 = obj.traj.r(:,t+1); else; v2 = zeros(3,1); end;

				obj.traj.ddr(:,t) = (-2*obj.traj.r(:,t)+v1+v2)/(dt)^2;
			end

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
				rot = [cos(angle(dv(1)+j*dv(2))),-sin(angle(dv(1)+j*dv(2)));sin(angle(dv(1)+j*dv(2))),cos(angle(dv(1)+j*dv(2)))];
				lines_init{i} = struct(); lines_init{i}.v1 = obj.v(:,idx_1); lines_init{i}.v2 = obj.v(:,idx_2);
				lines_init{i}.fc1 = rot*[-0.1;1]; lines_init{i}.fc2 = rot*[0.1;1];
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
			regions_init = {};
			for i = 1:obj.nv
				idx_1 = i;
				idx_2 = i+1;
				if idx_2 > obj.nv; idx_2 = 1; end;
				regions_init{i} = struct(); 
				res = inv([obj.v(:,idx_1)';obj.v(:,idx_2)'])*[1;1];
				regions_init{i}.A = -res'; regions_init{i}.b = -1;
			end 
			% regions during the motion
			obj.regions = {};
			for r = 1:obj.nv
				obj.regions{r} = struct();
				obj.regions{r}.A = zeros(size(regions_init{r}.A,1),size(regions_init{r}.A,2),NT);
				obj.regions{r}.b = zeros(size(regions_init{r}.b,1),NT);
				% reorients the planes through time
				for t = 1:NT
					th = obj.traj.r(3,t); trans = obj.traj.r(1:2,t);
					rotmat = [cos(th),-sin(th);sin(th),cos(th)];
					obj.regions{r}.A(:,:,t) = regions_init{r}.A*inv(rotmat);
					obj.regions{r}.b(:,t) = regions_init{r}.b + regions_init{r}.A*inv(rotmat)*trans;
				end
			end	
		end
	end

end