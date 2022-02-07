classdef cube_picking
	properties
		m = 0.01
		I = 1
		g = 9.8
		nv
		v
		regions
		lines
		traj
		env_regions
		ext_f
		env
	end

	methods
		function obj = cube_picking()
			obj.nv = 8;
			l = 0.1;
			obj.v = [-l/2,-l/2,-l/2; l/2,-l/2,-l/2; l/2,l/2,-l/2; -l/2,l/2,-l/2;...
					 -l/2,l/2,l/2; l/2,l/2,l/2; l/2,-l/2,l/2; -l/2,-l/2,l/2]';
			dt = 0.1;
			obj.I = obj.m*0.01^2/3*eye(3);

			% trajectory
			obj.traj.r = [zeros(2,5);linspace(sqrt(l/2),sqrt(l),5).^2;zeros(3,5)];
			obj.traj.dr = zeros(6,5);
			obj.traj.ddr = zeros(6,5);
			NT = 5;

			% derivatives
			for t = 2:NT-1
				obj.traj.dr(:,t) = (obj.traj.r(:,t+1)-obj.traj.r(:,t-1))/(dt);
			end

			for t = 2:NT-1
				obj.traj.ddr(:,t) = (-2*obj.traj.r(:,t)+obj.traj.r(:,t-1)+obj.traj.r(:,t+1))/(dt)^2;
			end

			% external force constraints
			ext_f = {};
			ext_f{1} = struct(); ext_f{1}.fc1 = [1,1,1]; ext_f{1}.fc2 = [1,-1,1];
			ext_f{1}.fc3 = [-1,1,1]; ext_f{1}.fc4 = [-1,-1,1];
 
			ext_f{2} = struct(); ext_f{2}.fc1 = [1,1,1]; ext_f{2}.fc2 = [1,-1,1];
			ext_f{2}.fc3 = [-1,1,1]; ext_f{2}.fc4 = [-1,-1,1];

			ext_f{3} = struct(); ext_f{3}.fc1 = [1,1,1]; ext_f{3}.fc2 = [1,-1,1];
			ext_f{3}.fc3 = [-1,1,1]; ext_f{3}.fc4 = [-1,-1,1];

			ext_f{4} = struct(); ext_f{4}.fc1 = [1,1,1]; ext_f{4}.fc2 = [1,-1,1];
			ext_f{4}.fc3 = [-1,1,1]; ext_f{4}.fc4 = [-1,-1,1];

			ext_f{5} = struct(); ext_f{5}.fc1 = [0,0,0]; ext_f{5}.fc2 = [0,0,0];
			ext_f{5}.fc3 = [0,0,0]; ext_f{5}.fc4 = [0,0,0];

			ext_f{6} = struct(); ext_f{6}.fc1 = [0,0,0]; ext_f{6}.fc2 = [0,0,0];
			ext_f{6}.fc3 = [0,0,0]; ext_f{6}.fc4 = [0,0,0];

			ext_f{7} = struct(); ext_f{7}.fc1 = [0,0,0]; ext_f{7}.fc2 = [0,0,0];
			ext_f{7}.fc3 = [0,0,0]; ext_f{7}.fc4 = [0,0,0];

			ext_f{8} = struct(); ext_f{8}.fc1 = [0,0,0]; ext_f{8}.fc2 = [0,0,0];
			ext_f{8}.fc3 = [0,0,0]; ext_f{8}.fc4 = [0,0,0];

			% constraints for the entire plan
			obj.ext_f = {};
			for v = 1:8
				obj.ext_f{v} = struct();
				obj.ext_f{v}.fc1 = zeros(3,5);
				obj.ext_f{v}.fc2 = zeros(3,5);
				obj.ext_f{v}.fc3 = zeros(3,5);
				obj.ext_f{v}.fc4 = zeros(3,5);
				obj.ext_f{v}.jac = zeros(3,6,5);
				for t = 1
					obj.ext_f{v}.fc1(:,t) = ext_f{v}.fc1;
					obj.ext_f{v}.fc2(:,t) = ext_f{v}.fc2;
					obj.ext_f{v}.fc3(:,t) = ext_f{v}.fc3;
					obj.ext_f{v}.fc4(:,t) = ext_f{v}.fc4;
					obj.ext_f{v}.jac(:,:,t) = [eye(3),zeros(3,3)];
				end
			end

			l = 0.1;
			% lines for reference
			lines_init = {};

			lines_init{1} = struct(); 
			lines_init{1}.v = [l/2,-l/2,-l/2; l/2,-l/2,l/2; -l/2,-l/2,l/2; -l/2,-l/2,-l/2]'
			lines_init{1}.fc1 = [1,1,1]'; lines_init{1}.fc2 = [-1,1,1]'; 
			lines_init{1}.fc3 = [-1,1,-1]'; lines_init{1}.fc4 = [1,1,-1]';

			lines_init{2} = struct(); 
			lines_init{2}.v = [l/2,l/2,l/2; -l/2,l/2,l/2; -l/2,-l/2,l/2; l/2,-l/2,l/2]'
			lines_init{2}.fc1 = [1,1,-1]'; lines_init{2}.fc2 = [-1,1,-1]'; 
			lines_init{2}.fc3 = [1,-1,-1]'; lines_init{2}.fc4 = [-1,-1,-1]';

			lines_init{3} = struct(); 
			lines_init{3}.v = [l/2,l/2,l/2; l/2,l/2,-l/2; -l/2,l/2,-l/2; -l/2,l/2,l/2]'
			lines_init{3}.fc1 = [1,-1,1]'; lines_init{3}.fc2 = [-1,-1,1]'; 
			lines_init{3}.fc3 = [-1,-1,-1]'; lines_init{3}.fc4 = [1,-1,-1]';

			lines_init{4} = struct(); 
			lines_init{4}.v = [l/2,l/2,-l/2; l/2,-l/2,-l/2; -l/2,-l/2,-l/2; -l/2,l/2,-l/2]'
			lines_init{4}.fc1 = [1,1,1]'; lines_init{4}.fc2 = [-1,1,1]'; 
			lines_init{4}.fc3 = [-1,1,1]'; lines_init{4}.fc4 = [-1,-1,1]';

			lines_init{5} = struct(); 
			lines_init{5}.v = [l/2,l/2,l/2; l/2,l/2,-l/2; l/2,-l/2,-l/2; l/2,-l/2,l/2]'
			lines_init{5}.fc1 = [-1,1,1]'; lines_init{5}.fc2 = [-1,-1,1]'; 
			lines_init{5}.fc3 = [-1,1,-1]'; lines_init{5}.fc4 = [-1,-1,-1]';

			lines_init{6} = struct(); 
			lines_init{6}.v = [-l/2,l/2,l/2; -l/2,-l/2,l/2; -l/2,-l/2,-l/2; -l/2,l/2,-l/2]'
			lines_init{6}.fc1 = [1,1,1]'; lines_init{6}.fc2 = [1,-1,1]'; 
			lines_init{6}.fc3 = [1,1,-1]'; lines_init{6}.fc4 = [1,-1,-1]';
			
			% lines for the entire plan
			obj.lines = {};

			% applies the transformation to each facet
			for l = 1:6
				obj.lines{l} = struct();
				obj.lines{l}.nv = 4;
				obj.lines{l}.v = zeros(3,4,5);
				obj.lines{l}.fc1 = zeros(3,5);
				obj.lines{l}.fc2 = zeros(3,5);
				obj.lines{l}.fc3 = zeros(3,5);
				obj.lines{l}.fc4 = zeros(3,5);

				for t = 1:5
					th = obj.traj.r(4:6,t); 
					trans = obj.traj.r(1:3,t);
					rotmat = eul2rotm(th');
					for v = 1:4
						obj.lines{l}.v(:,v,t) = trans + rotmat*lines_init{l}.v(:,v);
					end
					obj.lines{l}.fc1(:,t) = rotmat*lines_init{l}.fc1;
					obj.lines{l}.fc2(:,t) = rotmat*lines_init{l}.fc2;
					obj.lines{l}.fc3(:,t) = rotmat*lines_init{l}.fc3;
					obj.lines{l}.fc4(:,t) = rotmat*lines_init{l}.fc4;
				end
			end
			l = 0.1;
			% regions for reference
			regions_init = {};
			regions_init{1} = struct(); regions_init{1}.A = [ 1, 0, 0];  regions_init{1}.b = [-l/2];
			regions_init{2} = struct(); regions_init{2}.A = [-1, 0, 0];  regions_init{2}.b = [-l/2]; 
			regions_init{3} = struct(); regions_init{3}.A = [ 0,1, 0];  regions_init{3}.b = [-l/2]; 
			regions_init{4} = struct(); regions_init{4}.A = [ 0,-1, 0];  regions_init{4}.b = [-l/2]; 
			regions_init{5} = struct(); regions_init{5}.A = [ 0, 0, 1];  regions_init{5}.b = [-l/2]; 
			regions_init{6} = struct(); regions_init{6}.A = [ 0, 0, -1];  regions_init{6}.b = [-l/2]; 

			% regions during the motion
			obj.regions = {};
			for r = 1:6
				obj.regions{r} = struct();
				obj.regions{r}.A = zeros(size(regions_init{r}.A,1),size(regions_init{r}.A,2),5);
				obj.regions{r}.b = zeros(size(regions_init{r}.b,1),5);
				% reorients the planes through time
				for t = 1:5
					trans = obj.traj.r(1:3,t);
					rotmat = eul2rotm(obj.traj.r(4:6,t)');
					obj.regions{r}.A(:,:,t) = regions_init{r}.A*inv(rotmat);
					obj.regions{r}.b(:,t) = regions_init{r}.b + regions_init{r}.A*inv(rotmat)*trans;
				end
			end		

			% environment regiions
			obj.env_regions = {}
			obj.env_regions{1}.A = [0,0,-1]; obj.env_regions{1}.b = -0.001;	
			obj.env = {}; 
			obj.env{1}.x = [-1,1,1,-1,-1,1,1,-1]; 
			obj.env{1}.y = [-1,-1,0,0,-1,-1,0,0];
			obj.env{1}.z = [1,1,1,1,-1,-1,-1,-1];
			obj.env{1}.n = [0;0;1];
		end
	end

end