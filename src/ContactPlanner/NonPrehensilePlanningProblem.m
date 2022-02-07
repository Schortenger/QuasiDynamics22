classdef NonPrehensilePlanningProblem
	properties
	    object
	    vertices
	    options
	    manipulator
	    environment
	    M = 6
	    dt = 0.1
	    c_space
	    nodes
	    one_finger = false
	    add_cost = false
	    slices = 0
	    yumi_kin = true
	    L = []
	    R = []
	    node_path = {}
	end

	methods

		function obj = NonPrehensilePlanningProblem(object,environment,options)
			if nargin < 3
				options = struct('M',6,'dt',0.1,'slices',linspace(-pi/2,pi/2,9));
			end

			obj.object = object;
			obj.environment = environment;
			obj.options = options;
			obj.slices = options.slices
		end

		function obj = buildCspaceMap(obj)
			% constructs the c-slices map
			obj.c_space = {};
			for i = 1:size(obj.options.slices,2)
				% object shapes
				th = obj.options.slices(i);
				rot = [cos(th),-sin(th);sin(th),cos(th)];

				object = {};

				for j = 1:length(obj.object)
					object{end+1} = struct('v',-rot*obj.object{j}.v);
				end
				
				free_space = FreeSpaceAtSlice(obj.environment,object);

				polys = {};
				connectivity = [];

				% decomposes into convex polygons
				re = 0;
				for j = 1:length(free_space)
					cvx_decomp = getConvexDecomposition(free_space{j});
					
					for k = 1:length(cvx_decomp.polys)
						polys{end+1} = cvx_decomp.polys{k};
					end
					if re == 1;
						connectivity_1 = [connectivity, zeros(size(connectivity,1),size(cvx_decomp.neig,2))];
						connectivity_2 = [zeros(size(cvx_decomp.neig,1),size(connectivity,2)),cvx_decomp.neig];
						connectivity = [connectivity_1; connectivity_2];
					else
						connectivity = [cvx_decomp.neig];
						re = 1;
					end
				end

				% builds stucture
				obj.c_space{end+1} = struct('polygons',struct(),'connectivity',connectivity);
				obj.c_space{end}.polygons.i = zeros(1,length(polys));
				obj.c_space{end}.polygons.v = polys;
			end

			% checks connectivity between adjacent slices
			for i = 1:size(obj.options.slices,2)-1
				connectivity = zeros(length(obj.c_space{i}.polygons.v),length(obj.c_space{i+1}.polygons.v));
				for k = 1:length(obj.c_space{i}.polygons.v)
					for l = 1:length(obj.c_space{i+1}.polygons.v)
						if overlaps(polyshape(obj.c_space{i}.polygons.v{k}'),polyshape(obj.c_space{i+1}.polygons.v{l}'))
							connectivity(k,l) = 1;
						end
					end
				end
				% stores connectivity
				obj.c_space{i}.fwd_connectivity = connectivity;
				obj.c_space{i+1}.bwd_connectivity = connectivity';
			end

			obj.nodes = {};
			l = 1;
			for i = 1:length(obj.c_space)
				for j = 1:length(obj.c_space{i}.polygons.v)
					obj.c_space{i}.polygons.i(j) = l;
					l = l + 1;
				end
			end

			for i = 1:length(obj.c_space)
				for j = 1:length(obj.c_space{i}.polygons.v)
					% for each node
					node = struct('center',[0;0;0],'v', obj.c_space{i}.polygons.v{j},'f',0,'prnt',0,'succesors',0, 'feasibile', 1);
					for v = 1:size(obj.c_space{i}.polygons.v{j},2)
						x = obj.c_space{i}.polygons.v{j}(:,v)/size(obj.c_space{i}.polygons.v{j},2);
						node.center = node.center + [x;0];
					end

					node.center(3) = obj.options.slices(i);
					node.succesors = {};

					for k = 1:size(obj.c_space{i}.connectivity,2)
						if obj.c_space{i}.connectivity(j,k) == 1
							node.succesors{end+1} = obj.c_space{i}.polygons.i(k);
						end 
					end

					if i > 1
						for k = 1:size(obj.c_space{i}.bwd_connectivity,2)
							if obj.c_space{i}.bwd_connectivity(j,k) == 1
								node.succesors{end+1} = obj.c_space{i-1}.polygons.i(k);
							end
						end
					end

					if i < length(obj.c_space)
						for k = 1:size(obj.c_space{i}.fwd_connectivity,2)
							if obj.c_space{i}.fwd_connectivity(j,k) == 1
								node.succesors{end+1} = obj.c_space{i+1}.polygons.i(k);
							end 
						end
					end

					obj.nodes{end+1} = node;
				end
			end
		end

		function explored = searchRegionsDP(obj, q_start, q_goal)
			% finds regions
			start_node = -1;
			goal_node  = -1;
			for i = 1:length(obj.nodes)
				% vertives
				xv = obj.nodes{i}.v(1,:);
				yv = obj.nodes{i}.v(2,:);
				
				% finds finds the beginning and end points
				if (abs(q_start(3) - obj.nodes{i}.center(3)) < 1e-4 && inpolygon(q_start(1), q_start(2), xv, yv))
					start_node = i; 
				end
				
				if (abs(q_goal(3) - obj.nodes{i}.center(3)) < 1e-4 && inpolygon(q_goal(1), q_goal(2), xv, yv))
					goal_node  = i; 
				end
			end

			if start_node == -1 
				start_node = 1;
				disp('start not found') 
			end
			if goal_node == -1
				goal_node = length(obj.nodes);
				disp('goal not found')
			end

			start_node
			goal_node

			% assigns score to each node
			for i = 1:length(obj.nodes)
				f = obj.metric(obj.nodes{i},obj.nodes{goal_node});
				if size(find([obj.nodes{i}.succesors{:}] == goal_node),2) > 0
					obj.nodes{i}.connects = 1;
					obj.nodes{i}.f = f;
				else
					obj.nodes{i}.connects = 0;
					obj.nodes{i}.f = inf;
				end
			end

			obj.nodes{goal_node}.f = 0;
			obj.nodes{goal_node}.connects = 1;

			% value iteration
			for iter = 1:length(obj.nodes)
				for i = 1:length(obj.nodes)
					if obj.nodes{i}.connects == 1
						for j = 1:length(obj.nodes{i}.succesors)
							f_new = obj.nodes{i}.f + obj.metric(obj.nodes{i},obj.nodes{obj.nodes{i}.succesors{j}});
							if obj.nodes{obj.nodes{i}.succesors{j}}.connects == 1
								if obj.nodes{obj.nodes{i}.succesors{j}}.f > f_new
									obj.nodes{obj.nodes{i}.succesors{j}}.f = f_new;
								end
							else
								obj.nodes{obj.nodes{i}.succesors{j}}.connects = 1;
								obj.nodes{obj.nodes{i}.succesors{j}}.f = f_new;
							end
						end
					end
				end

				% if obj.nodes{start_node}.connects == 1
				% 	break;
				% end
			end

			if obj.nodes{start_node}.connects == 0
				disp('no connection');
				explored = {start_node};
				return
			end

			explored = {};
			current_node = start_node;

			% fins the shortest path
			while true
				if size(find([explored{:}] == current_node),2) <= 0
					explored{end+1} = current_node;
				end

				if current_node == goal_node
					break;
				end

				f_min = inf;
				next_node = 0;
				for i = 1:length(obj.nodes{current_node}.succesors)
					if obj.nodes{obj.nodes{current_node}.succesors{i}}.f < f_min
						if size(find([explored{:}] == obj.nodes{current_node}.succesors{i}),2) <= 0
							f_min = obj.nodes{obj.nodes{current_node}.succesors{i}}.f;
							next_node = obj.nodes{current_node}.succesors{i};
						end
					end
				end
				next_node
				
				current_node = next_node;
			end
		end

		function [q, p, f, f_ext] = solvePlanHierarchical(obj, q_start, q_goal)
			% finds regions
			start_node = -1;
			goal_node  = -1;
			for i = 1:length(obj.nodes)
				% vertives
				xv = obj.nodes{i}.v(1,:);
				yv = obj.nodes{i}.v(2,:);
				
				% finds finds the beginning and end points
				if (q_start(3) == obj.nodes{i}.center(3) && inpolygon(q_start(1), q_start(2), xv, yv))
					start_node = i; 
				end
				
				if (q_goal(3)  == obj.nodes{i}.center(3) && inpolygon(q_goal(1), q_goal(2), xv, yv))
					goal_node  = i; 
				end
			end

			if start_node == -1 
				start_node = 1;
				disp('start not found') 
				pause()
			end
			if goal_node == -1
				goal_node = length(obj.nodes);
				disp('goal not found')
				pause()
			end

			start_node
			goal_node

			% assigns score to each node
			for i = 1:length(obj.nodes)
				f = obj.metric(obj.nodes{i},obj.nodes{goal_node});
				if size(find([obj.nodes{i}.succesors{:}] == goal_node),2) > 0
					obj.nodes{i}.connects = 1;
					obj.nodes{i}.f = f;
				else
					obj.nodes{i}.connects = 0;
					obj.nodes{i}.f = inf;
				end
			end

			obj.nodes{goal_node}.f = 0;
			obj.nodes{goal_node}.connects = 1;

			% value iteration
			for iter = 1:length(obj.nodes)
				for i = 1:length(obj.nodes)
					if obj.nodes{i}.connects == 1
						for j = 1:length(obj.nodes{i}.succesors)
							f_new = obj.nodes{i}.f + obj.metric(obj.nodes{i},obj.nodes{obj.nodes{i}.succesors{j}});
							if obj.nodes{obj.nodes{i}.succesors{j}}.connects == 1
								if obj.nodes{obj.nodes{i}.succesors{j}}.f > f_new
									obj.nodes{obj.nodes{i}.succesors{j}}.f = f_new;
								end
							else
								obj.nodes{obj.nodes{i}.succesors{j}}.connects = 1;
								obj.nodes{obj.nodes{i}.succesors{j}}.f = f_new;
							end
						end
					end
				end

				% if obj.nodes{start_node}.connects == 1
				% 	break;
				% end
			end

			if obj.nodes{start_node}.connects == 0
				disp('no connection');
				explored = {start_node};
				return
			end

			explored = {};
			current_node = start_node;
			q_current = q_start;
			q = [];
			p = [];
			f = [];
			f_ext = [];
			R = [];
			L = [];
			lam = [];
			t_last = {};
			prev_t = 0;
			t = 0;

			% fins the shortest mechanically feasible path
			while true
				fprintf('\n\n\n\n\n\n STARTED IT \n\n\n\n\n\n')
				t
				if size(find([explored{:}] == current_node),2) <= 0
					explored{end+1} = current_node;
					t_last{end+1} = prev_t;
				end

				f_min = inf;
				next_node = 0;
				if current_node == goal_node
					if t > 0
						R_pass = R(:,:,:,end-t_last{end}+1:end);
						L_pass = L(:,:,:,end-t_last{end}+1:end);
						lam_pass = lam(:,:,:,end-t_last{end}+1:end);
						[q_, p_, f_, f_ext_, L_, R_, lam_, dp_, ddp_, flag] = solveTrajAndContactsNode(obj, q_current, current_node, current_node, L_pass, R_pass, lam_pass, dp_prev, ddp_prev, q_goal);
					else
						[q_, p_, f_, f_ext_, L_, R_, lam_, dp_, ddp_, flag] = solveTrajAndContactsNodeWithGoal(obj, q_current, current_node, current_node,q_goal);
					end

					if flag
						next_node = current_node;
						q_prev = q_;
						L_prev = L_;
						R_prev = R_;
						p_prev = p_;
						f_prev = f_;
						dp = dp_;
						ddp = ddp_;
						f_ext_prev = f_ext_;
						lam_prev = lam_;
					end
				else
					for i = 1:length(obj.nodes{current_node}.succesors)
						if obj.nodes{obj.nodes{current_node}.succesors{i}}.f < f_min
							if size(find([explored{:}] == obj.nodes{current_node}.succesors{i}),2) <= 0		
								if t > 0
									R_pass = R(:,:,:,end-t_last{end}+1:end);
									L_pass = L(:,:,:,end-t_last{end}+1:end);
									lam_pass = lam(:,:,:,end-t_last{end}+1:end);
									if current_node == goal_node
										[q_, p_, f_, f_ext_, L_, R_, lam_, dp_, ddp_, flag] = solveTrajAndContactsNode(obj, q_current, current_node, obj.nodes{current_node}.succesors{i}, L_pass, R_pass, lam_pass, dp_prev, ddp_prev, q_goal);
									else
										[q_, p_, f_, f_ext_, L_, R_, lam_, dp_, ddp_, flag] = solveTrajAndContactsNode(obj, q_current, current_node, obj.nodes{current_node}.succesors{i}, L_pass, R_pass, lam_pass, dp_prev, ddp_prev);
									end
								else
									[q_, p_, f_, f_ext_, L_, R_, lam_, dp_, ddp_, flag] = solveTrajAndContactsNode(obj, q_current, current_node, obj.nodes{current_node}.succesors{i});
								end

								if flag
									f_min = obj.nodes{obj.nodes{current_node}.succesors{i}}.f;
									next_node = obj.nodes{current_node}.succesors{i};
									q_prev = q_;
									L_prev = L_;
									R_prev = R_;
									p_prev = p_;
									f_prev = f_;
									dp = dp_;
									ddp = ddp_;
									f_ext_prev = f_ext_;
									lam_prev = lam_;
								end
							end
						end
					end
				end

				if next_node == 0
					next_node = explored{end-1};
					explored(end) = [];
					obj.nodes{current_node}.f = inf;

					q = q(:,1:end-t_last{end});
					p = p(:,:,:,1:end-t_last{end});
					f = f(:,:,:,1:end-t_last{end});
					f_ext = f_ext(:,:,1:end-t_last{end});
					R = R(:,:,:,1:end-t_last{end});
					L = L(:,:,:,1:end-t_last{end});
					lam = lam(:,:,:,1:end-t_last{end});
					t_last(end) = [];

					if length(explored) == 0
						t = 0;
					end
					t = t - 1;
				else
					if t == 0
						q = q_prev;
						p = p_prev;
						f = f_prev;
						f_ext = f_ext_prev;

						R = R_prev;
						L = L_prev;
						lam = lam_prev;

						prev_t = size(p_prev,4);
					else 
						q = cat(2, q, q_prev(:,2:end));
						p = cat(4, p, p_prev(:,:,:,2:end));
						f = cat(4, f, f_prev(:,:,:,2:end));
						f_ext = cat(3, f_ext, f_ext_(:,:,2:end));

						R = cat(4, R, R_prev(:,:,:,2:end));
						L = cat(4, L, L_prev(:,:,:,2:end));
						lam = cat(4, lam, lam_prev(:,:,:,2:end));

						prev_t = size(p_prev(:,:,:,2:end),4);
					end
					t = t + 1;
					dp_prev = dp;
					ddp_prev = ddp;
				end

				q_current = q(:,end)';
				
				if current_node == goal_node
					break;
				end

				current_node = next_node;
				fprintf('\n\n\n\n\n\n DID IT \n\n\n\n\n\n')
			end
		end

		function [q, p, f, f_ext, f_v, L, R, explored] = solvePlanBacktrack(obj, q_start, q_goal)
			% finds regions
			start_node = -1;
			goal_node  = -1;
			for i = 1:length(obj.nodes)
				% vertives
				xv = obj.nodes{i}.v(1,:);
				yv = obj.nodes{i}.v(2,:);
				
				% finds finds the beginning and end points
				if (abs(q_start(3) - obj.nodes{i}.center(3)) < 1e-1 && inpolygon(q_start(1), q_start(2), xv, yv))
					start_node = i; 
				end
				
				if (abs(q_goal(3) - obj.nodes{i}.center(3)) < 1e-1 && inpolygon(q_goal(1), q_goal(2), xv, yv))
					goal_node  = i; 
				end

				abs(q_goal(3) - obj.nodes{i}.center(3))
			end

			if start_node == -1 
				start_node = 1;
				disp('start not found') 
				pause()
			end
			if goal_node == -1
				goal_node = length(obj.nodes);
				disp('goal not found')
				pause()
			end

			if goal_node == start_node
				disp('start and goal are the same node, beware')
				pause()
			end

			start_node
			goal_node

			% assigns score to each node
			for i = 1:length(obj.nodes)
				fun = obj.metric(obj.nodes{i},obj.nodes{goal_node});
				if size(find([obj.nodes{i}.succesors{:}] == goal_node),2) > 0
					obj.nodes{i}.connects = 1;
					obj.nodes{i}.f = fun;
				else
					obj.nodes{i}.connects = 0;
					obj.nodes{i}.f = inf;
				end
			end

			obj.nodes{goal_node}.f = 0;
			obj.nodes{goal_node}.connects = 1;

			% value iteration
			max_iters = 100*length(obj.nodes);
			for iter = 1:max_iters
				old_nodes = obj.nodes;
				for i = 1:length(obj.nodes)
					if obj.nodes{i}.connects == 1
						for j = 1:length(obj.nodes{i}.succesors)
							f_new = obj.nodes{i}.f + obj.metric(obj.nodes{i},obj.nodes{obj.nodes{i}.succesors{j}});
							if obj.nodes{obj.nodes{i}.succesors{j}}.connects == 1
								if obj.nodes{obj.nodes{i}.succesors{j}}.f > f_new
									obj.nodes{obj.nodes{i}.succesors{j}}.f = f_new;
								end
							else
								obj.nodes{obj.nodes{i}.succesors{j}}.connects = 1;
								obj.nodes{obj.nodes{i}.succesors{j}}.f = f_new;
							end
						end
					end
				end

				difference_val = 0;
				for i = 1:length(obj.nodes)
					difference_val = difference_val + abs(obj.nodes{i}.f - old_nodes{i}.f);
				end

				if difference_val == 0
					disp('converged')
					iter
					% pause()
					break
				end
			end

			if obj.nodes{start_node}.connects == 0
				error('no connection');
				explored = {start_node};
				pause()
				return
			end

			explored = {};
			current_node = start_node;
			q_current = q_start';
			q = [q_start'];
			p = [];
			f = [];
			f_ext = [];
			f_v = [];
			R = [];
			L = [];
			lam = [];
			t_last = {};
			prev_t = 0;
			t = 0;

			% fins the shortest mechanically feasible path
			while true
				fprintf('\n\n\n\n\n\n STARTED IT \n\n\n\n\n\n')
				t
				if size(find([explored{:}] == current_node),2) <= 0
					explored{end+1} = current_node;
					t_last{end+1} = prev_t;
				end

				f_min = inf;
				next_node = 0;
				q_ref = q;
				if current_node == goal_node
					obj.add_cost = true;
					[q_, p_, f_, f_ext_, f_v_, L_, R_, flag] = solveTrajAndContactsBacktrack(obj, q, current_node, current_node, q_goal', 20);
					
					if flag
						next_node = current_node;
						q_prev = q_;
						p_prev = p_;
						f_prev = f_;
						L_prev = L_;
						R_prev = R_;
						f_v_prev = f_v_;
						f_ext_prev = f_ext_;
					end
					obj.add_cost = false;
					next_node = -1;
				else
					succesors_order = [];
					for i = 1:length(obj.nodes{current_node}.succesors)
						if obj.nodes{obj.nodes{current_node}.succesors{i}}.f < f_min
							if size(find([explored{:}] == obj.nodes{current_node}.succesors{i}),2) == 0		
								succesors_order = [obj.nodes{current_node}.succesors{i}, succesors_order];
								f_min = obj.nodes{obj.nodes{current_node}.succesors{i}}.f;
							end
						end
					end
					% succesors_order
					% pause()
					flag = 0;
					for i = 1:size(succesors_order,2)
						k = succesors_order(i);
						if flag == 0
							[q_, p_, f_, f_ext_, f_v_, L_, R_, flag] = solveTrajAndContactsBacktrack(obj, q, current_node, k);
							
							if flag
								% pause()
								f_min = obj.nodes{k}.f;
								next_node = k;
								q_prev = q_;
								p_prev = p_;
								f_prev = f_;
								L_prev = L_;
								R_prev = R_;
								f_v_prev = f_v_;
								f_ext_prev = f_ext_;
								diff_prev = size(q_prev,2) - size(q_ref,2);
							end
						end
					end
				end

				if next_node == 0
					next_node = explored{end-1};
					explored(end) = [];
					obj.nodes{current_node}.f = inf;

					q = q(:,1:end-t_last{end});
					p = p(:,:,:,1:end-t_last{end});
					f = f(:,:,:,1:end-t_last{end});
					f_ext = f_ext(:,:,1:end-t_last{end});
					f_v = f_v(:,:,1:end-t_last{end});

					L = L(:,:,:,1:end-t_last{end});
					R = R(:,:,:,1:end-t_last{end});

					t_last(end) = [];

					if length(explored) == 0
						t = 0;
					end
					t = t - 1;
				else
					if t == 0
						prev_t = size(p_prev,4);
					else 
						prev_t = diff_prev;
					end

					q = q_prev;
					p = p_prev;
					f = f_prev;
					L = L_prev;
					R = R_prev;
					f_ext = f_ext_prev;
					f_v = f_v_prev;

					t = t + 1;
				end

				if current_node == goal_node
					break;
				end

				current_node = next_node;
				fprintf('\n\n\n\n\n\n DID IT \n\n\n\n\n\n')
			end
		end

		function obj = searchRegionsAstar(obj,start_node, goal_node)
			% creates lists
			open_list = {start_node};
			closed_list = {};

			% performs A* search()
			searching = true;
			while searching
				% checks there are unexplored nodes
				if length(open_list) == 0; break; end;
				new_open_list = open_list;
				for i = length(open_list):1
					node = obj.nodes{open_list{i}};

					% pops last itme
					new_open_list = {new_open_list{1:i-1}};

					% iterates over succesors
					for j = 1:length(node.succesors)
						if node.succesors{j} == goal_node
							searching = false;
							break
						end

						g = obj.metric(obj.nodes{node.succesors{j}},node);
						h = obj.metric(obj.nodes{node.succesors{j}},obj.nodes{goal_node});
						f = g + h;

						obj.nodes{node.succesors{j}}.f = f;
						obj.nodes{node.succesors{j}}.prnt = open_list{i};
					end
					
					for j = 1:length(node.succesors)
						works = true;
						% checks other nodes of same level are all of higher f
						for k = 1:length(open_list)
							if node.succesors{j} == open_list{k}
							else
								if obj.nodes{open_list{k}}.f < obj.nodes{node.succesors{j}}.f
									if obj.nodes{open_list{k}}.prnt == obj.nodes{node.succesors{j}}.prnt
										works = false;
									end
								end
							end
						end

						for k = 1:length(closed_list)
							% checks other nodes of same level are all of higher f
							if node.succesors{j} == closed_list{k}
							else
								if obj.nodes{closed_list{k}}.f < obj.nodes{node.succesors{j}}.f
									if obj.nodes{closed_list{k}}.prnt == obj.nodes{node.succesors{j}}.prnt
										works = false;
									end
								end
							end
						end

						if works
							new_open_list = [node.succesors{j}, new_open_list];
						end
					end

					closed_list{end+1} = open_list{i};
				end
				open_list = new_open_list;
			end
		end

		function d = metric(obj,node1,node2)
			d_x  = node1.center(1) - node2.center(1);
			d_y  = node1.center(2) - node2.center(2);
			d_th = node1.center(3) - node2.center(3);

			d = abs(d_x)/2 + abs(d_y)/2 + abs(d_x)/pi;
		end

		function q = solveTraj(obj, q_start, q_goal, plan, plot_video)
			if nargin < 5
				plot_video = false;
			end

			if plot_video
				figure(770)
				for i = 1:length(plan)
					if i == 1
						q = [q_start'];
					else
						q = [q, obj.nodes{plan{i}}.center];
					end 

					if i < length(plan)
						p1 = polyshape(obj.nodes{plan{i}}.v');
						p2 = polyshape(obj.nodes{plan{i+1}}.v');

						p_i = intersect(p1,p2);
						if p_i.NumRegions == 0
							for k = 1:size(p1.Vertices,1)
								for j = 1:size(p2.Vertices,1)
									if p1.Vertices(k,1) == p2.Vertices(j,1)
										if p1.Vertices(k,2) == p2.Vertices(j,2)
											q_e = [p1.Vertices(k,:)'; obj.nodes{plan{i}}.center(3)];
											q = [q, q_e];
										end
									end
								end
							end
						else
							th_avg = obj.nodes{plan{i+1}}.center(3) + obj.nodes{plan{i+1}}.center(3);
							q_e = [p_i.Vertices(1,:)'; th_avg/2];
							q = [q, q_e];
						end
					end

					plot(polyshape(obj.nodes{plan{i}}.v'),'FaceColor','blue','FaceAlpha',i*0.99/length(plan));
					hold on;
				end
				q = [q, q_goal'];

				for i = 1:length(obj.environment)
					plot(obj.environment{i}(1,:),obj.environment{i}(2,:),'k','LineWidth',1)
					hold on
				end
				ylim([-0.5,1])
				xlim([-1,1])
				box on;
				hold off;

				pause()

				time = linspace(0,1,size(q,2));
				t1 = linspace(0,1,2*size(q,2));
				q = [interp1(time,q(1,:),t1); interp1(time,q(2,:),t1); interp1(time,q(3,:),t1)];

				for t = 1:size(q,2)
					h = figure(711);
					clf(h)

					% plots the environment
					for i = 1:length(obj.environment)
						plot(obj.environment{i}(1,:),obj.environment{i}(2,:),'k','LineWidth',1)
						hold on
					end

					th = q(3,t);
					rot = [cos(th),-sin(th);sin(th),cos(th)];
					object_rot = {};

					for i = 1:length(obj.object)
						object_rot{end+1} = struct('v',q(1:2,t) - rot*obj.object{i}.v);
					end
					for i = 1:length(object_rot)
						plot(polyshape(object_rot{i}.v'),'FaceColor','red','FaceAlpha',0.1);
						hold on	
					end
					ylim([-0.5,1])
					xlim([-1,1])
					box on;
					hold off;
					pause(0.1)
				end

				pause()

				close all;
			else
				for i = 1:length(plan)-1
					if i == 1
						q = [q_start'];
					else
						q = [q, obj.nodes{plan{i}}.center];
					end 

					if i < length(plan)
						p1 = polyshape(obj.nodes{plan{i}}.v');
						p2 = polyshape(obj.nodes{plan{i+1}}.v');

						p_i = intersect(p1,p2);
						if p_i.NumRegions == 0
							for k = 1:size(p1.Vertices,1)
								for j = 1:size(p2.Vertices,1)
									if p1.Vertices(k,1) == p2.Vertices(j,1)
										if p1.Vertices(k,2) == p2.Vertices(j,2)
											q_e = [p1.Vertices(k,:)'; obj.nodes{plan{i}}.center(3)];
											q = [q, q_e];
										end
									end
								end
							end
						else
							q_e = [p_i.Vertices(1,:)'; obj.nodes{plan{i}}.center(3)];
							q = [q, q_e];
						end
					end
				end
				q = [q, q_goal'];
			end
		end

		function q_cand = solveLocalTraj(obj, q_start, p_start, p1, p2, p_i, i, plan, trial)
			if p_i.NumRegions == 0
				for k = 1:size(p1.Vertices,1)
					for j = 1:size(p2.Vertices,1)
						if p1.Vertices(k,1) == p2.Vertices(j,1)
							if p1.Vertices(k,2) == p2.Vertices(j,2)
								q_e = [p1.Vertices(k,:)'; obj.nodes{plan{i}}.center(3)];
							end
						end
					end
				end
			else
				v = zeros(1,2);
				sv = 0;
				for kv = 1:size(p_i.Vertices,1)
					ran = rand(1);
					v = v + ran*p_i.Vertices(kv,:);
					sv = sv + ran;
				end
				q_e = [v'/sv; obj.nodes{plan{i}}.center(3)];
			end

			if i > 1
				q_mid = [];
				for tr = 1:trial
					v = zeros(1,2);
					sv = 0;
					for kv = 1:size(p1.Vertices,1)
						ran = rand(1);
						v = v + ran*p1.Vertices(kv,:);
						sv = sv + ran;
					end
					q_middle = [v'/sv; obj.nodes{plan{i}}.center(3)];
					q_mid = [q_mid, q_middle];
				end
				q_cand = [q_mid, q_e];
			else
				q_cand = [q_e];
			end
		end

		function q_cand = solveLocalTrajNode(obj, q_start, p_start, p1, p2, p_i, node, trial)
			
			% check if you are in a facet
			% in_facet = -1;
			% for i = 1:size(p1.Vertices,1)
			% 	ip1 = i+1;
			% 	if i == size(p1.Vertices,1)
			% 		ip1 = 1;
			% 	end
			% 	dq = q_start(1:2) - p1.Vertices(i,:)';
			% 	dv = p1.Vertices(ip1,:)' - p1.Vertices(i,:)';
			% 	if q_start(1) <= min(p1.Vertices(ip1,1),p1.Vertices(i,1)) && q_start(1) >= max(p1.Vertices(ip1,1),p1.Vertices(i,1))
			% 		if q_start(2) <= min(p1.Vertices(ip1,2),p1.Vertices(i,2)) && q_start(2) >= max(p1.Vertices(ip1,2),p1.Vertices(i,2))
			% 			if dq(2)*dv(1) - dv(2)*dq(1) == 0
			% 				in_facet = i;
			% 			end
			% 		end
			% 	end
			% end

			if p_i.NumRegions == 0
				for k = 1:size(p1.Vertices,1)
					for j = 1:size(p2.Vertices,1)
						if p1.Vertices(k,1) == p2.Vertices(j,1)
							if p1.Vertices(k,2) == p2.Vertices(j,2)
								q_e = [p1.Vertices(k,:)'; obj.nodes{node}.center(3)];
							end
						end
					end
				end
			else
				v = zeros(1,2);
				sv = 0;
				for kv = 1:size(p_i.Vertices,1)
					ran = rand(1);
					v = v + ran*p_i.Vertices(kv,:);
					sv = sv + ran;
				end
				if trial > 1 && rand(1) > 0.5
					q_e = [v'/sv; obj.nodes{node}.center(3)];
				else
					q_e = [p_poly_dist(q_start(1),q_start(2),p_i.Vertices(:,1),p_i.Vertices(:,2)); obj.nodes{node}.center(3)];
				end
			end

			q_mid = [];
			if trial < 2
				q_mid = (q_start + q_e)/2;
			else
				for tr = 1:1
					if rand(1) > 0.5
						v = zeros(1,2);
						sv = 0;
						for kv = 1:size(p1.Vertices,1)
							ran = rand(1);
							v = v + ran*p1.Vertices(kv,:);
							sv = sv + ran;
						end
						q_middle = [v'/sv; obj.nodes{node}.center(3)];
					elseif rand(1) > 0.25
						vi = randi([1,size(p1.Vertices,1)]);
						ran = rand(1);

						vip1 = vi + 1;
						if vi == size(p1.Vertices,1); vip1 = 1; end;

						v = p1.Vertices(vi,:)*ran + p1.Vertices(vip1,:)*(1-ran)
						q_middle = [v'; obj.nodes{node}.center(3)];
					else
						vi = randi([1,size(p1.Vertices,1)]);
						
						v = p1.Vertices(vi,:);
						q_middle = [v'; obj.nodes{node}.center(3)];
					end

					q_mid = [q_mid, q_middle];
				end
			end
			q_cand = [q_mid, q_e];
		end

		function [q, p, f, f_ext, f_v, L, R, flag] = solveTrajAndContactsBacktrack(obj, q_start, n1, n2, q_goal, MAX_TRIALS)

			if nargin > 4
				final = true;
			else
			 	final = false; 
			 	MAX_TRIALS = 20;
			end

			% MAX_TRIALS = 200;

			prev_q = q_start;
			q = [prev_q];
			i = 1;
			p = [];
			f = [];
			f_ext = [];
			f_v = [];
			iters = 1;
			trial = 1;

			solving = true;
			while solving
				try
					v1 = obj.nodes{n1}.v;
					v2 = obj.nodes{n2}.v;

					p1 = polyshape(obj.nodes{n1}.v');
					p2 = polyshape(obj.nodes{n2}.v');

					p_i = intersect(p1,p2);
					trial
					q_cand = obj.solveLocalTrajNode(prev_q(:,end), p(:,:,:,end), p1, p2, p_i, n1, trial);
					
					if final
						if trial == 1
							prev_q
							q_cand(:,end-1) = (prev_q(:,end) + q_goal)/2;
						end
						q_cand(:,end) = q_goal;
					end

					q = [prev_q, q_cand];
					q

					[p,f,f_ext,f_v,L,R,lam,dp,ddp] = obj.solveMIQP(q);
					flag = true;
					solving = false;
				catch
					trial = trial + 1;
					if trial > MAX_TRIALS
						flag = false;
						p = [];
						q = [];
						f = [];
						f_ext = [];
						f_v = [];
						solving = false;
					end
				end
			end
		end

		function [q, p, f, f_ext, f_v, L, R, lam, dp, ddp, flag] = solveTrajAndContactsNode(obj, q_start, n1, n2, L, R, lam, dp, ddp, q_goal)
			if nargin > 4
				prev_sol = true;
			else
			 	prev_sol = false; 
			end

			if nargin > 9
				final = true;
			else
			 	final = false; 
			end

			prev_q = q_start';
			q = [prev_q];
			i = 1;
			p = [];
			f = [];
			f_ext = [];
			iters = 1;
			trial = 1;

			MAX_TRIALS = 3;
			solving = true;
			while solving
				try
					v1 = obj.nodes{n1}.v;
					v2 = obj.nodes{n2}.v;

					p1 = polyshape(obj.nodes{n1}.v');
					p2 = polyshape(obj.nodes{n2}.v');

					p_i = intersect(p1,p2);
					q_cand = obj.solveLocalTrajNode(prev_q, p(:,:,:,end), p1, p2, p_i, n1, min(trial,2));
					if final
						if trial == 1
							prev_q
							q_goal
							q_cand
							q_cand(:,end-1) = (prev_q + q_goal)/2;
						end
						q_cand(:,end) = q_goal;
					end

					q_e = q_cand(:,end);

					if prev_sol
						[p,f,f_ext,L,R,lam,dp,ddp] = obj.solveMIQP([prev_q, q_cand], L, R, lam,dp,ddp);
					else
						[p,f,f_ext,L,R,lam,dp,ddp] = obj.solveMIQP([prev_q, q_cand]);
					end

					q = [q, q_cand];
					flag = true;
					solving = false;
				catch
					trial = trial + 1;
					if trial > MAX_TRIALS
						flag = false;
						p = [];
						q = [];
						f = [];
						f_ext = [];
						L = [];
						R = [];
						lab = [];
						solving = false;
					end
				end
			end
		end

		function [q, p, f, f_ext, f_v, L, R, lam, dp, ddp, flag] = solveTrajAndContactsNodeWithGoal(obj, q_start, n1, n2, q_goal)
			
			prev_sol = false; 
			final = true;

			prev_q = q_start';
			q = [prev_q];
			i = 1;
			p = [];
			f = [];
			f_ext = [];
			iters = 1;
			trial = 1;

			MAX_TRIALS = 10;
			solving = true;
			while solving
				try
					v1 = obj.nodes{n1}.v;
					v2 = obj.nodes{n2}.v;

					p1 = polyshape(obj.nodes{n1}.v');
					p2 = polyshape(obj.nodes{n2}.v');

					p_i = intersect(p1,p2);
					q_cand = obj.solveLocalTrajNode(prev_q, p(:,:,:,end), p1, p2, p_i, n1, min(trial,3));
					q_cand(:,end) = q_goal;
					
					q_e = q_cand(:,end);

					[p,f,f_ext,L,R,lam,dp,ddp] = obj.solveMIQP([prev_q, q_cand]);

					q = [q, q_cand];
					flag = true;
					solving = false;
				catch
					trial = trial + 1;
					if trial > MAX_TRIALS
						flag = false;
						p = [];
						q = [];
						f = [];
						f_ext = [];
						L = [];
						R = [];
						lab = [];
						solving = false;
					end
				end
			end
		end

		function [q, p, f, f_ext] = solveTrajAndContacts(obj, q_start, q_goal, plan)
			prev_q = q_start';
			prev_q_cand = []
			q = [prev_q];
			i = 1;

			p = [];
			f = [];
			f_ext = [];
			iters = 1;
			trial = 1;
			while i < length(plan)+1
				q_cand = [];
				try

					if i == length(plan)
						q_e = q_goal';

						p1 = polyshape(obj.nodes{plan{i}}.v');
						q_mid = [];
						for tr = 1:trial
							v = zeros(1,2);
							sv = 0;
							for kv = 1:size(p1.Vertices,1)
								ran = rand(1);
								v = v + ran*p1.Vertices(kv,:);
								sv = sv + ran;
							end
							q_mid = [q_mid, [v'/sv; obj.nodes{plan{i}}.center(3)]];
						end
						q_cand = [q_cand, q_mid, q_e];
					else
						v1 = obj.nodes{plan{i}}.v;
						v2 = obj.nodes{plan{i+1}}.v;

						p1 = polyshape(obj.nodes{plan{i}}.v');
						p2 = polyshape(obj.nodes{plan{i+1}}.v');

						p_i = intersect(p1,p2);
						q_cand = obj.solveLocalTraj(prev_q, p(:,:,:,end), p1, p2, p_i, i, plan, trial);
						q_e = q_cand(:,end);
					end

					if i == 1
						[p_,f_,f_ext_,L,R,lam] = obj.solveMIQP([prev_q, q_cand]);
					else
						[p_,f_,f_ext_,L,R,lam] = obj.solveMIQP([prev_q, q_cand],L, R,lam);
					end
					trial = 1;
					iters = 1;

					if i > 1
						p = cat(4, p, p_(:,:,:,2:end));
						f = cat(4, f, f_(:,:,:,2:end));
						f_ext = cat(3, f_ext, f_ext_(:,:,2:end));
					else
						p = cat(4, p, p_);
						f = cat(4, f, f_);
						f_ext = cat(3, f_ext, f_ext_);
					end

					prev_q = q_e;
					prev_q_cand = [q_cand];
					q = [q, q_cand];
					i = i + 1;
				catch
					% if trial < 3
					% 	trial = trial + 1;
					% else
					% 	trial = 1;
					% 	if iters < 5
					% 		i = i - 1;
					% 		q = q(:,1:end-size(prev_q_cand,2));
					% 		prev_q = q(:,end);
					% 		prev_q_cand = q(:,end);
					% 		p = p(:,:,:,1:end-size(prev_q_cand,2));
					% 		f = f(:,:,:,1:end-size(prev_q_cand,2));
					% 		f_ext = f_ext(:,:,1:end-size(prev_q_cand,2));
					% 		iters = iters + 1;
					% 	else
					% 		i = 1;
					% 		iters = 1;
					% 		p = [];
					% 		f = [];
					% 		f_ext = [];
					% 		prev_q = q_start';
					% 		prev_q_cand = []
					% 		q = [prev_q];
					% 	end
					% end
				end
			end
		end

		function [p, f, f_ext, f_v, L, R, lambda, dp_prev, ddp_prev] = solveMIQP(obj, q_traj, L, R, lambda,dp,ddp)
			if nargin < 3
				pre_sol = false;
			else
				pre_sol = true;
			end

			pre_sol = false;

			if pre_sol
				task = Task(obj.object, obj.vertices, obj.environment, q_traj,dp,ddp);
			else
				task = Task(obj.object, obj.vertices, obj.environment, q_traj);
			end
			miqp = CTO(task,2,1);
			miqp.McCormick = 1;
			miqp.M = 8;

			miqp = miqp.addDynamicsConstraints();
			miqp = miqp.addContactConstraints();
			miqp = miqp.addNonPenetrationConstraints();
			
			if obj.add_cost; miqp = miqp.addCostFunction(); end;

			if obj.one_finger; miqp = miqp.addOneFingerConditions(); end;
			if obj.yumi_kin; miqp = miqp.addYumiKinematicConstraints(); end;

			if pre_sol
				miqp = miqp.addInitialConditions(L, R, lambda);
			end

			miqp = miqp.solve();

			p = miqp.vars.p.value;
			f = miqp.vars.f.value;
			f_ext = miqp.vars.f_ext.value;
			f_v = miqp.vars.f_v.value;
			L = miqp.vars.L.value;
			R = miqp.vars.R.value;
			lambda = miqp.vars.lambda.value;
			dp_prev = task.traj.dr;
			ddp_prev = task.traj.ddr;
		end
	end
end