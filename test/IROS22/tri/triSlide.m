clear all; close all; clc;

% object = {struct('v',0.115*[-1, 1, 0; -1, -1, 1])};
% object_verts = 0.115*[-1, 1, 0; -1, -1, 1];
object = {struct('v',2*[0.1, -0.1, 0; 0.05, 0.05, -0.1])};
object_verts = 2*[0.1, -0.1, 0; 0.05, 0.05, -0.1];

% object = {struct('v',0.11*[-1, 1, 0; -1, -1, 0]),struct('v',0.11*[1, -1, 0; 1, 1, 0])};

env = {[1,1;0.5,1],[1,-1;1,1],[-1,-1;1,0],[-1,0.5;0,0],[0.5,0.5;0,0.5],[0.5,1;0.5,0.5]};

q_start = [0.4,0.201,pi/2];
q_goal =  [0.6,0.601,0];

options = struct('M',6,'dt',0.1,'slices',linspace(-pi/2,pi/2,5));
np = NonPrehensilePlanningProblem(object,env,options);
% np.one_finger = true;
np.vertices = object_verts;
np = np.buildCspaceMap();

plan = np.searchRegionsDP(q_start, q_goal);

% for th = linspace(-pi,0,9)

% 	object_1 = {};
% 	rot = [cos(th),-sin(th);sin(th),cos(th)];

% 	for i = 1:length(object)
% 		object_1{end+1} = struct('v',rot*object{i}.v);
% 	end

% 	figure(1)

% 	subplot(4,1,1)
% 	for i = 1:length(object)
% 		plot(polyshape(object_1{i}.v'),'FaceColor','red','FaceAlpha',0.1);
% 		hold on
% 	end
% 	ylim([-0.75,0.75])
% 	xlim([-1,1])
% 	box on;
% 	hold off;
% 	title('object at orietation')

% 	subplot(4,1,2)
% 	for i = 1:length(env)
% 		plot(env{i}(1,:),env{i}(2,:),'k','LineWidth',1)
% 		hold on
% 	end
% 	ylim([-0.5,1])
% 	xlim([-1,1])
% 	box on;
% 	hold off;
% 	title('workspace')

% 	subplot(4,1,3)

% 	for i = 1:length(env)
% 		plot(env{i}(1,:),env{i}(2,:),'k','LineWidth',1)
% 		hold on
% 	end
% 	cspace = FreeSpaceAtSlice(env,object_1);
% 	for i = 1:length(cspace)
% 		plot(polyshape(cspace{i}'))
% 		hold on
% 	end
% 	ylim([-0.5,1])
% 	xlim([-1,1])
% 	box on;
% 	hold off;
% 	title('Free-Space a C-slice')

% 	subplot(4,1,4)
% 	for i = 1:length(env)
% 		plot(env{i}(1,:),env{i}(2,:),'k','LineWidth',1)
% 		hold on
% 	end

% 	cspace = FreeSpaceAtSlice(env,object_1);
% 	for i = 1:length(cspace)
% 		decomp = getConvexDecomposition(cspace{i})
% 		cvx_dec = decomp.polys;
% 		for j = 1:length(cvx_dec)
% 			plot(polyshape(cvx_dec{j}'))
% 			hold on
% 		end
% 	end
% 	ylim([-0.5,1])
% 	xlim([-1,1])
% 	box on;
% 	hold off;
% 	title('segmented Free-Space')

% 	pause();
% end


% figure(2)
% for i = 1:length(plan)
% 	p1 = polyshape(np.nodes{plan{i}}.v');
% 	plot(p1)
% 	hold on;
% end

% pause()

[q, p, f, f_ext, f_v] = np.solvePlanBacktrack(q_start, q_goal)

% figure(2)
% plot(q(1,:),q(2,:),'*')
% hold on;
% pause()

task = Task(object, object_verts, env, q);
plan = struct();

plan.N_T = size(p,4);
plan.N_l = size(p,3);
plan.N_c = size(p,2);

plan.p = p;
plan.f = f;
plan.f_ext = f_ext;
plan.f_v = f_v;

animation_nonP(task,env,plan)