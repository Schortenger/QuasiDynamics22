clear all; close all; clc;

x = 0.2*cos(linspace(0,2*pi,7));
y = 0.2*sin(linspace(0,2*pi,7));
object = {struct('v',[x;y])};
object_verts =       [x;y];

env = {[1,1;0.5,1],[1,-0.7;1,1],[-0.7,-0.7;1,0],[-0.7,0.5;0,0],[0.5,0.5;0,0.5],[0.5,1;0.5,0.5]};

th_g = -pi/3;
q_start = [0.297     , 0.175+1e-6,0];
q_goal =  [0.9-0.1733, 0.7-0.0267+1e-6,th_g];
% q_goal =  [0.6,0.6,pi/2];
% q_goal = [0.425-1e-6,0.6,-pi/4];

options = struct('M',6,'dt',0.1,'slices',linspace(-th_g,th_g,5));
np = NonPrehensilePlanningProblem(object,env,options);
np.one_finger = true;
np.vertices = object_verts;
np = np.buildCspaceMap();

[q, p, f, f_ext, f_v] = np.solvePlanBacktrack(q_start, q_goal)

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