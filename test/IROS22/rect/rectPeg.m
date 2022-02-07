clear all; close all; clc;

% object = {struct('v',0.115*[-1, 1, 0; -1, -1, 1])};
% object_verts = 0.115*[-1, 1, 0; -1, -1, 1];
object = {struct('v',[-0.08, 0.08, 0.08, -0.08; -0.3, -0.3, 0.3, 0.3])};
object_verts = [-0.08, 0.08, 0.08, -0.08; -0.3, -0.3, 0.3, 0.3];

% object = {struct('v',0.11*[-1, 1, 0; -1, -1, 0]),struct('v',0.11*[1, -1, 0; 1, 1, 0])};

env = {[1,1;0,1],[1,-1;1,1],[-1,-1;1,0],[-1,-0.1001;0,0],[-0.1001,-0.1001;0,-0.5],[-0.1001,0.1001;-0.5,-0.5],[0.1001,0.1001;-0.5,0],[0.1001,1;0,0]};

% q_start = [0.6,0.1,-pi/2];
q_start =  [0,-0.2+1e-10,0];
q_goal =  [0.6,0.08+1e-10,pi/2];

np = NonPrehensilePlanningProblem(object,env);
% np.one_finger = true;
np.vertices = object_verts;
np = np.buildCspaceMap();

[q, p, f, f_ext] = np.solvePlanBacktrack(q_start, q_goal)

task = Task(object, object_verts, env, q);
plan = struct();

plan.N_T = size(p,4);
plan.N_l = size(p,3);
plan.N_c = size(p,2);

plan.p = p;
plan.f = f;
plan.f_ext = f_ext;

animation_nonP(task,env,plan)