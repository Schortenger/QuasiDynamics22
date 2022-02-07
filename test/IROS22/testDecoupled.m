clear all; close all; clc;

% object = {struct('v',0.25*[-1, 1, 0; -1, -1, 1])};
% object_verts = 0.25*[-1, 1, 0; -1, -1, 1];
object = {struct('v',0.24*[-1, 1, 1, -1; -1, -1, 1, 1])};
object_verts = 0.24*[-1, 1, 1, -1; -1, -1, 1, 1];

% object = {struct('v',0.15*[-1, 1, 0; -1, -1, 0]),struct('v',0.15*[1, -1, 0; 1, 1, 0])};
env = {[1,1;0,1],[1,0.25;1,1],[0.25,0;1,0.75],[0,-0.25;0.75,1],[-0.25,-1;1,1],[-1,-1;1,0],[-1,-0.25;0,0],[-0.25,0;0,0.25],[0,0.25;0.25,0],[0.25,1;0,0]};
% env = {[1,1;0,1],[1,-1;1,1],[-1,-1;1,0],[-1,1;0,0]};

q_start = [-0.5,0.5,-pi/4];
q_goal =  [-0.5,0.5,0];

np = NonPrehensilePlanningProblem(object,env);
np.vertices = object_verts;
np = np.buildCspaceMap();
plan = np.searchRegionsDP(q_start,q_goal)

% q = np.solveTraj(q_start,q_goal,plan);

[q,p,f,f_ext] = np.solveTrajAndContacts(q_start,q_goal,plan);

task = Task(object, object_verts, env, q);
plan = struct();

plan.N_T = size(p,4);
plan.N_l = size(p,3);
plan.N_c = size(p,2);

plan.p = p;
plan.f = f;
plan.f_ext = f_ext;

animation_nonP(task,env,plan)