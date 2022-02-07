clear all; close all; clc;

% object = {struct('v',0.115*[-1, 1, 0; -1, -1, 1])};
% object_verts = 0.115*[-1, 1, 0; -1, -1, 1];
object = {struct('v',2*[-0.1, 0.1, 0.1, -0.1; -0.1, -0.1, 0.1, 0.1])};
object_verts = 2*[-0.1, 0.1, 0.1, -0.1; -0.1, -0.1, 0.1, 0.1];

% object = {struct('v',0.11*[-1, 1, 0; -1, -1, 0]),struct('v',0.11*[1, -1, 0; 1, 1, 0])};

env = {[1,1;0.5,1],[1,-1;1,1],[-1,-1;1,0],[-1,0.5;0,0],[0.5,0.5;0,0.5],[0.5,1;0.5,0.5]};

q_start = [0,0.201,0];
q_goal =  [-0.4,0.201,pi/2];

np = NonPrehensilePlanningProblem(object,env);
% np.one_finger = true;
np.vertices = object_verts;
np = np.buildCspaceMap();

[q, p, f, f_ext, f_v, L, R, explored] = np.solvePlanBacktrack(q_start, q_goal)

task = Task(object, object_verts, env, q);
plan = struct();

plan.N_T = size(p,4);
plan.N_l = size(p,3);
plan.N_c = size(p,2);

plan.p = p;
plan.L = L;
plan.R = R;
plan.f = f;
plan.f_ext = f_ext;
plan.f_v = f_v;
explored;

% nlp = CTORefinementNLP(plan, task);
% nlp = nlp.addDecisionVariables();
% nlp = nlp.addDynamicsConstraints();
% nlp = nlp.addNonPenetrationConstraints();

% nlp = nlp.SolveOpti()

animation_nonP(task,env,plan)