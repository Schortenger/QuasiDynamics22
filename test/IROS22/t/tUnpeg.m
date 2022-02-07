clear all; close all; clc;

object = {struct('v',1.6*[-0.05,0.05,0.05,-0.05;-0.3,-0.3,0.0,0.0]),struct('v',1.6*[-0.2,0.2,0.2,-0.2;0.0,0.0,0.1,0.1])};
object_verts = 1.6*[-0.05,0.05,0.05,0.2,0.2,-0.2,-0.2,-0.05; -0.3, -0.3, 0.0, 0.0, 0.1, 0.1, 0.0, 0.0];

env = {[1,1;0,1],[1,-1;1,1],[-1,-1;1,0],[-1,-0.1101;0,0],[-0.1101,-0.1101;0,-0.5],[-0.1101,0.1101;-0.5,-0.5],[0.1101,0.1101;-0.5,0],[0.1101,1;0,0]};

% q_start = [0.6,0.1,-pi/2];
q_start = [-0.08,0.08+1e-10,pi/2];
q_goal =  [ 0.5,0.16+1e-10,0];

np = NonPrehensilePlanningProblem(object,env);
% np.one_finger = true;
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