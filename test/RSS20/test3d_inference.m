clear all; close all; clc;

% finds a task
task = cube_picking();

% infers the primitive
pre = MixedIntegerContactTrajOptProblem(task,2)
pre.McCormick = 4;
pre.M = 12;

% animation_shape(task,false)
% pause()

pre = pre.addDynamicsConstraints();
pre = pre.addContactConstraints();
pre = pre.addNonPenetrationConstraints();
% pre = pre.addCostFunction();

c = tic()
pre = pre.solve();
toc(c)

% animates and saves the video
% animation_contacts(task,pre,fals	e,false)

% SCS solver for MIQP