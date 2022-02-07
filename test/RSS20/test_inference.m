clear all; close all; clc;

% finds a task
task = square_sliding();

% infers the primitive
pre = MixedIntegerContactPlacementProblem(task,2,1)
pre.McCormick = 1;
pre.M = 12;

% animation_shape(task,true)
% pause()

pre = pre.addDynamicsConstraints();
pre = pre.addContactConstraints();
pre = pre.addNonPenetrationConstraints();
pre = pre.addCostFunction();

c = tic()
pre = pre.solve();
toc(c)

% animates and saves the video
animation_contacts(task,pre,false,false)

% SCS solver for MIQP