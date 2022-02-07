clear all; clc; close all;

% object_ref = {struct('v',0.25*[-1, 1, 0; -1, -1, 1])};
% object_ref = {struct('v',0.15*[-1, 1, 0; -1, -1, 0]),struct('v',0.15*[1, -1, 0; 1, 1, 0])};
% env = {[1,1;0,1],[1,0.25;1,1],[0.25,0;1,0.75],[0,-0.25;0.75,1],[-0.25,-1;1,1],[-1,-1;1,0],[-1,-0.25;0,0],[-0.25,0;0,0.25],[0,0.25;0.25,0],[0.25,1;0,0]};
% env = {[1,1;0,1],[1,-1;1,1],[-1,-1;1,0],[-1,1;0,0]};

object_ref = {struct('v',[-0.1, 0.1, 0.1, -0.1; -0.3, -0.3, 0.3, 0.3])};
env = {[1,1;0,1],[1,-1;1,1],[-1,-1;1,0],[-1,-0.5;0,0],[-0.5,-0.5;0,-0.5],[-0.5,0.5;-0.5,-0.5],[0.5,0.5;-0.5,0],[0.5,1;0,0]};

% animation
rang = -pi/4:pi/16:pi/4;
for th = rang

	object = {};
	rot = [cos(th),-sin(th);sin(th),cos(th)];

	for i = 1:length(object_ref)
		object{end+1} = struct('v',rot*object_ref{i}.v);
	end

	figure(1)

	subplot(4,1,1)
	for i = 1:length(object_ref)
		plot(polyshape(object{i}.v'),'FaceColor','red','FaceAlpha',0.1);
		hold on
	end
	ylim([-0.75,0.75])
	xlim([-1,1])
	box on;
	hold off;
	title('object at orietation')

	subplot(4,1,2)
	for i = 1:length(env)
		plot(env{i}(1,:),env{i}(2,:),'k','LineWidth',1)
		hold on
	end
	ylim([-0.5,1])
	xlim([-1,1])
	box on;
	hold off;
	title('workspace')

	subplot(4,1,3)

	for i = 1:length(env)
		plot(env{i}(1,:),env{i}(2,:),'k','LineWidth',1)
		hold on
	end
	cspace = FreeSpaceAtSlice(env,object);
	for i = 1:length(cspace)
		plot(polyshape(cspace{i}'))
		hold on
	end
	ylim([-0.5,1])
	xlim([-1,1])
	box on;
	hold off;
	title('Free-Space a C-slice')

	subplot(4,1,4)
	for i = 1:length(env)
		plot(env{i}(1,:),env{i}(2,:),'k','LineWidth',1)
		hold on
	end

	cspace = FreeSpaceAtSlice(env,object);
	for i = 1:length(cspace)
		decomp = getConvexDecomposition(cspace{i})
		cvx_dec = decomp.polys;
		for j = 1:length(cvx_dec)
			plot(polyshape(cvx_dec{j}'))
			hold on
		end
	end
	ylim([-0.5,1])
	xlim([-1,1])
	box on;
	hold off;
	title('segmented Free-Space')

	pause(0.01);
end

% c-slices
figure(2)
k = 0;
for th = rang
	k = k + 1;
	object = {};
	subplot(size(rang,2),1,k)
	
	rot = [cos(th),-sin(th);sin(th),cos(th)];
	for i = 1:length(object_ref)
		object{end+1} = struct('v',rot*object_ref{i}.v);
	end

	for i = 1:length(env)
		plot(env{i}(1,:),env{i}(2,:),'k','LineWidth',1)
		hold on
	end

	cspace = FreeSpaceAtSlice(env,object);
	for i = 1:length(cspace)
		decomp = getConvexDecomposition(cspace{i})
		cvx_dec = decomp.polys;
		for j = 1:length(cvx_dec)
			plot(polyshape(cvx_dec{j}'))
			hold on
		end
	end
	ylim([-0.5,1])
	xlim([-1,1])
	box on;
	hold off;

	pause(0.01);
end

figure(3)
k = 0;
for th = rang
	k = k + 1;
	object = {};
	subplot(size(rang,2),1,k)
	
	rot = [cos(th),-sin(th);sin(th),cos(th)];

	for i = 1:length(object_ref)
		object{end+1} = struct('v',rot*object_ref{i}.v + [0.5;0.5]);
	end

	for i = 1:length(object)
		plot(polyshape(object{i}.v'),'FaceColor','red','FaceAlpha',0.1);
		hold on
	end

	for i = 1:length(env)
		plot(env{i}(1,:),env{i}(2,:),'k','LineWidth',1)
		hold on
	end

	% cspace = FreeSpaceAtSlice(env,object);
	% for i = 1:length(cspace)
	% 	decomp = getConvexDecomposition(cspace{i})
	% 	cvx_dec = decomp.polys;
	% 	for j = 1:length(cvx_dec)
	% 		plot(polyshape(cvx_dec{j}'))
	% 		hold on
	% 	end
	% end
	ylim([-0.5,1])
	xlim([-1,1])
	box on;
	hold off;
	title('')
	pause(0.01);
end