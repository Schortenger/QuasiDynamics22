clear all; close all; clc;

% finds a task
for j = 29:40
	i = 1;
	data = {}
	while true
		try 
			begn = [0.2*rand()-0.1;0.2*rand()-0.1;(rand()-0.5)*pi/3];
			goal = [0.2*rand()-0.1;0.2*rand()-0.1;(rand()-0.5)*pi/3];

			% generates traj
			task = GenTraj(begn,goal);

			% infers the primitive
			pre = MixedIntegerContactPlacementProblem(task,2,1)
			pre.McCormick = 1;
			pre.M = 4;

			% animation_shape(task,false)
			% pause()

			pre = pre.addDynamicsConstraints();
			pre = pre.addContactConstraints();
			pre = pre.addNonPenetrationConstraints();
			pre = pre.addCostFunction();

			pre = pre.solve();

			% animates and saves the video
			% animation_contacts(task,pre,false)
			% pause()
			data{i}.pre = pre;

			x = [];
			y = [];
			for v = 1:task.nv
				new_vert = task.v(:,v);

				x = [x, new_vert(1)];
				y = [y, new_vert(2)];
			end

			f = figure('visible','off');
			% figure(1)
			pgon = polyshape(x,y);
			plot(pgon,'FaceAlpha',0,'FaceColor','white','EdgeColor','black')
			box on
			xlim([-0.05,0.05])
			ylim([-0.05,0.05])

			F = getframe(gcf);
			[X, ~] = frame2im(F);
			X = imresize(X,[50 50]);
			X = rgb2gray(X);
			
			data{i}.img = X;
		catch
			i = i - 1;
		end
		i = i + 1;
		if i > 200; break; end
	end

	ins = [];
	outs = [];

	for i = 1:length(data)
		img = reshape(rescale(double(data{i}.img),0.0,1.0),[1,50*50]);
		ins1 = [reshape(data{i}.pre.object.traj.r,[1,15]),...
				reshape(data{i}.pre.object.traj.dr,[1,15]),...
				reshape(data{i}.pre.object.traj.ddr,[1,15]),...
				img];

		outs1 = reshape(data{i}.pre.vars.p.value,[1,2*2*5]);

		ins = [ins; ins1];
		outs = [outs; outs1];
	end

	res = [ins,outs];

	csvwrite(strcat(strcat('./data/data_',num2str(j)),'_2f_pushing_randomized.csv'),res);
	save(strcat(strcat('./data/raw_',num2str(j)),'_2f_pushing_randomized.mat'),'data');
end