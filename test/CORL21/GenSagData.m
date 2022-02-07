clear all; close all; clc;

x_arr = {};
y_arr = {};

% v2 = [[cos(-pi/4:pi/16:pi/4), 0, cos(3*pi/4:pi/16:5*pi/4)]; [sin(-pi/4:pi/16:pi/4), 1.5, sin(3*pi/4:pi/16:5*pi/4)]]/20;

% x_arr{end+1} = v2(1,:);
% y_arr{end+1} = v2(2,:);
% x_arr{end+1} = v2(1,:);
% y_arr{end+1} = v2(2,:);
% x_arr{end+1} = v2(1,:);
% y_arr{end+1} = v2(2,:);

% x_arr{end+1} = [-1, 1, 1, 0, -1]/20;
% y_arr{end+1} = [-1, -1, 1, 2, 1]/20;
% x_arr{end+1} = [-1, 1, 1, 0, -1]/20;
% y_arr{end+1} = [-1, -1, 1, 2, 1]/20;
% x_arr{end+1} = [-1, 1, 1, 0, -1]/20;
% y_arr{end+1} = [-1, -1, 1, 2, 1]/20;

% v2 = [cos(pi/3), -sin(pi/3); sin(pi/3), cos(pi/3)]*[[cos(0:pi/16:pi), -1, 1]; [sin(0:pi/16:pi), -0.5, -0.5]]/20;

% x_arr{end+1} = v2(1,:);
% y_arr{end+1} = v2(2,:);
% x_arr{end+1} = v2(1,:);
% y_arr{end+1} = v2(2,:);
% x_arr{end+1} = v2(1,:);
% y_arr{end+1} = v2(2,:);

% v2 = [[cos(0:pi/16:(2*pi-pi/16))]; [sin(0:pi/16:(2*pi-pi/16))]]/20;

% x_arr{end+1} = v2(1,:);
% y_arr{end+1} = v2(2,:);
% x_arr{end+1} = v2(1,:);
% y_arr{end+1} = v2(2,:);
% x_arr{end+1} = v2(1,:);
% y_arr{end+1} = v2(2,:);

n_samples = 200;
idx_1 = 1;
n_sets = idx_1;

x_init = linspace(-0.1,0.1,n_samples);

% finds a task
for set_j = idx_1:n_sets
	i = 1;
	data = {};
	while true
		try
			begn = [0.2*(rand()-0.5);0.2*(rand()-0.5);0];
			% begn = [0.1;0.2*(rand()-0.5);0];
			goal = [0.2*(rand()-0.5);0.2*(rand()-0.5);pi/4*(rand())];

			% generates traj
			task = GenSagMotionNCVX(begn,goal,6);
			% task = GenSagMotionBottle(begn,goal,8);

			% x = x_arr{i};
			% y = y_arr{i};
			% task = GenSagMotionOmni(begn,goal,x_arr{i},y_arr{i})
			% task = GenSagMotionBlock(begn,goal,4);

			% pause()
			% infers the primitive
			pre = MixedIntegerContactPlacementProblem(task,2,1);
			pre.McCormick = 1;
			pre.M = 12;

			% animation_shape(task,false)
			% pause()

			pre = pre.addDynamicsConstraints();
			pre = pre.addContactConstraints();
			pre = pre.addNonPenetrationConstraints();
			pre = pre.addCostFunction();
			pre = pre.solve();

			% animates and saves the video
			% animation_contacts(task,pre,true,false,false)
			% pause()

			data{i}.pre = pre;

			data{i}.p = zeros(4, pre.N_T);
			data{i}.f = zeros(4, pre.N_T);

			% ground truth data
			for t = 1:pre.N_T
				data{i}.p([1,3],t) = pre.vars.p.value(:,1,1,t);
				data{i}.p([2,4],t) = pre.vars.p.value(:,1,2,t);

				data{i}.f([1,3],t) = pre.vars.f.value(:,1,1,t);
				data{i}.f([2,4],t) = pre.vars.f.value(:,1,2,t);
			end


			x = [];
			y = [];
			for v = 1:task.nv
				new_vert = task.v(:,v);

				x = [x, new_vert(1)];
				y = [y, new_vert(2)];
			end

			% store vertex and 	friction cone data
			data{i}.polygons = task.polygons;
			data{i}.v = zeros(8,pre.N_T) - 0.1;
			data{i}.fc = zeros(8,pre.N_T);
			data{i}.p_r = zeros(4, pre.N_T);
			
			for t = 1:pre.N_T

				data{i}.p_r([1,3],t) = pre.vars.p.value(:,1,1,t);
				data{i}.p_r([2,4],t) = pre.vars.p.value(:,1,2,t);
				for c = 1:pre.N_c
					for f = 1:pre.N_f
						if pre.vars.L.value(f,1,1,t) == 1
							% vertex data
							% data{i}.v([1,2],t) = pre.object.lines{f}.v1(:,t);
							% data{i}.v([3,4],t) = pre.object.lines{f}.v2(:,t);

							data{i}.v([1,2],t) = pre.vars.p.value(:,1,1,t);
							data{i}.v([3,4],t) = pre.vars.p.value(:,1,1,t);

							data{i}.fc([1,2],t) = pre.object.lines{f}.fc1(:,t);
							data{i}.fc([3,4],t) = pre.object.lines{f}.fc2(:,t);
						end

						if pre.vars.L.value(f,1,2,t) == 1
							% vertex data
							% data{i}.v([5,6],t) = pre.object.lines{f}.v1(:,t);
							% data{i}.v([7,8],t) = pre.object.lines{f}.v2(:,t);

							data{i}.v([5,6],t) = pre.vars.p.value(:,1,2,t);
							data{i}.v([7,8],t) = pre.vars.p.value(:,1,2,t);

							data{i}.fc([5,6],t) = pre.object.lines{f}.fc1(:,t);
							data{i}.fc([7,8],t) = pre.object.lines{f}.fc2(:,t);
						end
					end
					% no contact scenario
					if sum(pre.vars.L.value(:,1,1,t)) == 0
						data{i}.v([1,2],t) = pre.vars.p.value(:,1,1,t);
						data{i}.v([3,4],t) = pre.vars.p.value(:,1,1,t);
					end
					if sum(pre.vars.L.value(:,1,2,t)) == 0
						data{i}.v([5,6],t) = pre.vars.p.value(:,1,2,t);
						data{i}.v([7,8],t) = pre.vars.p.value(:,1,2,t);
					end
				end
			end

			% store external cpntact and its friction cone data
			data{i}.p_e = zeros(4,pre.N_T);
			data{i}.fc_e = zeros(8,pre.N_T);
			traj = pre.object.traj; 

			for t = 1:pre.N_T
				v_ex = task.ext_contact(1,t);
				if v_ex > 0
					rot_mat = [cos(traj.r(3,t)),-sin(traj.r(3,t));sin(traj.r(3,t)),cos(traj.r(3,t))];
						
					data{i}.p_e([1,2],t) = traj.r(1:2,t) + rot_mat*(pre.object.v(:,v_ex));
					data{i}.p_e([3,4],t) = traj.r(1:2,t) + rot_mat*(pre.object.v(:,v_ex));

					data{i}.fc_e([1,2],t) = pre.object.ext_f{v_ex}.fc1(:,t);
					data{i}.fc_e([3,4],t) = pre.object.ext_f{v_ex}.fc2(:,t);

					data{i}.fc_e([5,6],t) = pre.object.ext_f{v_ex}.fc1(:,t);
					data{i}.fc_e([7,8],t) = pre.object.ext_f{v_ex}.fc2(:,t);
				end
			end


			% Object mask
			data{i}.img = {};
			w  = 1;
			for th = linspace(-pi,pi,7)
				rot = [cos(th), sin(th); -sin(th), cos(th)];
				f = figure('visible','off');
				p = rot*[x;y];
				pgon = polyshape(p(1,:),p(2,:));

				plot(pgon,'FaceAlpha',1.0,'FaceColor','black','EdgeColor','black')
				axis off
				set(gcf,'color','w')
				xlim([-0.05,0.05])
				ylim([-0.05,0.05])

				F = getframe(gcf);
				close all;

				[Xor, ~] = frame2im(F);
				data{i}.img{w} = rgb2gray(imresize(Xor,[50 50]));
				w = w + 1;
			end

			% Environemnt mask
			
			f = figure('visible','off');

			for m = 1:length(task.env)
				x = task.env{m}.x;
				y = task.env{m}.y;
				pgon = polyshape(x,y);
				plot(pgon,'FaceAlpha',1.0,'FaceColor','black','EdgeColor','black');
				hold on;
			end

			pgon = polyshape(x,y);

			axis off;
			set(gcf,'color','w');
			xlim([-0.2,0.2]);
			ylim([-0.2,0.2]);

			F = getframe(gcf);
			close all;

			[Xor, ~] = frame2im(F);
			data{i}.env = rgb2gray(imresize(Xor,[50 50]));
			data{i}.video = generate_video_array(task);

			close all;


			f = figure('visible','off');
			pgon = polyshape(x,y);
			plot(pgon,'FaceAlpha',0.0,'FaceColor','black','EdgeColor','black')
			axis off;
			set(gcf,'color','w');
			xlim([-0.05,0.05]);
			ylim([-0.05,0.05]);

			F = getframe(gcf);
			close all;

			[Xor, ~] = frame2im(F);

			% get signed distance from map and inverse map
			field = double(bwdist(rgb2gray(Xor) > 100) - bwdist(rgb2gray(Xor) < 100));
			if isinf(field(1,1))
			    field = ones(size(field)) * 1000;
			end
			arr = imresize(field,[50 50]);
			data{i}.sdf = rescale(-reshape(arr',[1,50*50]),0.0,1.0);
		catch
			i = i - 1;
		end

		i = i + 1
		if i > n_samples; break; end
	end

	ins = [];
	outs = [];
	vids = [];
	pols = [];
	sag_params = [];

	% data to save:
	% inputs:
	% r = x_traj[i,0:15]
	% p_r = x[i,0:20]
	% fc = x[i,20:60]
	% p_e = x[i,60:80]
	% fc_e = x[i,80:120]
	% v = x[i,120:160]
	% outputs:
	% p = x[i,0:20]
	% f = x[i,20:40]

	% computes the number of triangulations
	nt = -1
	for i = 1:length(data)
		nt0 = length(data{i}.polygons);
		if  nt0 > nt
			nt = nt0;
		end
	end

	pols = ones(length(data),1 + 8*nt)*(-99.99);

	% fills the polygonal data
	for i = 1:length(data)
		pols(i,1) = length(data{i}.polygons);
		for n = 1:length(data{i}.polygons)
			for k = 1:3
				pols(i,1 + (n-1)*8 + (k-1)*2 + 1) = data{i}.polygons{n}.v(1,k);
				pols(i,1 + (n-1)*8 + (k-1)*2 + 2) = data{i}.polygons{n}.v(2,k);
			end
			pols(i,1 + (n-1)*8 + 7) = data{i}.polygons{n}.center(1);
			pols(i,1 + (n-1)*8 + 8) = data{i}.polygons{n}.center(2);
		end
	end


	for i = 1:length(data)
		img = []
		for k = 1:length(data{i}.img)
			img = [img, reshape(rescale(double(data{i}.img{k}),0.0,1.0)',[1,50*50])];
		end

		ins1 = [reshape(data{i}.pre.object.traj.r',[1,15]),...
				reshape(data{i}.pre.object.traj.dr',[1,15]),...
				reshape(data{i}.pre.object.traj.ddr',[1,15]),...
				reshape(data{i}.p_r',[1,20]),...
				reshape(data{i}.fc',[1,40]),...
				reshape(data{i}.p_e',[1,20]),...
				reshape(data{i}.fc_e',[1,40]),...
				reshape(data{i}.v',[1,40]),...
				img, reshape(rescale(double(data{i}.env),0.0,1.0)',[1,50*50])];

		outs1 = [reshape(data{i}.p',[1,20]),reshape(data{i}.f',[1,20])];
		params = [data{i}.pre.object.y_floor, data{i}.pre.object.th_floor];

		ins = [ins; ins1];
		outs = [outs; outs1];
		vids = [vids; data{i}.video];
		sag_params = [sag_params; params];
	end

	size(ins)

	res = [ins,outs];
	% j = 57;
	j = set_j;
	csvwrite(strcat(strcat('../../fair-quasidyn/data/data_sag_',num2str(j)),'_2f_sq.csv'),res);
	csvwrite(strcat(strcat('../../fair-quasidyn/data/vids_sag_',num2str(j)),'_2f_sq.csv'),vids);
	csvwrite(strcat(strcat('../../fair-quasidyn/data/polygons_sag_',num2str(j)),'_2f_sq.csv'),pols);
	csvwrite(strcat(strcat('../../fair-quasidyn/data/params_sag_',num2str(j)),'_2f_sq.csv'),sag_params);
	save(strcat(strcat('../data/AffordanceData/raw_sag_',num2str(j)),'_2f_sq.mat'),'data');
end