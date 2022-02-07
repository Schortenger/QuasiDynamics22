clear all; close all; clc;

set_j = 42

res0 = load(strcat(strcat('./../data/AffordanceData/raw_sag_',num2str(set_j)),'_2f_sq.mat'));

p_GT = csvread('../../fair-quasidyn/data/sim_p_GT4_sag_0_2f.csv');
p_Sim = csvread('../../fair-quasidyn/data/sim_p_full4_sag_0_2f.csv');
p_nn = csvread('../../fair-quasidyn/data/sim_p_nn4_sag_0_2f.csv');
p_cvx = csvread('../../fair-quasidyn/data/sim_p_mdr4_sag_0_2f.csv');
p_mdr = csvread('../../fair-quasidyn/data/sim_p_cvx4_sag_0_2f.csv');

r_GT = csvread('../../fair-quasidyn/data/sim_r_GT4_sag_0_2f.csv');
r_Sim = csvread('../../fair-quasidyn/data/sim_r_full4_sag_0_2f.csv');
r_nn = csvread('../../fair-quasidyn/data/sim_r_nn4_sag_0_2f.csv');
r_cvx = csvread('../../fair-quasidyn/data/sim_r_mdr4_sag_0_2f.csv');
r_mdr = csvread('../../fair-quasidyn/data/sim_r_cvx4_sag_0_2f.csv');

data = csvread(strcat(strcat('../../fair-quasidyn/data/data_sag_',num2str(set_j)),'_2f_sq.csv'));

for i = 1:20
	pre = res0.data{i}.pre;
	task = pre.object;
	pre.vars.f.value = zeros(2,1,2,5);
	

	r_ref = [task.traj.r(1,:); task.traj.r(2,:); task.traj.r(3,:)];

	i

	% plots
	% figure(1)
	% animation_sim(task,pre,1,false,false,true,true,strcat(strcat('../videos/',num2str(i)),'vid_sag_RSS'))
	% title('GT')
	animation_sim(task,pre,2,true,false,true,false,strcat(strcat('../videos/',num2str(i)),'vid_sag_GT'))
	close all;
	title('GT')

	traj_1 = [p_GT(i, 1:5); p_GT(i, 11:15)];
	traj_2 = [p_GT(i, 6:10); p_GT(i, 16:20)];
	r = [r_GT(i, 1:5); r_GT(i, 6:10); r_GT(i, 11:15)];

	r_err_miqp = norm(r_ref - [r_GT(i, 1:5); r_GT(i, 6:10); r_GT(i, 11:15)])

	p_data = zeros(2,1,2,5);
	p_data(:,1,1,:) = traj_1;
	p_data(:,1,2,:) = traj_2;

	task.traj.r = r;
	pre.vars.p.value = p_data;

	% plots
	% figure(2)
	animation_sim(task,pre,2,true,false,true,true,strcat(strcat('../videos/',num2str(i)),'vid_sag_MIQP'))
	close all;
	title('MIQP')

	traj_1 = [p_Sim(i, 1:5); p_Sim(i, 11:15)];
	traj_2 = [p_Sim(i, 6:10); p_Sim(i, 16:20)];
	r = [r_Sim(i, 1:5); r_Sim(i, 6:10); r_Sim(i, 11:15)];

	r_err_full = norm(r_ref - [r_Sim(i, 1:5); r_Sim(i, 6:10); r_Sim(i, 11:15)])

	p_data = zeros(2,1,2,5);
	p_data(:,1,1,:) = traj_1;
	p_data(:,1,2,:) = traj_2;

	task.traj.r = r;
	pre.vars.p.value = p_data;

	% plots
	% figure(3)
	animation_sim(task,pre,3,true,false,true,true,strcat(strcat('../videos/',num2str(i)),'vid_sag_DDM'))
	close all;
	title('SIM')

	traj_1 = [p_mdr(i, 1:5); p_mdr(i, 11:15)];
	traj_2 = [p_mdr(i, 6:10); p_mdr(i, 16:20)];
	r = [r_mdr(i, 1:5); r_mdr(i, 6:10); r_mdr(i, 11:15)];

	r_err_MDR = norm(r_ref - [r_mdr(i, 1:5); r_mdr(i, 6:10); r_mdr(i, 11:15)])

	p_data = zeros(2,1,2,5);
	p_data(:,1,1,:) = traj_1;
	p_data(:,1,2,:) = traj_2;

	task.traj.r = r;
	pre.vars.p.value = p_data;

	% plots
	% figure(3)
	animation_sim(task,pre,4,true,false,true,true,strcat(strcat('../videos/',num2str(i)),'vid_sag_MDR'))
	close all;
	title('MDR')

	traj_1 = [p_cvx(i, 1:5); p_cvx(i, 11:15)];
	traj_2 = [p_cvx(i, 6:10); p_cvx(i, 16:20)];
	r = [r_cvx(i, 1:5); r_cvx(i, 6:10); r_cvx(i, 11:15)];

	r_err_cvx = norm(r_ref - [r_cvx(i, 1:5); r_cvx(i, 6:10); r_cvx(i, 11:15)])

	p_data = zeros(2,1,2,5);
	p_data(:,1,1,:) = traj_1;
	p_data(:,1,2,:) = traj_2;

	task.traj.r = r;
	pre.vars.p.value = p_data;

	% plots
	% figure(3)
	animation_sim(task,pre,5,true,false,true,true,strcat(strcat('../videos/',num2str(i)),'vid_sag_CVX'))
	close all;
	title('CVX')

	traj_1 = [p_nn(i, 1:5); p_nn(i, 11:15)];
	traj_2 = [p_nn(i, 6:10); p_nn(i, 16:20)];
	r = [r_nn(i, 1:5); r_nn(i, 6:10); r_nn(i, 11:15)];

	r_err_NN = norm(r_ref - [r_nn(i, 1:5); r_nn(i, 6:10); r_nn(i, 11:15)])

	p_data = zeros(2,1,2,5);
	p_data(:,1,1,:) = traj_1;
	p_data(:,1,2,:) = traj_2;

	task.traj.r = r;
	pre.vars.p.value = p_data;

	% plots
	% figure(3)
	animation_sim(task,pre,6,true,false,true,true,strcat(strcat('../videos/',num2str(i)),'vid_sag_NN'))
	title('NN')

	close all;
end