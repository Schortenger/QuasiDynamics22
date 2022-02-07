clear all; close all; clc;

set_j = 99

res0 = load(strcat(strcat('./../data/AffordanceData/raw_',num2str(set_j)),'_2f_sq.mat'));

p_GT = csvread('../../fair-quasidyn/data/sim_p_GT_0_2f.csv');
p_Sim = csvread('../../fair-quasidyn/data/sim_p_full_0_2f.csv');
p_nn = csvread('../../fair-quasidyn/data/sim_p_nn_0_2f.csv');
p_cvx = csvread('../../fair-quasidyn/data/sim_p_mdr_0_2f.csv');
p_mdr = csvread('../../fair-quasidyn/data/sim_p_cvx_0_2f.csv');

r_GT = csvread('../../fair-quasidyn/data/sim_r_GT_0_2f.csv');
r_Sim = csvread('../../fair-quasidyn/data/sim_r_full_0_2f.csv');
r_nn = csvread('../../fair-quasidyn/data/sim_r_nn_0_2f.csv');
r_cvx = csvread('../../fair-quasidyn/data/sim_r_mdr_0_2f.csv');
r_mdr = csvread('../../fair-quasidyn/data/sim_r_cvx_0_2f.csv');

data = csvread(strcat(strcat('../../fair-quasidyn/data/data_',num2str(set_j)),'_2f_sq.csv'));

for i = 1:5
	pre = res0.data{i}.pre;
	task = pre.object;
	pre.vars.f.value = zeros(2,1,2,5);

	i

	% plots
	% figure(1)
	animation_sim(task,pre,1,true,false,true,false,strcat(strcat('../videos/',num2str(i)),'vid_GT'))
	title('GT')

	traj_1 = [p_GT(i, 1:5); p_GT(i, 11:15)];
	traj_2 = [p_GT(i, 6:10); p_GT(i, 16:20)];
	r = [r_GT(i, 1:5); r_GT(i, 6:10); asin((1/0.03)*r_GT(i, 11:15))];
	r_ref = [r_GT(i, 1:5); r_GT(i, 6:10); r_GT(i, 11:15)];

	p_data = zeros(2,1,2,5);
	p_data(:,1,1,:) = traj_1;
	p_data(:,1,2,:) = traj_2;

	task.traj.r = r;
	pre.vars.p.value = p_data;

	% plots
	% figure(2)
	animation_sim(task,pre,2,true,false,true,true,strcat(strcat('../videos/',num2str(i)),'vid_MIQP'))
	title('MIQP')

	traj_1 = [p_Sim(i, 1:5); p_Sim(i, 11:15)];
	traj_2 = [p_Sim(i, 6:10); p_Sim(i, 16:20)];
	r = [r_Sim(i, 1:5); r_Sim(i, 6:10); asin((1/0.03)*r_Sim(i, 11:15))];

	r_err1 = norm(r_ref - [r_Sim(i, 1:5); r_Sim(i, 6:10); r_Sim(i, 11:15)])

	p_data = zeros(2,1,2,5);
	p_data(:,1,1,:) = traj_1;
	p_data(:,1,2,:) = traj_2;

	task.traj.r = r;
	pre.vars.p.value = p_data;

	% plots
	% figure(3)
	animation_sim(task,pre,3,true,false,true,true,strcat(strcat('../videos/',num2str(i)),'vid_DDM'))
	title('SIM')

	traj_1 = [p_mdr(i, 1:5); p_mdr(i, 11:15)];
	traj_2 = [p_mdr(i, 6:10); p_mdr(i, 16:20)];
	r = [r_mdr(i, 1:5); r_mdr(i, 6:10); asin((1/0.03)*r_mdr(i, 11:15))];

	r_err2 = norm(r_ref - [r_mdr(i, 1:5); r_mdr(i, 6:10); r_mdr(i, 11:15)])

	p_data = zeros(2,1,2,5);
	p_data(:,1,1,:) = traj_1;
	p_data(:,1,2,:) = traj_2;

	task.traj.r = r;
	pre.vars.p.value = p_data;

	% plots
	% figure(3)
	animation_sim(task,pre,3,true,false,true,true,strcat(strcat('../videos/',num2str(i)),'vid_MDR'))
	title('MDR')

	traj_1 = [p_cvx(i, 1:5); p_cvx(i, 11:15)];
	traj_2 = [p_cvx(i, 6:10); p_cvx(i, 16:20)];
	r = [r_cvx(i, 1:5); r_cvx(i, 6:10); asin((1/0.03)*r_cvx(i, 11:15))];

	r_err3 = norm(r_ref - [r_cvx(i, 1:5); r_cvx(i, 6:10); r_cvx(i, 11:15)])

	p_data = zeros(2,1,2,5);
	p_data(:,1,1,:) = traj_1;
	p_data(:,1,2,:) = traj_2;

	task.traj.r = r;
	pre.vars.p.value = p_data;

	% plots
	% figure(3)
	animation_sim(task,pre,3,true,false,true,true,strcat(strcat('../videos/',num2str(i)),'vid_CVX'))
	title('CVX')

	traj_1 = [p_nn(i, 1:5); p_nn(i, 11:15)];
	traj_2 = [p_nn(i, 6:10); p_nn(i, 16:20)];
	r = [r_nn(i, 1:5); r_nn(i, 6:10); asin((1/0.03)*r_nn(i, 11:15))];

	r_err4 = norm(r_ref - [r_nn(i, 1:5); r_nn(i, 6:10); r_nn(i, 11:15)])

	p_data = zeros(2,1,2,5);
	p_data(:,1,1,:) = traj_1;
	p_data(:,1,2,:) = traj_2;

	task.traj.r = r;
	pre.vars.p.value = p_data;

	% plots
	% figure(3)
	animation_sim(task,pre,3,true,false,true,true,strcat(strcat('../videos/',num2str(i)),'vid_NN'))
	title('CVX')

	close all;
end