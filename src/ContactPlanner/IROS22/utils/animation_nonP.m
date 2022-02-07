function animation_nonP(object, environment, plan)

	% record video?
	close all;
	if nargin < 4; vid_on = true; end;
	if nargin < 5; trans = false; end;
	if nargin < 6; play_anim = false; end;
	% Animates the execution of the contact optimization plan
	drawArrow = @(x,y) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0,'linewidth',2,'color','r')   
	drawArrow2 = @(x,y) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0,'linewidth',2,'color',[86 178 29]/255)    
	drawArrow3 = @(x,y) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0,'linewidth',2,'color',[0 133 215]/255)
	drawArrow4 = @(x,y) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0,'linewidth',2,'color',[133 0 215]/255)
	drawArrow5 = @(x,y) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0,'linewidth',2,'color','k')    

	% reads all relevant parameters of the problem
	N_l = plan.N_l;
	N_c = plan.N_c;
	N_T = plan.N_T;
	r = object.traj.r(1:2,:);
	dr = object.traj.dr(1:2,:);
	ddr = object.traj.ddr(1:2,:);

	th = object.traj.r(3,:);
	verts = object.v;

	% draws the key frames
	h0 = figure(210)
	clf(h0);
	hold on;

	env = zeros(2,length(environment)+1);

	for i = 1:length(environment)
		% save environment vector
		env(:,i) = environment{i}(:,1);
	end
	env(:,end) = environment{end}(:,2);

	for t = 1:N_T
		% check if object is penetrating the environment
		penetrates = false;
		min_pen = -1;
		normal_pen = zeros(2,1);
		for i = 1:length(environment)
			rot = [cos(th(t)),-sin(th(t));sin(th(t)),cos(th(t))];
			tran = r(:,t);
			for v = 1:object.nv
				new_vert = tran + rot*verts(:,v);

				v1 = new_vert 			 - environment{i}(:,1);
				v2 = environment{i}(:,2) - environment{i}(:,1);
				proj1_2 = v1'*v2/(norm(v2) + 1e-6);
				vn = [0,-1;1,0]*v2/(norm(v2) + 1e-6);
				[in,on_] = inpolygon(new_vert(1),new_vert(2),[env(1,:)],[env(2,:)]);
				inside = 1 - max(in,on_);

				if proj1_2 > 0 && proj1_2 < norm(v2) && inside
					if v1'*vn < 0
						% inside
						% v1
						% v2
						% vn
						% new_vert
						if v1'*vn > min_pen
							normal_pen = vn*v1'*vn;

							% figure(1984)

							% plot(env(1,:),env(2,:)) % polygon
							% hold on;
							% plot(new_vert(1),new_vert(2),'r+')

							% pause()

							penetrates = true;
							t
							v1
							v2
							new_vert
							min_pen = v1'*vn
						end
						% t
						% pause()
					end
				end
			end
		end
		if penetrates
			r(:,t) = r(:,t) - normal_pen;

			for c = 1:N_c
				for l = 1:N_l
					plan.p(:,c,l,t,1) = plan.p(:,c,l,t,1) - normal_pen;
				end
			end
		end
	end

	% draws the environment
	pgon = polyshape(env(1,:)',env(2,:)');
	plot(pgon,'FaceAlpha',0.2,'FaceColor',[100 100 100]/255,'EdgeColor','black')
	hold on;

	for t = 1:N_T
		% transformation matrix
		rot = [cos(th(t)),-sin(th(t));sin(th(t)),cos(th(t))];
		tran = r(:,t);

		% draws the regions
		v = -0.2:0.002:0.2;  % plotting range from -5 to 5
		[x, y] = meshgrid(v);  % Get 2-D mesh for x and y based on r
		[x, y] = meshgrid(v);  % get 2-D mesh for x and y
		cond = ones(length(v)); % Initialize
		for re = 1:length(object.regions)
			A = object.regions{re}.A;
			b = object.regions{re}.b;
			% condition = (A(1,1,t)*x + A(1,2,t)*y < b(t));
			% cond(condition) = 0;
		end
		% surf(x, y, cond)
		% view(0,90)
		hold on
		% pause();		

		for v = 1:object.n_env_v
			new_vert = env(:,v);

			fv_1 = zeros(1,N_T);
			fv_2 = zeros(1,N_T);

			fv_1(1,:) = reshape(plan.f_v(1,v,:),[1,N_T]);
			fv_2(1,:) = reshape(plan.f_v(2,v,:),[1,N_T]);

			x_ = [new_vert(1), new_vert(1) + fv_1(t)/4]; 
			y_ = [new_vert(2), new_vert(2) + fv_2(t)/4];
			drawArrow5(x_,y_);
			hold on;
		end

		% draws the polygon
		x = [];
		y = [];
		for v = 1:object.nv
			new_vert = tran + rot*verts(:,v);

			x = [x, new_vert(1)];
			y = [y, new_vert(2)];

			fext_1 = zeros(1,N_T);
			fext_2 = zeros(1,N_T);

			fext_1(1,:) = reshape(plan.f_ext(1,v,:),[1,N_T]);
			fext_2(1,:) = reshape(plan.f_ext(2,v,:),[1,N_T]);

			x_ = [new_vert(1), new_vert(1) + fext_1(t)/4]; 
			y_ = [new_vert(2), new_vert(2) + fext_2(t)/4];
			drawArrow3(x_,y_);
			hold on;
		end

		pgon = polyshape(x,y);
		plot(pgon,'FaceAlpha',0.9,'FaceColor',[221 217 195]/255,'EdgeColor','black')
		hold on;

		for l = 1:N_l
			for c = 1:N_c

				p1 = zeros(1,N_T); p2 = zeros(1,N_T);
				f1 = zeros(1,N_T); f2 = zeros(1,N_T);

				p1(1,:) = plan.p(1,c,l,:,1);
				p2(1,:) = plan.p(2,c,l,:,1);

				f1(1,:) = plan.f(1,c,l,:,1);
				f2(1,:) = plan.f(2,c,l,:,1);

				viscircles([p1(1,t),p2(1,t)],0.002);
				hold on;

				x = [p1(t), p1(t) + f1(t)/8]; 
				y = [p2(t), p2(t) + f2(t)/8];
				
				drawArrow2(x,y);
				hold on;

				x = [r(1,t), r(1,t) + 2*dr(1,t)]; 
				y = [r(2,t), r(2,t) + 2*dr(2,t)];
				% drawArrow3(x,y);
				% hold on;
			end
		end

		xlim([-1,1]);
		ylim([-0.5,1]);
	end

	xlim([-1,1]);
	ylim([-0.5,1]);

	for t = 1:N_T
		for l = 1:N_l
			for c = 1:N_c

				p1 = zeros(1,N_T); p2 = zeros(1,N_T);
				f1 = zeros(1,N_T); f2 = zeros(1,N_T);

				p1(1,:) = plan.p(1,c,l,:,1);
				p2(1,:) = plan.p(2,c,l,:,1);

				f1(1,:) = plan.f(1,c,l,:,1);
				f2(1,:) = plan.f(2,c,l,:,1);

				viscircles([p1(1,t),p2(1,t)],0.002);
				hold on;
			end
		end
	end
	if trans == 1; col = 0.2; end;
	if trans == 0; col = 1.0; end;
	set(gca,'color',[col col col])
	set(gca,'XTickLabel',[]);
	set(gca,'YTickLabel',[]);
	box on

	if play_anim
		pause()
		mult = 10;
		time = linspace(0,1,N_T);
		t1 = linspace(0,1,mult*N_T);

		r = [interp1(time,r(1,:),t1); interp1(time,r(2,:),t1)];
		dr = [interp1(time,dr(1,:),t1); interp1(time,dr(2,:),t1)];
		ddr = [interp1(time,ddr(1,:),t1); interp1(time,ddr(2,:),t1)];

		th = interp1(time,th,t1);

		if vid_on
			name = input('movie file name: ','s')

			writerObj = VideoWriter(name);
			writerObj.FrameRate = 10;
			open(writerObj);
		end

		% does the animation
		h = figure(420)
		clf(h);
		hold on;
		for t = 1:mult*N_T
			% pause()
			% figure(420)
			clf(h)
			% figure(420)
			% transformation matrix
			rot = [cos(th(t)),-sin(th(t));sin(th(t)),cos(th(t))];
			tran = r(:,t);

			pgon = polyshape(env(1,:)',env(2,:)');
			plot(pgon,'FaceAlpha',0.2,'FaceColor',[100 100 100]/255,'EdgeColor','black')
			hold on;

			% draws the polygon
			x = [];
			y = [];
			for v = 1:object.nv
				new_vert = tran + rot*verts(:,v);

				x = [x, new_vert(1)];
				y = [y, new_vert(2)];
			end

			pgon = polyshape(x,y);
			plot(pgon,'FaceAlpha',0.9,'FaceColor',[221 217 195]/255,'EdgeColor','black')
			hold on;

			for v = 1:object.n_env_v
				new_vert = env(:,v);
				
				fv_1 = zeros(1,N_T);
				fv_2 = zeros(1,N_T);

				fv_1(1,:) = plan.f_v(1,v,:);
				fv_2(1,:) = plan.f_v(2,v,:);

				fv_1 = interp1(time,fv_1(1,:),t1);
				fv_2 = interp1(time,fv_2(1,:),t1);

				if abs(plan.f_v(1,v,ceil(t/mult))) < 1e-3
					fv_1 = fv_1*0;
				end

				if abs(plan.f_v(2,v,ceil(t/mult))) < 1e-3
					fv_2 = fv_2*0;
				end

				if ceil(t/mult) < N_T
					if abs(plan.f_v(1,v,ceil(t/mult)+1)) < 1e-3
						fv_1 = fv_1*0;
					end

					if abs(plan.f_v(2,v,ceil(t/mult)+1)) < 1e-3
						fv_2 = fv_2*0;
					end
				end

				if floor(t/mult) > 0
					if abs(plan.f_v(1,v,floor(t/mult))) < 1e-3
						fv_1 = fv_1*0;
					end

					if abs(plan.f_v(2,v,floor(t/mult))) < 1e-3
						fv_2 = fv_2*0;
					end
				end
				
				x_ = [new_vert(1), new_vert(1) + fv_1(t)/2]; 
				y_ = [new_vert(2), new_vert(2) + fv_2(t)/2];
				drawArrow5(x_,y_);
				hold on;
			end

			for v = 1:object.nv
				new_vert = tran + rot*verts(:,v);
				
				fext_1 = zeros(1,N_T);
				fext_2 = zeros(1,N_T);

				fext_1(1,:) = plan.f_ext(1,v,:);
				fext_2(1,:) = plan.f_ext(2,v,:);

				fext_1 = interp1(time,fext_1(1,:),t1);
				fext_2 = interp1(time,fext_2(1,:),t1);

				if abs(plan.f_ext(1,v,ceil(t/mult))) < 1e-3
					fext_1 = fext_1*0;
				end

				if abs(plan.f_ext(2,v,ceil(t/mult))) < 1e-3
					fext_2 = fext_2*0;
				end

				if ceil(t/mult) < N_T
					if abs(plan.f_ext(1,v,ceil(t/mult)+1)) < 1e-3
						fext_1 = fext_1*0;
					end

					if abs(plan.f_ext(2,v,ceil(t/mult)+1)) < 1e-3
						fext_2 = fext_2*0;
					end
				end

				if floor(t/mult) > 0
					if abs(plan.f_ext(1,v,floor(t/mult))) < 1e-3
						fext_1 = fext_1*0;
					end

					if abs(plan.f_ext(2,v,floor(t/mult))) < 1e-3
						fext_2 = fext_2*0;
					end
				end


				x_ = [new_vert(1), new_vert(1) + fext_1(t)/2]; 
				y_ = [new_vert(2), new_vert(2) + fext_2(t)/2];
				drawArrow(x_,y_);
				hold on;
			end

			% draws the floor
			for i = 1:length(object.env)
				x = object.env{i}.x;
				y = object.env{i}.y;
				pgon = polyshape(x,y);
				plot(pgon,'FaceAlpha',0.9,'FaceColor',[100 100 100]/255,'EdgeColor','black')
				hold on;
			end

			for l = 1:N_l
				for c = 1:N_c
					% interpolates force
					p1 = zeros(1,N_T); p2 = zeros(1,N_T);
					f1 = zeros(1,N_T); f2 = zeros(1,N_T);
					p1(1,:) = plan.p(1,c,l,:);
					p2(1,:) = plan.p(2,c,l,:);

					f1(1,:) = plan.f(1,c,l,:);
					f2(1,:) = plan.f(2,c,l,:);

					p1 = interp1(time,p1,t1);
					p2 = interp1(time,p2,t1);

					viscircles([p1(1,t),p2(1,t)],0.002);
					hold on;

					f1 = interp1(time,f1,t1);
					f2 = interp1(time,f2,t1);

					% pause()

					if abs(plan.f(1,c,l,ceil(t/mult))) < 1e-3
						f1 = f1*0;
					end

					if abs(plan.f(2,c,l,ceil(t/mult))) < 1e-3
						f2 = f2*0;
					end

					if ceil(t/mult) < N_T
						if abs(plan.f(1,c,l,ceil(t/mult)+1)) < 1e-3
							f1 = f1*0;
						end

						if abs(plan.f(2,c,l,ceil(t/mult)+1)) < 1e-3
							f2 = f2*0;
						end
					end

					if floor(t/mult) > 0
						if abs(plan.f(1,c,l,floor(t/mult))) < 1e-3
							f1 = f1*0;
						end

						if abs(plan.f(2,c,l,floor(t/mult))) < 1e-3
							f2 = f2*0;
						end
					end

					x = [p1(t), p1(t) + f1(t)/2]; 
					y = [p2(t), p2(t) + f2(t)/2];
					
					drawArrow2(x,y);
					hold on;

					x = [r(1,t), r(1,t) + object.m*ddr(1,t)/2]; 
					y = [r(2,t), r(2,t) + object.m*ddr(2,t)/2  + object.m*9.8/2];
					% drawArrow(x,y);
					% hold on;

					% pause()
				end
			end

			xlim([-1,1]);
			ylim([-0.5,1]);
			if trans == true; col = 0.2; end
			if trans == false; col = 1.0; end
			col = 0.2;
			set(gca,'color',[1 1 1])
			set(gca,'XTickLabel',[]);
			set(gca,'YTickLabel',[]);
			% axis off

			hold on;
			% pause()

			pause(0.0001);
			if vid_on
				frame = getframe(gcf);
				writeVideo(writerObj, frame);
			end
		end
		if vid_on; close(writerObj); end;
	end
end