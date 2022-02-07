function animation_sim(object, plan, fig_n, vid_on, trans, play_anim, add_finger, vid_name)

	% record video?
	if nargin < 3; vid_on = false; end;
	if nargin < 4; trans = false; end;
	if nargin < 5; play_anim = false; end;
	if nargin < 6; add_finger = false; end;
	if nargin < 7; vid_name = './vid_unlabeled'; end;
	% Animates the execution of the contact optimization plan
	drawArrow = @(x,y) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0,'linewidth',2,'color','r');
	drawArrow2 = @(x,y) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0,'linewidth',2,'color',[86 178 29]/255);   
	drawArrow3 = @(x,y) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0,'linewidth',2,'color',[0 133 215]/255);

    vid_on = 1;
    play_anim = 1;
    
	% reads all relevant parameters of the problem
	N_l = plan.N_l;
	N_c = plan.N_c;
	N_T = plan.N_T;
	r = object.traj.r(1:2,:);
	dr = object.traj.dr(1:2,:);
	th = object.traj.r(3,:);
	verts = object.v;

	% draws the key frames
	if true
		h0 = figure('visible','off');
		clf(h0);
		hold on;

		x = [];
		y = [];
		for v = 1:object.nv
			new_vert = verts(:,v);

			x = [x, new_vert(1)];
			y = [y, new_vert(2)];
		end

		pgon = polyshape(x,y);
		plot(pgon,'FaceAlpha',0.9,'FaceColor',[239 188 64]/255,'EdgeColor','black','LineWidth',3)
		hold on;

		exportgraphics(h0,strcat(vid_name,'_shape.pdf'),'ContentType','vector')

		h0 = figure('visible','off');
		clf(h0);
		hold on;

		% draws the floor
		for i = 1:length(object.env)
			x = object.env{i}.x;
			y = object.env{i}.y;
			pgon = polyshape(x,y);
			plot(pgon,'FaceAlpha',0.9,'FaceColor',[70 161 248]/255,'EdgeColor','black','LineWidth',3)
			hold on;
		end

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

			% draws the polygon
			x = [];
			y = [];
			for v = 1:object.nv
				new_vert = tran + rot*verts(:,v);

				x = [x, new_vert(1)];
				y = [y, new_vert(2)];

				fext_1 = zeros(1,N_T);
				fext_2 = zeros(1,N_T);

				fext_1(1,:) = plan.vars.f_ext.value(1,v,:,1);
				fext_2(1,:) = plan.vars.f_ext.value(2,v,:,1);

				x_ = [new_vert(1), new_vert(1) + fext_1(t)/4]; 
				y_ = [new_vert(2), new_vert(2) + fext_2(t)/4];
				% drawArrow3(x_,y_);
				hold on;
			end

			pgon = polyshape(x,y);
			plot(pgon,'FaceAlpha',0.9*t/N_T,'FaceColor',[239 188 64]/255,'EdgeColor','black','LineWidth',3)
			hold on;

			for l = 1:N_l
				for c = 1:N_c

					p1 = zeros(1,N_T); p2 = zeros(1,N_T);
					f1 = zeros(1,N_T); f2 = zeros(1,N_T);

					p1(1,:) = plan.vars.p.value(1,c,l,:,1);
					p2(1,:) = plan.vars.p.value(2,c,l,:,1);

					f1(1,:) = plan.vars.f.value(1,c,l,:,1);
					f2(1,:) = plan.vars.f.value(2,c,l,:,1);
					if add_finger
						% viscircles([p1(1,t),p2(1,t)],0.005);
						range_t = 0.05:0.5:2*pi;
						x1 = 0.005*cos(range_t);
						y1 = 0.005*sin(range_t);
						pgon = polyshape(x1 + p1(1,t),y1 + p2(1,t));
						plot(pgon,'FaceAlpha',t/N_T,'FaceColor',[1 0 0],'EdgeColor','red','LineWidth',0.05)
					end
					hold on;

					x = [p1(t), p1(t) + f1(t)/8]; 
					y = [p2(t), p2(t) + f2(t)/8];
					
					% drawArrow2(x,y);
					hold on;

					x = [r(1,t), r(1,t) + 2*dr(1,t)]; 
					y = [r(2,t), r(2,t) + 2*dr(2,t)];
					% drawArrow3(x,y);
					hold on;
				end
			end

			xlim([-0.2,0.2]);
			ylim([-0.1,0.2]);
		end

		xlim([-0.2,0.2]);
		ylim([-0.1,0.2]);

		% for t = 1:N_T
		% 	for l = 1:N_l
		% 		for c = 1:N_c

		% 			p1 = zeros(1,N_T); p2 = zeros(1,N_T);
		% 			f1 = zeros(1,N_T); f2 = zeros(1,N_T);

		% 			p1(1,:) = plan.vars.p.value(1,c,l,:,1);
		% 			p2(1,:) = plan.vars.p.value(2,c,l,:,1);

		% 			f1(1,:) = plan.vars.f.value(1,c,l,:,1);
		% 			f2(1,:) = plan.vars.f.value(2,c,l,:,1);

		% 			if add_finger
		% 				viscircles([p1(1,t),p2(1,t)],0.001);
		% 			end
		% 			hold on;
		% 		end
		% 	end
		% end
		if trans == 1; col = 0.2; end
		if trans == 0; col = 1.0; end
		set(gca,'color',[col col col])
		set(gca,'XTickLabel',[]);
		set(gca,'YTickLabel',[]);
		box on

		exportgraphics(h0,strcat(vid_name,'.pdf'),'ContentType','vector')
		savefig(h0,strcat(vid_name,'.fig'))

		pause(10)
	end 

	if play_anim
		% pause()

		time = linspace(0,1,N_T);
		t1 = linspace(0,1,10*N_T);
		r = [interp1(time,r(1,:),t1); interp1(time,r(2,:),t1)];
		dr = [interp1(time,dr(1,:),t1); interp1(time,dr(2,:),t1)];
		th = interp1(time,th,t1);

		if vid_on
			if nargin < 7
				name = input('movie file name: ','s')
			else
				name = vid_name;
			end

			writerObj = VideoWriter(name,'MPEG-4');
			writerObj.FrameRate = 10;
			open(writerObj);
		end

		% does the animation
		for t = 1:10*N_T
			h = figure('Name','video','visible','off');
			clf(h)
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

			% draws the polygon
			x = [];
			y = [];
			for v = 1:object.nv
				new_vert = tran + rot*verts(:,v);

				x = [x, new_vert(1)];
				y = [y, new_vert(2)];
			end
			pgon = polyshape(x,y);
			plot(pgon,'FaceAlpha',1,'FaceColor',[239 188 64]/255,'EdgeColor','black','LineWidth',3)
			hold on;

			% for v = 1:object.nv
			% 	new_vert = tran + rot*verts(:,v);
				
			% 	fext_1 = zeros(1,N_T);
			% 	fext_2 = zeros(1,N_T);

			% 	fext_1(1,:) = plan.vars.f_ext.value(1,v,:);
			% 	fext_2(1,:) = plan.vars.f_ext.value(2,v,:);

			% 	fext_1 = interp1(time,fext_1(1,:),t1);
			% 	fext_2 = interp1(time,fext_2(1,:),t1);
			% 	if t > 10
			% 		if fext_1(t-10) == 0
			% 			fext_1 = fext_1*0
			% 		end
			% 		if fext_2(t-10) == 0
			% 			fext_2 = fext_2*0
			% 		end
			% 	end
			% 	x_ = [new_vert(1), new_vert(1) + fext_1(t)/2]; 
			% 	y_ = [new_vert(2), new_vert(2) + fext_2(t)/2];
			% 	% drawArrow(x_,y_);
			% 	hold on;
			% end

			% draws the floor
			for i = 1:length(object.env)
				x = object.env{i}.x;
				y = object.env{i}.y;
				pgon = polyshape(x,y);
				plot(pgon,'FaceAlpha',0.9,'FaceColor',[70 161 248]/255,'EdgeColor','black','LineWidth',3)
				hold on;
			end

			for l = 1:N_l
				for c = 1:N_c

					p1 = zeros(1,N_T); p2 = zeros(1,N_T);
					f1 = zeros(1,N_T); f2 = zeros(1,N_T);
					p1(1,:) = plan.vars.p.value(1,c,l,:);
					p2(1,:) = plan.vars.p.value(2,c,l,:);

					% f1(1,:) = plan.vars.f.value(1,c,l,:);
					% f2(1,:) = plan.vars.f.value(2,c,l,:);

					p1 = interp1(time,p1,t1);
					p2 = interp1(time,p2,t1);
					if add_finger
						viscircles([p1(1,t),p2(1,t)],0.005);
						hold on;
					end

					% f1 = interp1(time,f1,t1);
					% f2 = interp1(time,f2,t1);

					% x = [p1(t), p1(t) + f1(t)/2]; 
					% y = [p2(t), p2(t) + f2(t)/2];
					
					% % drawArrow2(x,y);
					% hold on;

					% x = [r(1,t), r(1,t) + 2*dr(1,t)]; 
					% y = [r(2,t), r(2,t) + 2*dr(2,t)];
					% % drawArrow3(x,y);
					% hold on;
				end
			end

			xlim([-0.2,0.2]);
			ylim([-0.2,0.2]);
			if trans == true; col = 0.2; end
			if trans == false; col = 1.0; end
			axis off
			set(gca,'color',[0 0 0])
			
			hold on;
			if vid_on
				frame = getframe(gcf);
				writeVideo(writerObj, frame);
			end
		end
		if vid_on; close(writerObj); end;
	end
end