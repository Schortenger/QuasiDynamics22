function vid_array = generate_video_array(object)

	% record video?
	if nargin < 2; vid_on = false; end;

	% Animates the execution of the contact optimization plan
	drawArrow = @(x,y) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0,'linewidth',1,'color','r')   
	drawArrow2 = @(x,y) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0,'linewidth',1,'color','b')    
	drawArrow3 = @(x,y) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0,'linewidth',1,'color','g')    

	% reads all relevant parameters of the problem
	r = object.traj.r(1:2,:);
	dr = object.traj.dr(1:2,:);
	th = object.traj.r(3,:);
	N_T = size(r,2);
	verts = object.v;

	vid_array = [];

	% draws the animation
	for t = 1:N_T
		h0 = figure('visible','off');
		clf(h0);
		hold on;

		for i = 1:length(object.env)
			x = object.env{i}.x;
			y = object.env{i}.y;
			pgon = polyshape(x,y);
			plot(pgon,'EdgeColor','red','FaceColor','red')
			hold on;
		end

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
		end

		pgon = polyshape(x,y);
		plot(pgon,'FaceAlpha',1.0,'FaceColor','blue','EdgeColor','blue')
		hold on;

		xlim([-0.2,0.2]);
		ylim([-0.2,0.2]);

		axis off
		% set(gcf,'color','w')
		set(gca,'color',[0 0 0])


		F = getframe(gcf);
		[X, ~] = frame2im(F);
		X = imresize(X,[100 100]);

		% figure(100)
		% image(X)
		% pause()

		vid_array = [vid_array, reshape(rescale(double(X),0.0,1.0),[1,3*100*100])];
	end
end