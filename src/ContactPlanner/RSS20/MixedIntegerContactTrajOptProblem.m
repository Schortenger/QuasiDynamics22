classdef MixedIntegerContactTrajOptProblem < Quad_MixedIntegerConvexProgram
% Developed by Bernardo Aceituno-C (MIT MCube Lab)
  properties
    object
    N_c
    N_l
    N_T
    N_f
    N_r
    N_env
    M = 6
    idx = 0
    idx_sos2 = 0
    dt = 0.1
    McCormick = 1
    modes = 0
  end

  methods
    function obj = MixedIntegerContactTrajOptProblem(object,n_c)
      % Constructs the optimization problem and declares the variables for each contact
      % (This version operates in general 3D environments)
      % @param object: structure with elements:
      %               - object.v        : a 3XNv vector with the position of each vertex wrt to shape    
      %               - object.regions  : a set of convex regions that represent the complement 
      %                                   space of the shape through time
      %               - object.lines    : a set of line segments covering the shape
      %               - object.nv       : number of vertices in all the polygons of the shape
      %               - object.traj     : a 3xN_T vector with the object CoM = [x(t),y(t),th(t)]'
      %               - object.ext_f    : a N_v cell with the contact constraints through time
      % @param n_links: links on the gripper

      assert(nargin >= 1)

      if nargin < 2; n_c = 2; end;

      % parses the inputs
      obj.object = object;
      obj.N_T = size(object.traj.r,2);
      obj.N_l = 1;
      obj.N_c = n_c;
      obj.N_f = length(object.lines);
      obj.N_r = length(object.regions);

      obj.N_env = length(object.env);

      % contact locatops
      obj = obj.addVariable('p', 'C', [3, obj.N_c, obj.N_l, obj.N_T], -inf, inf);
      obj = obj.addVariable('dp', 'C', [3, obj.N_c, obj.N_l, obj.N_T], -inf, inf);
      obj = obj.addVariable('ddp', 'C', [3, obj.N_c, obj.N_l, obj.N_T], -inf, inf);

      obj = obj.addVariable('rel_ddp', 'C', [3, obj.N_c, obj.N_l, obj.N_T], -inf, inf);

      % obj = obj.addVariable('ddth_aux', 'C', [1, obj.N_T], -inf, inf);
      
      % time integration
      obj = obj.addVariable('dr', 'C', [6, obj.N_T], -inf, inf);
      obj = obj.addVariable('ddr', 'C', [6, obj.N_T], -inf, inf);
      obj = obj.addVariable('dt', 'C', [1, obj.N_T], 0, 1);

      % contact forces
      obj = obj.addVariable('f', 'C', [3, obj.N_c, obj.N_l, obj.N_T], -inf, inf);
      obj = obj.addVariable('alpha_f', 'C', [1, obj.N_c, obj.N_l, obj.N_T], 0, inf);
      obj = obj.addVariable('alpha_min', 'C', [1, obj.N_T], 0, inf);
      % obj = obj.addVariable('alpha_ext', 'C', [1, obj.object.nv, obj.N_T], 0, inf);
      obj = obj.addVariable('tau_aux', 'C', [6, obj.N_c, obj.N_l, obj.N_T], -inf, inf);
      obj = obj.addVariable('f_ext', 'C', [3, obj.object.nv, obj.N_T], -1, 1);
    end

    function obj = addSOS2Constraint(x,y,x_rang,y_rang)
      % size of the range
      n = length(x_rang);

      % weights and binary variables
      lambda = strcat('sos2_',string(obj.idx_sos2));
      binvar = strcat('binsos2_',string(obj.idx_sos2));
      obj.idx_sos2 = obj.idx_sos2 + 1;

      % adds the varialbes for the problem
      obj = obj.addVariable(lambda, 'C', [n,1], 0, 1);
      obj = obj.addVariable(binvar, 'B', [n-1,1], 0, 1);

      % all weights must add to 1
      Aeq = sparse(2, obj.nv)
      beq = [1;1];

      Aeq(1,obj.vars.(lambda).i(:,1)) = 1;
      Aeq(2,obj.vars.(binvar).i(:,1)) = 1;

      obj = obj.addLinearConstraints([], [], Aeq, beq);

      % value assignment
      Aeq = sparse(2, obj.nv)
      beq = zeros(2,1);

      Aeq(1,x) = -1;
      Aeq(1,obj.vars.(lambda).i(:,1)) = x_rang;
      
      Aeq(2,y) = -1;
      Aeq(2,obj.vars.(lambda).i(:,1)) = y_rang;

      obj = obj.addLinearConstraints([], [], Aeq, beq);

      % at most two consecutive weights can add to 1
      for i = 1:n-1
        Ai = sparse(2, obj.nv)
        bi = [1;-1]+3;

        Ai(1,obj.vars.(lambda).i(i:i+1,1)) = 1;
        Ai(2,obj.vars.(lambda).i(i:i+1,1)) = -1;

        Ai(:,obj.vars.(binvar).i(i,1)) = 3;

        obj = obj.addLinearConstraints(Ai, bi, [], []);
      end
    end

    function obj = addBilinearSOS2Constraint(obj,w,x,y,x_ra,y_ra)
      % receives the indexes for variables w, x and y.
      % adds constrait w ~= x*y

      if nargin < 6
        y_ra = [-1,1];
      end

      if nargin < 5
        x_ra = [-1,1];
      end

      M = obj.M+1;

      % adds the binary variables for assignment
      vars_1 = strcat('mc_1_',string(obj.idx));
      vars_2 = strcat('mc_2_',string(obj.idx));

      gamma_ = strcat('gamma_',string(obj.idx));

      alpha_ = strcat('alpha_',string(obj.idx));
      beta_ = strcat('beta_',string(obj.idx));

      obj.idx = obj.idx + 1;
      
      obj = obj.addVariable(vars_1,'B', [M-1, 1], 0, 1);
      obj = obj.addVariable(vars_2,'B', [M-1, 1], 0, 1); 
      
      obj = obj.addVariable(gamma_,'C', [M, M], 0, 1);
      obj = obj.addVariable(alpha_,'C', [M, 1], 0, 1);
      obj = obj.addVariable(beta_,'C', [M, 1], 0, 1);

      % segmentation of x
      x_r = linspace(x_ra(1),x_ra(2),M);
      y_r = linspace(y_ra(1),y_ra(2),M);

      % adds the assignment constraints
      Aeq = sparse(2, obj.nv);
      beq = zeros(2,1);

      Aeq(1,x) = -1;
      Aeq(1,obj.vars.(alpha_).i(:,1)) = x_r(:);
      
      Aeq(2,y) = -1;
      Aeq(2,obj.vars.(beta_).i(:,1)) = y_r(:);
      
      obj = obj.addLinearConstraints([], [], Aeq, beq);

      % adds the approximation constraints
      Aeq = sparse(1, obj.nv);
      beq = 0;

      Aeq(:,w) = -1;

      for i = 1:M
        for j = 1:M
          Aeq(:,obj.vars.(gamma_).i(i,j)) = x_r(i)*y_r(j);
        end
      end

      obj = obj.addLinearConstraints([], [], Aeq, beq);

      % SOS2 constraint in alpha and beta
      Aeq = sparse(5, obj.nv);
      beq = ones(5,1);

      Aeq(1,obj.vars.(vars_1).i(:,1)) = 1;
      Aeq(2,obj.vars.(alpha_).i(:,1)) = 1;

      Aeq(3,obj.vars.(vars_2).i(:,1)) = 1;
      Aeq(4,obj.vars.(beta_).i(:,1)) = 1;

      Aeq(5,obj.vars.(gamma_).i(:,:)) = 1;

      obj = obj.addLinearConstraints([], [], Aeq, beq);

      for i = 1:M-1
        Ai = sparse(2, obj.nv);
        bi = [1;-1]+3;

        Ai(1,obj.vars.(alpha_).i(i:i+1,1)) = 1;
        Ai(2,obj.vars.(alpha_).i(i:i+1,1)) = -1;

        Ai(:,obj.vars.(vars_1).i(i,1)) = 3;

        obj = obj.addLinearConstraints(Ai, bi, [], []);

        Ai = sparse(2, obj.nv);
        bi = [1;-1]+3;

        Ai(1,obj.vars.(beta_).i(i:i+1,1)) = 1;
        Ai(2,obj.vars.(beta_).i(i:i+1,1)) = -1;

        Ai(:,obj.vars.(vars_2).i(i,1)) = 3;

        obj = obj.addLinearConstraints(Ai, bi, [], []);
      end

      % composition of gamma
      for i = 1:M
        % sum of columns adds to alphas
        Aeq = sparse(1,obj.nv);
        beq = 0;

        Aeq(1,obj.vars.(gamma_).i(i,:)) = 1;
        Aeq(1,obj.vars.(alpha_).i(i,1)) = -1;
        
        obj = obj.addLinearConstraints([], [], Aeq, beq);

        % sum of rows adds to betas
        Aeq = sparse(1,obj.nv);
        beq = 0;

        Aeq(1,obj.vars.(gamma_).i(:,i)) = 1;
        Aeq(1,obj.vars.(beta_).i(i,1)) = -1;
        
        obj = obj.addLinearConstraints([], [], Aeq, beq);
      end
    end

    function obj = addMcCormickEnvelopeConstraints(obj,w,x,y,x_r)
      if nargin < 5
        x_r = [-1,1];
      end
      % receives the indexes for variables w, x and y.
      % adds constrait w ~= x*y
      M = obj.M;

      % adds the binary variables for assignment
      vars_1 = strcat('mc1_',string(obj.idx));
      vars_2 = strcat('mc2_',string(obj.idx));
      obj.idx = obj.idx + 1;
      obj = obj.addVariable(vars_1, 'B', [M, 1], 0, 1);
      obj = obj.addVariable(vars_2, 'B', [M, 1], 0, 1);

      % segmentation of x
      x_ra = linspace(x_r(1),x_r(2),M+1);
      x_u = x_ra(2:end);
      x_l = x_ra(1:end-1);

      % the equality must be in one of the regions
      Aeq = sparse(1, obj.nv);
      beq = 1;

      Aeq(1,obj.vars.(vars_1).i(:,1)) = 1;
      Aeq(1,obj.vars.(vars_2).i(:,1)) = 1;

      obj = obj.addLinearConstraints([], [], Aeq, beq);

      % large K
      K = 10;

      % adds the envelope constraints for v >= 0
      for i = 1:M
        H1 = [0, x_l(i), -1;... 
              1, x_u(i), -1;... 
             -1, -x_l(i), 1;... 
              0, -x_u(i), 1;... 
              1, 0, 0;... 
             -1, 0, 0];

        Ai = sparse(6, obj.nv);
        bi = [0; x_u(i); -x_l(i); 0; x_u(i); -x_l(i)] + K;

        Ai(:,[x,y,w]) = H1;
        Ai(:,obj.vars.(vars_1).i(i,1)) = K;

        obj = obj.addLinearConstraints(Ai, bi, [], []);
      end

      % adds the envelope constraints for v < 0
      for i = 1:M
        H2 = [0, -x_l(i), -1;... 
             -1, -x_u(i), -1;... 
              1, x_l(i), 1;...
              0, x_u(i), 1;...
             -1, 0, 0;...
              1, 0, 0];
        
        Ai = sparse(6, obj.nv);
        bi = [0; x_u(i); -x_l(i); 0; x_u(i); -x_l(i)] + K;

        Ai(:,[x,y,w]) = H2;
        Ai(:,obj.vars.(vars_2).i(i,1)) = K;

        obj = obj.addLinearConstraints(Ai, bi, [], []);
      end
    end

    function obj = addDynamicsConstraints(obj)
      % constrains the object dynamics to be consistent with contact forces
      % and gravity

      % parameters
      m = obj.object.m; I = obj.object.I;
      g = obj.object.g; traj = obj.object.traj;

      M = 100;

      % constrains translational dynamics
      for t = 1:obj.N_T
        % computes relative acceleration
        for c = 1:obj.N_c
          for l = 1:obj.N_l
            Aeq = sparse(3,obj.nv);
            beq = traj.ddr(1:3,t);

            rot_mat = eul2rotm([traj.r(4,t),traj.r(5,t),traj.r(6,t)]);

            Aeq(:,obj.vars.ddp.i(:,c,l,t)) = inv(rot_mat);
            Aeq(:,obj.vars.rel_ddp.i(:,c,l,t)) = -eye(3);

            obj = obj.addLinearConstraints([], [], Aeq, beq);
          end
        end

        Aeq = sparse(3,obj.nv);
        beq = m*traj.ddr(1:3,t) + [0;0;m*g];

        % adds environment force constraints
        for v = 1:obj.object.nv
          Aeq(:,obj.vars.f_ext.i(:,v,t)) = eye(3);
        end

        % adds contact force constraints
        for c = 1:obj.N_c
          for l = 1:obj.N_l
            Aeq(:,obj.vars.f.i(:,c,l,t)) = eye(3);
          end
        end

        obj = obj.addLinearConstraints([], [], Aeq, beq);
      end

      % constrains rotational dynamics
      for t = 1:obj.N_T
        Aeq = sparse(3,obj.nv);
        beq = I*traj.ddr(4:6,t);

        % adds environment torque constraints 
        for v = 1:obj.object.nv          
          rot_mat = eul2rotm([traj.r(4,t),traj.r(5,t),traj.r(6,t)]);
          dp = rot_mat*(obj.object.v(:,v));
          Aeq(:,obj.vars.f_ext.i(:,v,t)) = [0,-dp(3),dp(2);dp(3),0,-dp(1);-dp(2),dp(1),0];
        end

        % adds contact torque constraints
        for c = 1:obj.N_c
          for l = 1:obj.N_l
            Aeq(1,obj.vars.tau_aux.i(1,c,l,t)) = 1;
            Aeq(1,obj.vars.tau_aux.i(2,c,l,t)) = -1;

            Aeq(1,obj.vars.f.i(2,c,l,t)) = traj.r(3,t);
            Aeq(1,obj.vars.f.i(3,c,l,t)) = -traj.r(2,t);

            Aeq(2,obj.vars.tau_aux.i(3,c,l,t)) = 1;
            Aeq(2,obj.vars.tau_aux.i(4,c,l,t)) = -1;

            Aeq(2,obj.vars.f.i(3,c,l,t)) = traj.r(1,t);
            Aeq(2,obj.vars.f.i(1,c,l,t)) = -traj.r(3,t);

            Aeq(3,obj.vars.tau_aux.i(5,c,l,t)) = 1;
            Aeq(3,obj.vars.tau_aux.i(6,c,l,t)) = -1;

            Aeq(3,obj.vars.f.i(1,c,l,t)) = traj.r(2,t);
            Aeq(3,obj.vars.f.i(2,c,l,t)) = -traj.r(1,t);
          end
        end

        % Aeq(:,obj.vars.ddth_aux.i(1,t)) = 1;

        obj = obj.addLinearConstraints([], [], Aeq, beq);

        % computes bilinear terms
        for c = 1:obj.N_c
          for l = 1:obj.N_l
            w1 = obj.vars.tau_aux.i(1,c,l,t);
            w2 = obj.vars.tau_aux.i(2,c,l,t);
            w3 = obj.vars.tau_aux.i(3,c,l,t);
            w4 = obj.vars.tau_aux.i(4,c,l,t);
            w5 = obj.vars.tau_aux.i(5,c,l,t);
            w6 = obj.vars.tau_aux.i(6,c,l,t);

            x1 = obj.vars.p.i(1,c,l,t);
            x2 = obj.vars.p.i(2,c,l,t);
            x3 = obj.vars.p.i(3,c,l,t);

            f1 = obj.vars.f.i(1,c,l,t);
            f2 = obj.vars.f.i(2,c,l,t); 
            f3 = obj.vars.f.i(3,c,l,t); 

            if obj.McCormick
              obj = obj.addMcCormickEnvelopeConstraints(w1,x3,f2,[-.2,.2]);
              obj = obj.addMcCormickEnvelopeConstraints(w2,x2,f3,[-.2,.2]);

              obj = obj.addMcCormickEnvelopeConstraints(w3,x3,f1,[-.2,.2]);
              obj = obj.addMcCormickEnvelopeConstraints(w4,x1,f3,[-.2,.2]);

              obj = obj.addMcCormickEnvelopeConstraints(w5,x1,f2,[-.2,.2]);
              obj = obj.addMcCormickEnvelopeConstraints(w6,x2,f1,[-.2,.2]);
            else
              % obj = obj.addBilinearSOS2Constraint(w1,x1,y1,[-0.1,0.1],[-0.06,0.06]);
              % obj = obj.addBilinearSOS2Constraint(w2,x2,y2,[0,0.2],[-0.03,0.03]);
              error('not implemented yet')
            end
          end
        end
      end

      % performs euler integration of the motion of each contact
      for t = 2:obj.N_T
        for c = 1:obj.N_c
          for l = 1:obj.N_l
            % integrates velocity
            Aeq = sparse(3,obj.nv);
            beq = zeros(3,1);

            Aeq(:,obj.vars.p.i(:,c,l,t)) = eye(3);
            Aeq(:,obj.vars.p.i(:,c,l,t-1)) = -eye(3);

            Aeq(:,obj.vars.dp.i(:,c,l,t)) = -obj.dt*eye(3);

            obj = obj.addLinearConstraints([], [], Aeq, beq);

            % integrates acceleration
            Aeq = sparse(3,obj.nv);
            beq = zeros(3,1);

            Aeq(:,obj.vars.dp.i(:,c,l,t)) = eye(3);
            Aeq(:,obj.vars.dp.i(:,c,l,t-1)) = -eye(3);

            Aeq(:,obj.vars.ddp.i(:,c,l,t)) = -obj.dt*eye(3);

            obj = obj.addLinearConstraints([], [], Aeq, beq);
          end
        end
      end

      % % defines initial conditions 
      Aeq = sparse(6,obj.nv);
      beq = zeros(6,1);
      
      Aeq(1:3,obj.vars.ddp.i(:,c,l,1)) = eye(3);
      Aeq(4:6,obj.vars.dp.i(:,c,l,1)) = eye(3);
      
      obj = obj.addLinearConstraints([], [], Aeq, beq);
    end

    function obj = addContactConstraints(obj)
      % Constrains reaction forces to be active only during contact
      % and forces to lie in the respective friction cone
      obj = obj.addVariable('L', 'B', [obj.N_f, obj.N_c, obj.N_l, obj.N_T], 0, 1);

      for f = 1:obj.N_f
        lambda = strcat('lambda_',string(f));
        obj = obj.addVariable(lambda, 'C', [obj.object.lines{f}.nv, obj.N_c, obj.N_l, obj.N_T], 0, 1);
      end

      obj = obj.addVariable('weight', 'C', [4, obj.N_c, obj.N_l, obj.N_T], 0, .1);
      obj = obj.addVariable('weight1', 'C', [4, obj.N_c, obj.N_l, obj.N_T], 0, .1);
      obj = obj.addVariable('weight_ext', 'C', [4, obj.object.nv, obj.N_T], 0, 0.5);

      % big-M
      M = 10;

      % activation of forces when in contact
      for t = 1:obj.N_T
        % computes the minimum margin
        for v = 1:obj.object.nv
          Ai = sparse(6,obj.nv);
          bi = zeros(6,1);

          Ai(1:3,obj.vars.f_ext.i(:,v,t)) = eye(3);
          Ai(1:3,obj.vars.weight_ext.i(1,v,t)) = -obj.object.ext_f{v}.fc1(:,t);
          Ai(1:3,obj.vars.weight_ext.i(2,v,t)) = -obj.object.ext_f{v}.fc2(:,t);
          Ai(1:3,obj.vars.weight_ext.i(3,v,t)) = -obj.object.ext_f{v}.fc3(:,t);
          Ai(1:3,obj.vars.weight_ext.i(4,v,t)) = -obj.object.ext_f{v}.fc4(:,t);

          Ai(4:6,obj.vars.f_ext.i(:,v,t)) = -eye(3);
          Ai(4:6,obj.vars.weight_ext.i(1,v,t)) = obj.object.ext_f{v}.fc1(:,t);
          Ai(4:6,obj.vars.weight_ext.i(2,v,t)) = obj.object.ext_f{v}.fc2(:,t);
          Ai(4:6,obj.vars.weight_ext.i(3,v,t)) = obj.object.ext_f{v}.fc3(:,t);
          Ai(4:6,obj.vars.weight_ext.i(4,v,t)) = obj.object.ext_f{v}.fc4(:,t);
          obj = obj.addLinearConstraints(Ai, bi, [], []);
        end

        for c = 1:obj.N_c
          for l = 1:obj.N_l
            Ai = sparse(1,obj.nv);
            bi = M;
            Ai(1,obj.vars.alpha_min.i(1,t)) = 1;
            Ai(1,obj.vars.alpha_f.i(1,c,l,t)) = -1;
            Ai(1,obj.vars.L.i(:,c,l,t)) = M;
            obj = obj.addLinearConstraints(Ai, bi, [], []);
 
            % complementarity constraint
            Aeq = sparse(1,obj.nv);
            beq = ones(1,1);

            for f = 1:obj.N_f
              lambda = strcat('lambda_',string(f));
              Aeq(1,obj.vars.(lambda).i(:,c,l,t)) = 1;
            end

            obj = obj.addLinearConstraints([], [], Aeq, beq);

            % if in contact, lambda does not change between time-steps
            if t < obj.N_T
              for f = 1:obj.N_f
                lambda = strcat('lambda_',string(f));
                Ai = sparse(2*obj.object.lines{f}.nv,obj.nv);
                bi = M*ones(2*obj.object.lines{f}.nv,1);
                Ai(1:obj.object.lines{f}.nv,obj.vars.(lambda).i(:,c,l,t)) = eye(obj.object.lines{f}.nv);
                Ai(1:obj.object.lines{f}.nv,obj.vars.(lambda).i(:,c,l,t+1)) = -eye(obj.object.lines{f}.nv);

                Ai(obj.object.lines{f}.nv+1:end,obj.vars.(lambda).i(:,c,l,t)) = -eye(obj.object.lines{f}.nv);
                Ai(obj.object.lines{f}.nv+1:end,obj.vars.(lambda).i(:,c,l,t+1)) = eye(obj.object.lines{f}.nv);

                Ai(1:4,obj.vars.L.i(:,c,l,t)) = M;

                obj = obj.addLinearConstraints(Ai, bi, [], []);
              end
            end

            % forces are zero unless in contact with a facet
            Ai = sparse(6,obj.nv);
            bi = zeros(6,1);

            Ai(1:3,obj.vars.f.i(:,c,l,t)) = eye(3);
            Ai(4:6,obj.vars.f.i(:,c,l,t)) = -eye(3);

            Ai(1:6,obj.vars.L.i(:,c,l,t)) = -M;

            obj = obj.addLinearConstraints(Ai, bi, [], []);

            % can only contact in one facet
            Ai = sparse(1,obj.nv);
            bi = ones(1,1);
            Ai(1,obj.vars.L.i(:,c,l,t)) = 1;
            obj = obj.addLinearConstraints(Ai, bi, [], []);

            for f = 1:obj.N_f
              % facet assignment constraint
              Ai = sparse(6,obj.nv);
              bi = M*ones(6,1);

              lambda = strcat('lambda_',string(f));

              Ai(1:3,obj.vars.p.i(:,c,l,t)) = eye(3);
              for v = 1:obj.object.lines{f}.nv
                Ai(1:3,obj.vars.(lambda).i(v,c,l,t)) = -obj.object.lines{f}.v(:,v,t);
              end

              Ai(4:6,obj.vars.p.i(:,c,l,t)) = -eye(3);
              for v = 1:obj.object.lines{f}.nv
                Ai(4:6,obj.vars.(lambda).i(v,c,l,t)) = obj.object.lines{f}.v(:,v,t);
              end
                            
              Ai(:,obj.vars.L.i(f,c,l,t)) = M;

              obj = obj.addLinearConstraints(Ai, bi, [], []);
              
              if t < obj.N_T
                Ai = sparse(1,obj.nv);
                bi = 0;
                
                Ai(1,obj.vars.L.i(f,c,l,t)) = -1;
                Ai(1,obj.vars.L.i(f,c,l,t+1)) = 1;

                obj = obj.addLinearConstraints(Ai, bi, [], []);
              end
              if obj.modes == 0
                % friction cone constraint
                Ai = sparse(6,obj.nv);
                bi = M*ones(6,1);

                Ai(1:3,obj.vars.f.i(:,c,l,t)) = eye(3);
                % Ai(1:2,obj.vars.alpha_f.i(1,c,l,t)) = -(obj.object.lines{f}.fc1(:,t)+obj.object.lines{f}.fc2(:,t))/2;
                Ai(1:3,obj.vars.weight1.i(1,c,l,t)) = -obj.object.lines{f}.fc1(:,t);
                Ai(1:3,obj.vars.weight1.i(2,c,l,t)) = -obj.object.lines{f}.fc2(:,t);
                Ai(1:3,obj.vars.weight1.i(3,c,l,t)) = -obj.object.lines{f}.fc3(:,t);
                Ai(1:3,obj.vars.weight1.i(4,c,l,t)) = -obj.object.lines{f}.fc4(:,t);
                
                Ai(4:6,obj.vars.f.i(:,c,l,t)) = -eye(3);
                % Ai(3:4,obj.vars.alpha_f.i(1,c,l,t)) = (obj.object.lines{f}.fc1(:,t)+obj.object.lines{f}.fc2(:,t))/2;
                Ai(4:6,obj.vars.weight1.i(1,c,l,t)) = obj.object.lines{f}.fc1(:,t);
                Ai(4:6,obj.vars.weight1.i(2,c,l,t)) = obj.object.lines{f}.fc2(:,t);
                Ai(4:6,obj.vars.weight1.i(3,c,l,t)) = obj.object.lines{f}.fc3(:,t);
                Ai(4:6,obj.vars.weight1.i(4,c,l,t)) = obj.object.lines{f}.fc4(:,t);
                
                Ai(:,obj.vars.L.i(f,c,l,t)) = M;

                obj = obj.addLinearConstraints(Ai, bi, [], []);

                Ai = sparse(6,obj.nv);
                bi = M*ones(6,1);

                Ai(1:3,obj.vars.f.i(:,c,l,t)) = eye(3);
                Ai(1:3,obj.vars.alpha_f.i(1,c,l,t)) = -(obj.object.lines{f}.fc1(:,t)+obj.object.lines{f}.fc2(:,t)+obj.object.lines{f}.fc3(:,t)+obj.object.lines{f}.fc4(:,t))/4;
                Ai(1:3,obj.vars.weight.i(1,c,l,t)) = -obj.object.lines{f}.fc1(:,t);
                Ai(1:3,obj.vars.weight.i(2,c,l,t)) = -obj.object.lines{f}.fc2(:,t);
                Ai(1:3,obj.vars.weight.i(3,c,l,t)) = -obj.object.lines{f}.fc3(:,t);
                Ai(1:3,obj.vars.weight.i(4,c,l,t)) = -obj.object.lines{f}.fc4(:,t);
                
                Ai(4:6,obj.vars.f.i(:,c,l,t)) = -eye(3);
                Ai(4:6,obj.vars.alpha_f.i(1,c,l,t)) = (obj.object.lines{f}.fc1(:,t)+obj.object.lines{f}.fc2(:,t)+obj.object.lines{f}.fc3(:,t)+obj.object.lines{f}.fc4(:,t))/4;
                Ai(4:6,obj.vars.weight.i(1,c,l,t)) = obj.object.lines{f}.fc1(:,t);
                Ai(4:6,obj.vars.weight.i(2,c,l,t)) = obj.object.lines{f}.fc2(:,t);
                Ai(4:6,obj.vars.weight.i(3,c,l,t)) = obj.object.lines{f}.fc3(:,t);
                Ai(4:6,obj.vars.weight.i(4,c,l,t)) = obj.object.lines{f}.fc4(:,t);
                
                Ai(:,obj.vars.L.i(f,c,l,t)) = M;

                obj = obj.addLinearConstraints(Ai, bi, [], []);
              else 
                % TODO in 3D
              end
            end
          end 
        end 
      end
    end

    function obj = addNonPenetrationConstraints(obj)
      % Constrains each implicit contact point to lie in a convex region of 
      % space covering the the complement of the object.

      % region assignment constraints
      obj = obj.addVariable('R', 'B', [obj.N_r, obj.N_c, obj.N_l, obj.N_T], 0, 1);
      regions = obj.object.regions;
      K = 100;

      env_regions = obj.object.env_regions;

      for t = 1:obj.N_T
        for l = 1:obj.N_l
          for c = 1:obj.N_c
            % non-floor penetration
            for r = 1:length(env_regions)
              Ai = sparse(size(env_regions{r}.b,1),obj.nv);
              bi = env_regions{r}.b;
              
              Ai(:,obj.vars.p.i(:,c,l,t)) = env_regions{r}.A;
              
              obj = obj.addLinearConstraints(Ai, bi,[], []);
            end
            % assignment constraints
            for r = 1:obj.N_r
              % loads the matrices
              A = regions{r}.A(:,:,t);
              b = regions{r}.b(:,t);            

              Ai = sparse(size(b,1),obj.nv);
              bi = b + K;
              
              Ai(:,obj.vars.p.i(:,c,l,t)) = A;
              Ai(:,obj.vars.R.i(r,c,l,t)) = K;
              
              obj = obj.addLinearConstraints(Ai, bi,[], []);

              if t < obj.N_T
                Ai = sparse(size(b,1),obj.nv);
                bi = regions{r}.b(:,t+1) + K;
                
                Ai(:,obj.vars.p.i(:,c,l,t+1)) = regions{r}.A(:,:,t+1);
                Ai(:,obj.vars.R.i(r,c,l,t)) = K;
                
                % obj = obj.addLinearConstraints(Ai, bi,[], []);
              end

              if c < obj.N_c
                Ai = sparse(size(b,1),obj.nv);
                bi = b + K;
                
                Ai(:,obj.vars.p.i(:,c,l,t)) = A;
                Ai(:,obj.vars.R.i(r,c+1,l,t)) = K;
                
                % obj = obj.addLinearConstraints(Ai, bi,[], []);
              end
            end

            % each point must be in a region
            Aeq = sparse(1,obj.nv);
            beq = 1;
            
            Aeq(1,obj.vars.R.i(:,c,l,t)) = 1;

            obj = obj.addLinearConstraints([], [], Aeq, beq);
          end
        end
      end
    end

    function obj = addCostFunction(obj)
      % minimizes the following metrics:
      % 1. Acceleration of the contacts
      % 2. Mode switches
      % TODO: maximize stability of grasp
      for t = 1:obj.N_T
        % min_marging
        qi = zeros(1,obj.nv);
        qi(1,obj.vars.alpha_min.i(:,t)) = -100;
        obj = obj.addCost([],qi',[]);
        
        for v = 1:obj.object.nv 
          Qi = sparse(obj.nv,obj.nv);
          Qi(obj.vars.f_ext.i(:,v,t),obj.vars.f_ext.i(:,v,t)) = eye(3);
          % obj = obj.addCost(Qi,[],[]);
        end
        for c = 1:obj.N_c
          for l = 1:obj.N_l  
            Qi = sparse(obj.nv,obj.nv);
            Qi(obj.vars.ddp.i(:,c,l,t),obj.vars.ddp.i(:,c,l,t)) = eye(3);
            obj = obj.addCost(Qi,[],[]);

            Qi = sparse(obj.nv,obj.nv);
            Qi(obj.vars.f.i(:,c,l,t),obj.vars.f.i(:,c,l,t)) = eye(3);
            obj = obj.addCost(Qi,[],[]);
          end
        end
      end
    end

    % end of methods
  end
end