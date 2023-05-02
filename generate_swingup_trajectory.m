function [T, U_sol, X_sol, U_guess, X_guess] = generate_swingup_trajectory(p)
    % Generate an optimal open-loop trajectory for the swing-up problem
    % for the cartpole system described by input parameter struct.
    % Inputs:
    %    p is parameter struct specifying scalar system parameters:
    %        This struct must contain the following scalar values:
    %        num_sample (=n), time_final, pos_final, pos_min, pos_max,
    %        force_min, force_max, mass_cart, mass_pole, len_pole,
    %        grav_accel.
    % Outputs:
    %    T is [n, 1] time vector of sampled points in time.
    %    U_sol is [n, 1] vector of control effort exerted at each time in T.
    %    X_sol is [n, 4] matrix of the 4 state variables at each time in T:
    %        Its columns are linear position of cart, angular position of pole,
    %        linear velocity of cart, and angular velocity of pole respectively.
    %    U_guess is [n, 1] vector of initial guess of control effort over time.
    %    X_guess is [n, 4] matrix of initial guess of state state over time:
    %        It has the same structure as X_sol.


    % ---- UNPACK SYSTEM PARAMETERS ----

    n  = p.num_samples;
    m1 = p.mass_cart;
    m2 = p.mass_pole;
    L  = p.len_pole;
    g  = p.grav_accel;
    d  = p.pos_final;
    
    % Compute time vector and timestep from system parameters 
    T  = linspace(0, p.time_final, n)';
    dt = T(2) - T(1);


    % --- SPECIFY STATE BOUNDARY CONDITIONS ----

    % These initial and final states are inherent to this swingup problem.
    % Cart starts at linear position 0 and ends at d specified in parameters.
    % Pendulum starts hanging with angle=0 and ends upright at (pi) radians.
    % System starts & ends at rest, so initial and final velocities are zero.

    %  state:      lin_pos   ang_pos   lin_vel   ang_vel
    %  units:      [m]       [rad]     [m/s]     [rad/s]
    state_init  = [ 0;        0;        0;        0 ];
    state_final = [ d;        pi;       0;        0 ];


    % ---- LOWER & UPPER BOUNDS ON DECISION VARS ----

    % Constrain decision vector with simple inequalities between a minimum
    % and maximum allowable value. We use these to constrain the cart to 
    % the length of the track and to enforce our limits on control effort.
    % We also use these bounds to "sandwich" elements to be exactly equal
    % to a particular value. Specifically, we constrain the system's initial
    % and final states to those specified above!

    % Control (force) values must simply respect max force bounds
    max_control = p.force_max * ones(n,1);
    min_control = p.force_min * ones(n,1);

    % Linear position: Cart must stay on track & start+end at specified positions
    max_lin_pos = [state_init(1);  p.pos_max*ones(n-2,1);  state_final(1)];
    min_lin_pos = [state_init(1);  p.pos_min*ones(n-2,1);  state_final(1)];

    % All other state variables must simply start+end at specified values
    max_ang_pos = [state_init(2);   Inf(n-2,1);  state_final(2)];
    min_ang_pos = [state_init(2);  -Inf(n-2,1);  state_final(2)];

    max_lin_vel = [state_init(3);   Inf(n-2,1);  state_final(3)];
    min_lin_vel = [state_init(3);  -Inf(n-2,1);  state_final(3)];

    max_ang_vel = [state_init(4);   Inf(n-2,1);  state_final(4)];
    min_ang_vel = [state_init(4);  -Inf(n-2,1);  state_final(4)];

    % These bounding vectors get concatenated into full upper & lower bound
    % vectors which have the same structure as the full decision vector.
    % decision:        control;     lin_pos;     ang_pos;     lin_vel;     ang_vel
    upper_bound = [max_control; max_lin_pos; max_ang_pos; max_lin_vel; max_ang_vel];
    lower_bound = [min_control; min_lin_pos; min_ang_pos; min_lin_vel; min_ang_vel];
    

    % ---- INITIALIZE DECISION VECTOR ----

    % Initialize control and state vectors with no control effort and
    % linear interpolation between initial+final linear & angular positions
    control_guess = zeros(n, 1);
    lin_pos_guess = linspace(state_init(1), state_final(1), n)';
    ang_pos_guess = linspace(state_init(2), state_final(2), n)';
    lin_vel_guess = zeros(n, 1);
    ang_vel_guess = zeros(n, 1);

    % Concatenate guess vectors into full decision vector guess
    guess = [ control_guess; lin_pos_guess; ang_pos_guess;
              lin_vel_guess; ang_vel_guess ];


    % ======== RUN OPTIMIZATION ========

    % Set fmincon to run with interior point "FeasibilityMode", which
    % prefers maintaining feasibility over obtaining strict optimality
    % in objective function. This is important in our case because
    % respecting the dynamics of the system is most crucial for
    % achieving a trajectory which works in practice.

    options = optimoptions('fmincon', ...
              'Algorithm',             'interior-point', ...
              'EnableFeasibilityMode',  true, ...
              'SubproblemAlgorithm',   'cg', ...
              'MaxFunctionEvaluations', 1e6);

    % Run fmincon and hopefully obtain solution!

    SOLUTION = fmincon(@cost_function, guess, [], [], [], [], ...
                lower_bound, upper_bound, @nonlinear_dynamics, options);


    % ---- FORMAT & RETURN OUTPUTS ----

    % Unpack solution vector into control vector & state matrix.
    U_sol = SOLUTION(1:n);
    X_sol = reshape( SOLUTION(n+1:end), [n,4] );

    % Return initial solution guesses for comparison to actual solution.
    U_guess = control_guess;
    X_guess = reshape( guess(n+1:end), [n,4] );




    function cost = cost_function(decision_vec)
        % Compute cost of the candidate solution encoded in the decision vector.
        % Inputs:
        %    decision_vec is [n*5, 1] vector of candidate optimization solution.
        % Outputs:
        %    cost is scalar cost of solution, which fmincon seeks to minimize.

        % Extract control sequence from decision vector
        control_vec = decision_vec(1:n);
        
        % Take trapezoidal sum over elementwise-squared control vector
        % representing control effort exerted at each sampled point.
        % This approximates the integral of the squared control effort in
        % the continuous case. We can safely ignore the duration of the
        % timestep, as this only scales the cost by a constant.

        cost = trapz(control_vec.^2);
    end


    function [c,ceq] = nonlinear_dynamics(decision_vec)
        % Compute nonlinear functions of decision vector which solver will
        % attempt to constrain. The solver tries to satisfy ceq=0 and c<=0.
        % We use the equality constraint on ceq to constrain the state and
        % control at our sampled points such that they respect the nonlinear
        % dynamics of our physical cartpole system.
        % This function's interface is dictated by fmincon's <nonlcon> method.
        % Inputs:
        %    decision_vec is [n*5, 1] vector of candidate optimization solution.
        % Outputs:
        %    c is null vector (as required by fmincon).
        %    ceq is [(n-1)*4, 1] vector which fmincon will constrain to 0.
        
        c = [];    % We do not use c, which expresses nonlinear inequalities

        % The output ceq will have one expression per state variable per sample
        % interval. Preallocate space for (4) state variables * (n-1) intervals
        ceq = zeros(4*(n-1), 1);
        
        for sample = 1:n
            % Extract current state & control from decision vector
            control   = decision_vec(sample);
            state     = decision_vec(sample+[n, n*2, n*3, n*4]);

            % Take time derivative of state
            state_dot = differentiate_state(state, control);

            % First interval gets calculated at second sample point
            if sample >= 2

                % Compute index for allocation in ceq vector
                ind = (sample-2)*4 + 1;
                
                % At the most basic level, our dynamic constraints can be
                % expressed as Euler's method: the next state must be equal
                % to the current state plus the derivative times the timstep.
                % Since we assume the dynamics are linear between sample points,
                % each sample has a sharp corner, so we approximate the derivative
                % at this point as the average of the derivatives before & after
                % each sample point. Thus our constraint is:
                %    state = state_prev + dt*(state_dot + state_dot_prev)/2
                % Since fmincon requires we express this constraint via a
                % vector that it will drive to zero, we simply subtract the
                % right hand side from both sides of the equation so it equals 0.

                ceq(ind:ind+3) = state - state_prev ...
                                   - dt*(state_dot + state_dot_prev)/2;
            end

            % Store state & derivative for use in next interval calculation
            state_prev     = state;
            state_dot_prev = state_dot;
        end
    end


    function state_dot = differentiate_state(state, control)
        % Compute derivative of cartpole system state with respect to time
        % as a function of current state & control.
        % Inputs:
        %    state is [4,1] vector of cartpole system state
        %    control is scalar value of control effort
        % Outputs:
        %    state_dot is [4,1] vector of state derivative

        % Unpack state variables from state vector
        % lin_pos = state(1);   % Linear position is unused in dynamics
        ang_pos = state(2);
        lin_vel = state(3);
        ang_vel = state(4);

        % System dynamics govern accelerations as function of state and
        % control. These equations due to [Kelly], [Tedrake]

        lin_accel =  ( L*m2*sin(ang_pos)*ang_vel^2 ...
                       + control + m2*g*cos(ang_pos)*sin(ang_pos) ) ...
                      / ( m1 + m2*sin(ang_pos)^2 );

        ang_accel = -( L*m2*cos(ang_pos)*sin(ang_pos)*ang_vel^2 ...
                       + control*cos(ang_pos) + (m1+m2)*g*sin(ang_pos) ) ...
                      / ( L*m1 + L*m2*sin(ang_pos)^2 );

        % Pack derivatives of state variables into state derivative vector
        %  state:    lin_pos;  ang_pos;   lin_vel;    ang_vel
        state_dot = [lin_vel;  ang_vel;  lin_accel;  ang_accel];
    end

end