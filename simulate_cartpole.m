function [T_sim, X_sim] = simulate_cartpole(param, control_seq)
    % Forward-simulate the cartpole system for a given open-loop control
    % sequence using ode45.
    % Inputs:
    %    param is parameter struct specifying scalar system parameters:
    %        This struct must contain the following scalar values:
    %        num_sample (=n), time_final, mass_cart, mass_pole,
    %        len_pole, grav_accel.
    %    control_seq is [n, 1] vector specifying control value at each
    %        sampled point in time, as determined in param.
    % Outputs:
    %    T_sim is vector of time values for resulting simulation
    %    X_sim is matrix of simulated state values over time

    % Unpack parameters from struct
    m1 = param.mass_cart;
    m2 = param.mass_pole;
    L  = param.len_pole;
    g  = param.grav_accel;

    % Prepare time vector
    n  = param.num_samples;
    T  = linspace(0, param.time_final, n)';

    % Simulate for extra 20% of time to see if system stabilizes
    t_span = [0, T(end)*1.2];
    
    % Initial system state is at rest
    state_init = [0; 0; 0; 0];

    % Set ode45 to low error tolerance for precise simulation
    options = odeset('RelTol',1e-6);

    % Run ode45 simulation!
    [T_sim, X_sim] = ode45(@rate_func, t_span, state_init, options);


    function state_dot = rate_func(t, state)
        % Compute derivative of cartpole system state with respect to time
        % as a function of current state, as the rate_func used by ode45.
        % Control input to the derivative is interpolated from control_seq.
        % Inputs:
        %    state is [4,1] vector of cartpole system state
        % Outputs:
        %    state_dot is [4,1] vector of state derivative

        % Get control at this moment in time by linearly interpolating along
        % control vector, if in the control timespan. Otherwise no control.
        if t <= T(end)
            control = 1.0*interp1(T, control_seq, t);
        else
            control = 0;
        end

        %lin_pos = state(1);   % Linear position is unused in dynamics
        ang_pos = state(2);
        lin_vel = state(3);
        ang_vel = state(4);

        % System dynamics govern accelerations as function of state
        % These equations due to [Kelly], [Tedrake]

        lin_accel =  ( L*m2*sin(ang_pos)*ang_vel^2 ...
                       + control + m2*g*cos(ang_pos)*sin(ang_pos) ) ...
                      / ( m1 + m2*sin(ang_pos)^2 );

        ang_accel = -( L*m2*cos(ang_pos)*sin(ang_pos)*ang_vel^2 ...
                       + control*cos(ang_pos) + (m1+m2)*g*sin(ang_pos) ) ...
                      / ( L*m1 + L*m2*sin(ang_pos)^2 );

        % Pack derivatives of state variables into derivative state vector
        %  state  = [lin_pos; ang_pos; lin_vel;   ang_vel]
        state_dot = [lin_vel; ang_vel; lin_accel; ang_accel];
    end

end
