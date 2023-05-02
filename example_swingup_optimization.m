% example_swingup_optimization.m
% Luke Raus 2023
%
% This script demonstrates setting up the parameters of a cartpole system,
% generating an optimal trajectory using generate_swingup_trajectory(),
% and comparing this approximated trajectory against an external simulation
% using the trajectory's open-loop control sequence.


% ---- SET SYSTEM PARAMETERS ----

params = struct;

params.num_samples = 50;

params.time_final = 2;      % [sec]
params.pos_final  = 1;      % [m]

params.pos_min    = -2;     % [m]
params.pos_max    =  2;     % [m]

params.force_min  = -20;    % [N]
params.force_max  =  20;    % [N]

params.mass_cart  = 1;      % [kg]
params.mass_pole  = 0.3;    % [kg]
params.len_pole   = 0.5;    % [m]
params.grav_accel = 9.81;   % [m/s^2]


% ---- RUN OPTIMIZATION ----

[T, U_sol, X_sol, U_guess, X_guess] = generate_swingup_trajectory(params);

% ---- PLOT OPTIMIZATION RESULTS ----

figure
tiledlayout(3,1)
title('Optimal swing-up trajectory generated')
nexttile;  plot(T, U_sol, '.-');         xlabel('Time [sec]'); ylabel('Control force [N]')
nexttile;  plot(T, X_sol(:,1),   'm.-'); xlabel('Time [sec]'); ylabel('Linear pos [m]')
 hold on;  plot(T, X_guess(:,1), 'm:');  hold off
 legend('Optimized', 'Initial guess')
nexttile;  plot(T, X_sol(:,2),   'r.-'); xlabel('Time [sec]'); ylabel('Angular pos [rad]')
 hold on;  plot(T, X_guess(:,2), 'r:');  hold off


% ---- RUN EXTERNAL SIMULATION ----

[T_sim, X_sim] = simulate_cartpole(params, U_sol);

% ---- PLOT SIMULATION RESULTS ----

figure
tiledlayout(2,1)
title('Comparison of trajectory generated with external simulation')
nexttile;  plot(T_sim, X_sim(:,1), 'm');  xlabel('Simulation time [sec]'); ylabel('Linear pos [m]')
 hold on;  plot(T, X_sol(:,1),     'm:'); hold off
 legend('Simulated', 'Optimized')
nexttile;  plot(T_sim, X_sim(:,2), 'r');  xlabel('Simulation time [sec]'); ylabel('Angular pos [rad]')
 hold on;  plot(T, X_sol(:,2),     'r:'); hold off
