%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: PRGalgorithm.m
%
% Description: This script solves the Predictive Reference Governor (PRG)
% problem using CasADI with the IPOPT solver. By running this file,
% the optimal virtual reference is obtained and applied to the closed-loop
% model.
%
% Requirements:
%   - CasADI toolbox for MATLAB: https://web.casadi.org/get
%
% Repository: https://github.com/juandiegozambrano/perovskite-evaporation-benchmark
% Version: 1.0
% Date: 02-10-2025
% Author: J.D. Zambrano-Torres
%
% If you use this script or its data, please cite:
% E. Masero, J.D. Zambrano-Torres, J. Vollbrecht, J.M. Maestre (2026).
% "A Benchmark on Perovskite Thin-Film Deposition via Thermal Evaporation
% for Photovoltaic Solar Cell Manufacturing Systems." https://doi.org/xxxxxxx
%
% License: MIT License
% Contact: spjuandiego@outlook.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

addpath(fullfile(pwd,'CASADIfiles'))    % Add folder containing CasADI functions to MATLAB path
import casadi.*                         % Import CasADi package

% State space representation
load('ModelCL_H.mat'); % Load the closed-loop model and initial state
A = modelH_ss.A;       % State-space A matrix
B = modelH_ss.B;       % State-space B matrix
C = modelH_ss.C;       % State-space C matrix
D = modelH_ss.D;       % State-space D matrix

nx = size(A,1);        % Number of states
nu = size(B,2);        % Number of inputs

% Constraints and parameters 
w_max = 0.65;         % Maximum virtual reference value
w_min = 0;            % Minimum virtual reference value
pmx   = 60;           % Maximum converter output (power constraint)

Ts = 1;               % Sampling time [s]
N_pred = 120;         % Prediction horizon [discrete time steps]
Qry = 1300;           % Tracking weight
Qrw = 25;             % Convergence weight
Qwy = 10;             % Shape weight
Qdw = 15;             % Damping weight
N_sim = 2000;         % Total simulation steps

ref = 0.5;                      % Desired rate reference
x = initial_state_closedLoop';  % Initial state vector

% History vectors to plot 
x3_history = zeros(N_sim+1, 1);  % History of power (x3)
c1_history = zeros(N_sim+1, 1);  % History of output rate
u_history  = zeros(N_sim, 1);    % History of virtual reference (control input)

% Initial values 
x3_history(1) = x(3);        
c1_history(1) = C(1,:)*x;   

% Define the optimization problem (created once) 
opti = casadi.Opti();     % Instantiate an optimization problem

% Decision variables
X = opti.variable(nx, N_pred+1);  % Future predicted states
U = opti.variable(nu, N_pred);    % Future predicted inputs (virtual reference w)

% Parameter for current state
x_param = opti.parameter(nx, 1);   % Parameter to update current state at each iteration

% Initial condition constraint
opti.subject_to(X(:,1) == x_param);  % Enforce equality between initial predicted and actual state

% System dynamics
for j = 1:N_pred             
    opti.subject_to(X(:, j+1) == A*X(:, j) + B*U(:, j));            
end

% Inequality constraints
opti.subject_to(0 <= X(3,:) <= pmx);        % Power constraint (0 <= x3 <= 60)
opti.subject_to(w_min <= U(1,:) <= w_max);  % Virtual reference limits (0 <= U <= 0.65)

% Cost function definition
cost = 0;                                                          
for j = 1:N_pred-1              
    outs = C(1,:)*X(:,j);     % System output (rate)
    outs2 = U(:,j);           % Virtual reference w(n|k)
    outs3 = U(:,j+1);         % Virtual reference w(n+1|k)

    % Quadratic cost components: tracking, convergence, shaping, and damping
    cost = cost + Qry*(ref - outs)^2    + Qrw*(ref - outs2)^2 ...
                + Qwy*(outs2 - outs)^2  + Qdw*(outs3 - outs2)^2;
end
opti.minimize(cost);    % Define objective for minimization

% Solver configuration with warm start
p_opts = struct('print_time', 0,'verbose',false);                   % Problem-level options (silent mode)
s_opts = struct('print_level',0,'warm_start_init_point','yes');     % Solver-level options (enable warm start)
opti.solver('ipopt', p_opts, s_opts);                               % Assign IPOPT solver to problem

%% PRG simulation 
x = x(:);                               % Ensure column vector format
last_sol_U = zeros(nu, N_pred);         % Initialize previous optimal input trajectory (for warm start)
last_sol_X = repmat(x, 1, N_pred+1);    % Initialize previous state trajectory (for warm start)

for k = 1:N_sim                                                    
    opti.set_value(x_param, x);         % Update the parameter with the current state

    % Warm start from previous solution
    opti.set_initial(U, last_sol_U);    % Set previous solution as initial guess for inputs
    opti.set_initial(X, last_sol_X);    % Set previous solution as initial guess for states

   
    % Solve optimization problem
    try                                  % Extract optimal control and update warm start data
        sol = opti.solve();              
        u_opt       = sol.value(U(:,1)); % Extract optimal optimal virtual reference w*(k|k) 
        last_sol_U  = sol.value(U);      % Save full solution U for warm start
        last_sol_X  = sol.value(X);      % Save full solution X for warm start
    catch
        u_opt = last_sol_U(:,1);           % Reuse last valid control

    end

    % Soft start: fix virtual reference at max for first 10 steps
    if k <= 10               % Apply soft start to avoid transients
        u_opt = w_max;       % Use maximum virtual reference
    end

    % Update system state
    x = A*x + B*u_opt;                

    % Store results
    u_history(k)    = u_opt;     % Save applied input (w)
    x3_history(k+1) = x(3);      % Save updated power state (P)
    c1_history(k+1) = C(1,:)*x;  % Save updated output rate (r)
end

x3_history(1:10) = pmx;   % Adjust initial display values for consistency

% Save results
save('PRGresults.mat','u_history','x3_history','c1_history');       

%% Plot results 

N_sim = 2000;       % Total simulation steps
t  = 0:N_sim;       % Time vector for states/outputs
ts = 0:N_sim-1;     % Time vector for control inputs

%%% Power input and deposition rate
figure;                                                  

% Subplot 1: output rate, r(k)
subplot(3,1,1)                                                      
plot(t, c1_history, 'b-', 'LineWidth', 2); hold on;   
yline(ref, 'k--', 'LineWidth', 2); 
ylim([0 0.6]);         % Set y-axis limits  
xlabel('Time (s)','Interpreter', 'latex'); 
ylabel('Rate (\AA/s)','Interpreter', 'latex'); grid on;    
title('PRG Rate Control');   

% Subplot 2: power input, P(k)
subplot(3,1,2)                                                     
plot(t, x3_history, 'b-', 'LineWidth', 2); hold on;                 
yline(60, 'k--','LineWidth', 2);                                           
yline(0, 'k--','LineWidth', 2);                                            
xlabel('Time (s)','Interpreter', 'latex'); 
ylabel('Power input (\%)','Interpreter', 'latex');
grid on;                  
title('PRG Power Input');                              

% Subplot 3: virtual reference, w(k)
subplot(3,1,3)                                                      
plot(ts, u_history, 'b-', 'LineWidth', 2); hold on;                 
yline(w_max, 'k--', 'LineWidth', 2);                                         
yline(w_min, 'k--', 'LineWidth', 2);                                         
xlabel('Time (s)','Interpreter', 'latex'); 
ylabel('Virtual Ref. (\AA/s)','Interpreter', 'latex'); grid on;     
title('PRG Virtual Reference');                               

%% Virtual reference plot 
t = 1:1:2000;
a1 = u_history(1:end);

figure;
hold on;
cRG  = [0.00, 0.45, 0.70]; % color azul

p1 = plot(t, u_history, 'Color',cRG, 'LineWidth', 2);      % Virtual reference
p2 = plot(t, 0.5 * ones(size(t)), 'k--', 'LineWidth', 2);  % Global reference

xlabel('Time (s)', 'Interpreter', 'latex');
ylabel('Rate~$\mathrm{(\AA/s)}$', 'Interpreter', 'latex');
ylim([0.3 0.65]);

grid on;
legend({'Virtual reference, $w$', 'Global reference, $r^\mathrm{ref}$'}, ...
       'Interpreter', 'latex', 'Location', 'northeast');
set(gca, 'FontSize', 16);
set(gca, 'Box', 'on');
set(gcf, 'Color', 'w');

disp("Figure saved")
savefig(gcf, 'PRG_virtual_ref.fig');  % Save in .fig