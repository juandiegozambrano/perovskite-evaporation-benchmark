%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: MPCalgorithm.m
%
% Description:
%   This script implements a Model Predictive Control (MPC) strategy 
%   formulated as a Quadratic Program (QP) and solved using MATLAB's 
%   quadprog solver. The objective is to compute, at each iteration, 
%   the optimal power input sequence that minimizes the future tracking 
%   error of the deposition rate under actuator constraints.
%
% MPC formulation:
%   - Decision vector:   Z = [x(0..N); e(0..N); u(0..N-1)]
%   - Optimization:      min Z   (1/2)*Z'*H*Z
%                        subject to   Aeq*Z = beq (model + error dynamics)
%                                      A*Z <= b     (input bounds)
%   - System dynamics:   x(k+1) = AG*x(k) + BG*u(k)
%   - Error definition:  e(k) = r_ref – r(k), where r(k) = CG(2,:)*x(k)
%
% Requirements:
%   - MATLAB Optimization Toolbox (for quadprog)
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

clear all; close all; clc

%% MPC Problem Setup 

% Deposition rate reference
r_ref = 0.5; 

% Optimization solver options
options = optimset('Display','off', 'TolFun', 1e-8,'MaxIter', 10000,...
                    'TolConSQP', 1e-6);

% MPC parameters
mpciterations = 2000;   % Simulation time in seconds, since Ts=1s
T = 120;                % Prediction horizon [seconds]
Ts = 1;                 % Sampling time
N = T/Ts;               % Prediction horizon [discrete time steps]

% Control-oriented model considers only rate as output
load("ModelOL_G.mat") % Load the model G
Ad = AG;        % System matrix 
Bd = BG;        % Input matrix 
Cd = CG(1,:);   % Output matrix (only rate is used, i.e., C2)    
Dd = 0;

n = length(Ad(1,:)); % Number of states in x(k)
m = length(Bd(1,:)); % Number of inputs in u(k)

% Optional terminal state (disabled here)
xTerm = zeros(4,1);

% Cost weights
Q = 0*eye(n);           % Weight on original state x(k) (not used here)
Qe = 100*eye(N+1);      % Weight on error e(k) over horizon
R = 0.0000001*eye(m);   % Weight on control input effort

% Initial conditions
tm = 0;                 % Time counter 
ym = 0;                 % Measured output 
xm = initial_state_G';  % Initial state of model G at time k_m. That is, considering post-threshold dynamics 
                        % [-30599.445570938115 -0 0 0.11696176707252705]'; 

% Initial guesses of optimization variables 
u0 = repmat(zeros(m,1),N,1);    % input sequence
x0 = repmat(xm,N+1,1);          % state sequence (order: [x_1(0);x_2(0);x_1(1);x_2(1);...])

% Storage for logging
t = []; x = []; u = []; y = []; setp = [];


%%  Constraints Setup 

%%% Power input bounds: 0 <= u(k) <= 60 
% Inequality constraint in the form  H_u * u(k) <= k_u
H_u = [1; -1];
k_u = [60; 0];
A1 = kron(eye(N), H_u);
A = [zeros(length(A1(:,1)),(n+1)*(N+1)), A1];  % Pad with zeros to align in Z = [x,e,u]
b = repmat(k_u,N,1);

s1 = 0; % Apply terminal constraint if s1 = 1 

%%% State dynamics and error 
% Equality constraint matrices in the form Aeq * [x; e] = beq

% Preallocate matrices Aeq1 and beq1
Aeq1 = zeros(n*(N+1), n*(N+1)+m*N);  
beq = zeros((n+1)*(N+1+s1), 1);     

% State dynamic constraints x(k+1) = Ad*x(k) + Bd*u(k)  for k = 0,...,N-1
for k = 0:N-1
    Aeq1(n*k+1:n*(k+1),1:n*(N+1)) = [zeros(n,n*k), Ad, -eye(n), zeros(n,n*(N-1-k))];
    Aeq1(n*k+1:n*(k+1),n*(N+1)+1:end) = [zeros(n,m*k), Bd, zeros(n,m*(N-1-k))];
end
clear('k')

% Now expand to include the error variables e(k) and enforce e(k) = r_ref - Cd*x(k)
Aeq2 = [Aeq1(:,1:n*(N+1)), zeros(length(Aeq1(:,1)),N+1), Aeq1(:,1+n*(N+1):end)];
Aeq = [Aeq2; [kron(eye(N+1),Cd), eye(N+1), zeros(N+1,N)]]; 

l = length(Aeq(:,1));
ll = length(Aeq(1,:)); 

%%% Set initial conditions: x(0) = xm, and e(0) = r_ref
% Equality constraint 
Aeq(1+l:n+l,1:ll) = [eye(n),zeros(n,ll-n)];
beq(1+l:n+l) = xm;                             
beq(1+n*(N+1):(n+1)*(N+1)) = r_ref*ones(N+1,1); 

%%% Optional terminal constraint: X(N) = xTerm
if s1 == 1
    v = 4;
    l = length(Aeq(:,1));
    Aeq(l+1:l+n,:) = [zeros(n,n*N),eye(n),zeros(n,(m*N)+(N+1))];
    beq(l+1:l+n) = xTerm;
else
    v = 0;
end

%% Cost funtion

% Build quadratic weighting matrices across horizon
Qstack = []; Rstack = [];
for k = 1:N
    Qstack = blkdiag(Qstack,Ts*Q);  
    Rstack = blkdiag(Rstack,Ts*R);   
end
clear('k')

% Terminal state cost is zero
Qstack = blkdiag(Qstack,zeros(n));

% Assemble overall quadratic matrix H and linear term f
% H corresponds to blockdiag (Qstack over x, Qe over e, Rstack over u)
H = blkdiag(Qstack,Qe,Rstack);
f = zeros(length(Aeq(1,:)),1); % No linear term

%% MPC simulation 

% Print Header
fprintf('   k  |      u(k)        x(1)        x(2)        x(3)        x(4)     Time \n');
fprintf('---------------------------------------------------------------------------\n');
l = length(Aeq(:,1))-n;
u_OL(1) = 60;

for ii = 1:mpciterations

    % Update initial condition in constraint (based on current xm)
    beq((1+l-v):n-v+l) = xm;
    t_Start = tic;

    % Solve QP optimization problem
    % Z* = argmin_Z (1/2*Z'HZ) 
    % subject to: Aeq*Z = beq, A*Z <= b
    Z_opti_vector = quadprog(H,f,A,b,Aeq,beq,[],[],[],options);
    
    % Extract predicted trajectories from solution vector
    x_OL_tilde = Z_opti_vector(1:n*(N+1),1);                % Extract state 
    e_OL_tilde = Z_opti_vector(1+n*(N+1):(n+1)*(N+1),1);    % Extract error 
    u_OL_tilde = Z_opti_vector((n+1)*(N+1)+1:end,1);        % Extract control action 
    x_OL = reshape(x_OL_tilde, n, N+1);
    u_OL = u_OL_tilde';
    t_Elapsed = toc(t_Start);

    % Save current results 
    t = [t, tm];
    x = [x, xm];
    u = [u, u_OL(1)];
    y = [y, ym];
    setp = [setp, r_ref];

    % Apply first control move to the system
    xm = Ad*xm + Bd*u_OL(1:m);  % state evolution: xm = Ad*xm + Bd*u
    ym = Cd*xm;                 % rate dynamics: r = Cd*x
    tm = tm + Ts;

    % Warm-start next iteration with shifted sequence
    x0 = [x_OL_tilde(n+1:end); zeros(n,1)];
    u0 = [u_OL_tilde(m+1:end); zeros(m,1)];

    % Debug print
    fprintf(' %3d  | %+11.6f %+11.6f %+11.6f %+11.6f %+11.6f %+6.3f\n', ii, u(end),...
        x(1,end), x(2,end),x(3,end), x(4,end),t_Elapsed);
end

% Save results
save('MPCresults.mat','t','x','u','y','setp');


%% Generate graphical results

f = figure;
set(f, 'Color', 'w');

% Subplot 1: System Response
subplot(2,1,1);
plot(t, setp, 'k--', 'LineWidth', 2); hold on;
plot(t, y, 'LineWidth', 2);

title('MPC Rate Control', ...
      'Interpreter', 'latex', 'FontSize', 16)
xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('Rate (\AA/s)', 'Interpreter', 'latex', 'FontSize', 14)

legend({'Setpoint $r^\mathrm{ref}$', 'Rate $r$'}, ...
       'Interpreter', 'latex', 'FontSize', 14, 'Location', 'best')

ylim([0 0.6]);         % Set y-axis limits             
%xlim([min(t) max(t)]); % Fit x-axis to data range
grid on; box on;

set(gca, 'FontSize', 13, ...
         'LineWidth', 1.2, ...
         'TickLabelInterpreter', 'latex');

% Subplot 2: Control Input
subplot(2,1,2);
stairs(t, u, 'LineWidth', 2);

title('MPC Power Input', ...
      'Interpreter', 'latex', 'FontSize', 16)
xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('Power input (\%) ', 'Interpreter', 'latex', 'FontSize', 14)

grid on; box on;
set(gca, 'FontSize', 13, ...
         'LineWidth', 1.2, ...
         'TickLabelInterpreter', 'latex');


% Save figure
disp("Figure saved")
savefig(f, 'MPC_RateControl.fig');           
