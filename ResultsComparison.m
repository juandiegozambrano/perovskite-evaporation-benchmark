%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: ResultsComparison.m
%
% Description: This script generates graphical results and performance 
% metrics ISE (Integral Square Error) and ISCO (Integral Square Control
% Output) for the three controllers:
%
%   - PID (From dataset 2)
%   - MPC for rate control
%   - PRG for rate control
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

%% Load data
load('PIDresults.mat');  % Variables: tpid, rpid, upid
load('MPCresults.mat');  % Variables: t, y, u
load('PRGresults.mat');  % Variables: c1_history, u_history, x3_history

% Define time vectors 
time = 1:2000;          
          
% Define tracking error 
ref = 0.5;                              % Rate setpoint
errorPID = ref - rpid;                  % PID rate error
errorMPC = ref - y;                     % MPC rate error
errorPRG = ref - c1_history(1:2000);    % PRG rate error

% Define power input vectors
uPID = upid;                 
uMPC = u;                    
uPRG = x3_history(1:2000);  

%% Performance metrics ISE and ISCO  
ISEa = trapz(t, errorPID.^2);
ISEb = trapz(t, errorMPC.^2);
ISEc = trapz(t, errorPRG.^2);

ISCOa = trapz(t, uPID.^2);
ISCOb = trapz(t, uMPC.^2);
ISCOc = trapz(t, uPRG.^2);

fprintf('\n--- PERFORMANCE INDICES ---\n');
fprintf('PID  -> ISE = %.5f, ISCO = %.5f\n', ISEa, ISCOa);
fprintf('MPC  -> ISE = %.5f, ISCO = %.5f\n', ISEb, ISCOb);
fprintf('PRG  -> ISE = %.5f, ISCO = %.5f\n', ISEc, ISCOc);


%% Graphical comparison of results

% Define plot colors
cPID = [0.85, 0.37, 0.00]; % orange
cMPC = [0.00, 0.60, 0.00]; % green 
cPRG  = [0.00, 0.45, 0.70]; % blue

% Start figure
fig = figure;
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

% Subplot 1: Rate Tracking Error (Å/s) 
subplot(2,1,1);
hold on; grid on; box on;

p1 = plot(time, errorPID, '-',  'LineWidth', 2.5, 'Color', cPID); hold on
p2 = plot(time, errorMPC, ':',  'LineWidth', 2.5, 'Color', cMPC);
p3 = plot(time, errorPRG, '-.', 'LineWidth', 2.5, 'Color', cPRG);
plot(time, zeros(1,2000), '--k', 'LineWidth', 2);

ylim([-0.1,0.6])
ylabel({'Rate error ($\mathrm{\AA}$/s)'}, 'Interpreter', 'latex', 'FontSize', 16);
xlabel({'Time (s)'}, 'Interpreter', 'latex', 'FontSize', 16);
grid on
legend([p1,p2,p3], {'PID', 'MPC', 'RG-MPC'}, ...
    'Location', 'northeast','Orientation','vertical', 'Interpreter', 'latex', 'FontSize', 16);
set(gca, 'FontSize', 16, 'Box', 'on');
set(gcf, 'Color', 'w');

%%% Subplot 2: Power Input (%) 
subplot(2,1,2);
hold on; grid on; box on;

p4 = plot(time, uPID, '-',  'LineWidth', 2.5, 'Color', cPID); hold on
p5 = plot(time, uMPC, ':',  'LineWidth', 2.5, 'Color', cMPC);
p6 = plot(time, uPRG, '-.', 'LineWidth', 2.5, 'Color', cPRG);
plot(time, 60*ones(1,2000), '--k', 'LineWidth', 2);    % P upper limit
plot(time, zeros(1,2000), '--k', 'LineWidth', 2);      % P lower limit

ylim([-10,70])
ylabel({'Power (\%)'}, 'Interpreter', 'latex', 'FontSize', 16);
xlabel({'Time (s)'}, 'Interpreter', 'latex', 'FontSize', 16);
grid on
set(gca, 'FontSize', 16, 'Box', 'on');
set(gcf, 'Color', 'w');

%% Save figure 
disp("Figure saved")
savefig(gcf, 'ResultsComparison_PID_MPC_PRG.fig');     
