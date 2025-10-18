%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: ModelClosedLoop.m
%
% Description: This script constructs the closed-loop model considering the 
% cascade PID control scheme. By running this file, you obtain:
%   - The state-space models of both PID controllers.
%   - The final closed-loop state-space model H.
%   - The initial state for simulating the closed loop.
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

%% PID Controller Parameters

% Outer PID1 parameters (rate loop)
kp1 = 630;
ki1 = 1/200;
kd1 = 0;

% Inner PID2 parameters (temperature loop)
kp2 = 1.91;
ki2 = 1/393;
kd2 = 1/47;
denp = kp2 + kp2*kd2-1;

%% State-space representation of the two PIDs

% Outer PID1 (Ac1, Bc1, Cc1, Dc1)
Ac1 = [0 kp1*ki1; 0 1];
Bc1 = [kp1*ki1; 1];
Cc1 = [1 0];
Dc1 = kp1;

% Inner PID2 (Ac2, Bc2, Cc2, Dc2)
Ac2 = [0 kp2*ki2; 0 1];
Bc2 = [kp2*ki2-kp2*kd2; 1];
Cc2 = [1 0];
Dc2 = kp2 + kp2*kd2;

%% Load process model
load("ModelOL_G.mat");

AG;             % Process AG matrix
BG;             % Process BG matrix
C2 = CG(1,:);   % Actually it is matrix C2 (rate) 
C1 = CG(2,:);   % Actually it is matrix C1 (temperature) 
DG;             % Process DG matrix
D1 = DG(2);     % This is D1 (power-temp)

%% Agreggated closed-loop system state
% Note that this state can alternatively be estimated using the observer in ComputeGainObserver.m

initial_state_pid1 = [315 0];   % PID_1 initial state: xc1(0) = [v_1(k); sum of e_1 from t=0 to k]
initial_state_pid2 = [60 0];    % PID_2 initial state: xc2(0) = [v_2(k); sum of e_2 from t=0 to k]
initial_state_closedLoop = [initial_state_pid1, initial_state_pid2, initial_state_G];  % Full aggregated state

%% Closed-loop state-space matrices (plant with cascade PID control scheme)
% For equations details, please refer to the book chapter.
beta = 1/(1+Dc2*D1);       
sigma = beta*Bc2*D1*Dc2*Dc1*C2 + beta*Bc2*D1*Dc2*C1;

AH = [Ac1 zeros(2,2) -Bc1*C2;
      Bc2*Cc1-beta*Bc2*D1*Dc2*Cc1,  Ac2-beta*Bc2*D1*Cc2,  -Bc2*Dc1*C2-Bc2*C1+sigma;
      beta*BG*Dc2*Cc1,    beta*BG*Cc2,     AG-beta*BG*Dc2*Dc1*C2-beta*BG*Dc2*C1];

BH = [Bc1;
      Bc2*Dc1-beta*Bc2*D1*Dc2*Dc1;
      beta*BG*Dc2*Dc1];

CH = [beta*DG*Dc2*Cc1, beta*DG*Cc2, CG-beta*DG*Dc2*Dc1*C2-beta*DG*Dc2*C1];

DH = beta*DG*Dc2*Dc1;

modelH_ss = ss(AH,BH,CH,DH);  % Closed-loop model in state space 

disp('Saved Closed-Loop Model H');
save('ModelCL_H.mat','modelH_ss','initial_state_closedLoop');
