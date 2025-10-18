%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: ModelIdentification.m
%
% Description: This script identified the thermal evaporation model for the
% PbI2. By running this file, you obtain:
%   - The subsystem G1 (power input - temperature)
%   - The subsystem G2 (temperature - rate):
%   - Post-threshold model of the process G (power input - rate): AG,
%   BG,CG,DG
%   - Initial states of G after volatility threshold: initial_state_G
%
% Requirements:
%   - System Identification Toolbox (for iddata, tfest, ssest, merge)
%   - Control System Toolbox (for tf, ss, tfdata, ssdata, balred)
%   - Signal Processing Toolbox (for butter, filtfilt)
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

close all; clc; 

%% Datasets for system G1 (Input: Power Supply, Output: Temperature)

% Dataset 1A
d1u1 = d1s1au; 
d1y1 = d1s2au - d1s2au(1);          % Remove initial temperature offset
data1 = iddata(d1y1, d1u1, Ts1);    % Create iddata object for system ID

% Dataset 2A
d2u1 = d2s1au; 
d2y1 = d2s2au - d2s2au(1); % Remove initial offset
data2 = iddata(d2y1, d2u1, Ts2);

% % Dataset 3A
% d3u1 = d3s1au; 
% d3y1 = d3s2au - d3s2au(1); % Remove initial offset
% data3 = iddata(d3y1, d3u1, Ts3);
% 
% % Dataset 4A
% d4u1 = d4s1au; 
% d4y1 = d4s2au - d4s2au(1); % Remove initial offset
% data4 = iddata(d4y1, d4u1, Ts4);
% 
% % Dataset 5A
% d5u1 = d5s1au; 
% d5y1 = d5s2au - d5s2au(1); % Remove initial offset
% data5 = iddata(d5y1, d5u1, Ts5);

% Dataset 6A
d6u1 = d6s1au; 
d6y1 = d6s2au - d6s2au(1); % Remove initial offset
data6 = iddata(d6y1, d6u1, Ts6);

%% Transfer Function Identification for G1

np = 2;  % Number of poles in the model
nz = 2;  % Number of zeros
nk = 0;  % Input-output delay

% Preprocess input using a low-pass Butterworth filter
fs = 1 / Ts1;                               % Sampling frequency
fc = 0.42;                                  % Cutoff frequency in Hz
[b, a] = butter(4, fc / (fs / 2), 'low');   % 4th-order low-pass filter

% Filter training datasets
data1.InputData = filtfilt(b, a, data1.InputData); 
data2.InputData = filtfilt(b, a, data2.InputData);
data6.InputData = filtfilt(b, a, data6.InputData);

% Merge selected datasets into one training dataset
data_train = merge(data1, data2, data6);

% Set transfer function estimation options
options = tfestOptions;
options.Focus = 'stability';            % Prioritize model stability
options.InitialCondition = 'estimate';  % Estimate initial conditions

% Estimate discrete-time transfer function G1
d1systf1 = tfest(data_train, np, nz, nk, 'Ts', Ts1, options); 
[d1num, d1den] = tfdata(d1systf1, 'v'); % Extract numerator and denominator

% % Convert transfer function to discrete-time state-space
% Gz1 = tf(d1num, d1den, Ts1);
% [Ad1, Bd1, Cd1, Dd1] = ssdata(Gz1);

%% Datasets for system G2 (Input: Temperature, Output: Deposition Rate)

% Dataset 1B
d1y1f = d1y1(775:3069); 
d1y2f = d1y2(775:3069); 
data7 = iddata(d1y2f, d1y1f, Ts1);

% Dataset 2B
d2y1f = d2y1(759:6063); 
d2y2f = d2y2(759:6063); 
data8 = iddata(d2y2f, d2y1f, Ts2);

% Dataset 3B
d3y1f = d3y1(736:6063); 
d3y2f = d3y2(736:6063); 
data10 = iddata(d3y2f, d3y1f, Ts3);

% Dataset 4B
d4y1f = d4y1(736:6063); 
d4y2f = d4y2(736:6063); 
data11 = iddata(d4y2f, d4y1f, Ts4);

% Dataset 6B
d6y1f = d6y1(753:10708); 
d6y2f = d6y2(753:10708); 
data9 = iddata(d6y2f, d6y1f, Ts6);

%% Transfer Function Identification for G2

np = 3;  % Number of poles
nz = 1;  % Number of zeros
nk = 0;  % Delay

% Merge selected datasets for training
data_train = merge(data7, data8, data9);

% TF estimation options
options = tfestOptions;
options.Focus = 'stability';            % Prioritize stability
options.InitialCondition = 'estimate';  % Estimate initial conditions
options.InputOffset = [-d1s2au(1) -d2s2au(1) -d6s2au(1)]; % Remove input offsets

% Estimate discrete-time transfer function G2
d1systf2 = tfest(data_train, np, nz, nk, 'Ts', Ts1, options); 
[d1num2, d1den2] = tfdata(d1systf2, 'v');

%% Convert G1 and G2 Transfer Functions to State-Space

% Subsystem G1
G1_tf = tf(d1num, d1den, Ts1); 
G1_ss = ss(G1_tf);              % Convert to state-space
G1_ss_red = balred(G1_ss,1);    % Reduce G1 to first-order system via 
% balanced truncation to improve the accuracy by converting the output 
% into an expression that depends on the input value (D1̸=0).

% Subsystem G2
G2_tf = tf(d1num2, d1den2, Ts1);    % Transfer function
G2_ss = ss(G2_tf);                  % State-space


%% Overall process model G

Ag1 = G1_ss_red.A;  % State matrix of reduced-order G1
Bg1 = G1_ss_red.B;  % Input matrix of G1
Cg1 = G1_ss_red.C;  % Output matrix of G1
Dg1 = G1_ss_red.D;  % Direct transmission term of G1

Ag2 = G2_ss.A;  % State matrix of G2
Bg2 = G2_ss.B;  % Input matrix of G2
Cg2 = G2_ss.C;  % Output matrix of G2
Dg2 = G2_ss.D;  % Direct transmission term of G2

ug1 = 50;           % Constant input (Power Supply)
xg1 = [];           % Initialize state vector
yg1 = [];           % Initialize output vector

xg1(1) = 0;         % Initial state
yg1(1) = 0;         % Initial output

% Simulate G1 response for k_m = 756 samples
% NOTE: The control model considers only post-threshold dynamics,
% starting once the volatility temperature is reached at sample k_m = 756 
for i = 1:756
    xg1(i+1) = Ag1*xg1(i) + Bg1*ug1;                % State update
    yg1(i+1) = Cg1*xg1(i) + Dg1*ug1 + d2s2au(1);    % Output with offset
end

% Display final state and output of G1
xg1(end); % -3.0577e+04
yg1(end); % 235.3242

% Combined state-space matrices for obtaining the overall process model G:
Ag_combined = [Ag1, zeros(1,3); Bg2 * Cg1, Ag2];    % State matrix AG      
Bg_combined = [Bg1; Bg2 * Dg1];                     % Input matrix BG        
Cg_combined = [0, Cg2];                             % Matrix CG 
Dg_combined = Dg2;                                  % Matrix DG                

% Adjust Original Simulation Data
data_original = d2s1sim.Data;          % Original input data
d3ori = d2s3sim.Data;                  % Original output 2 (deposition rate)
d2ori = d2s2sim.Data;                  % Original output 1 (temperature)
time_original = d2s1sim.Time;          % Original input time
time3_ori = d2s3sim.Time;              % Original time for output 2
time2_ori = d2s2sim.Time;              % Original time for output 1

% Take adjusted data starting from sample 757
% This ensures the model focuses on post-threshold dynamics
data_adjusted = data_original(757:end);
time_adjusted = time_original(757:end) - time_original(757); % Re-origin time

data3_adjusted = d3ori(757:end); 
time3_adjusted = time3_ori(757:end) - time_original(757);

data2_adjusted = d2ori(757:end); 
time2_adjusted = time2_ori(757:end) - time_original(757);

% Create new timeseries objects
d2s1sim_ad = timeseries(data_adjusted, time_adjusted);
d2s3sim_ad = timeseries(data3_adjusted, time3_adjusted);
d2s2sim_ad = timeseries(data2_adjusted, time2_adjusted);

% Compute initial states for G2 using steady-state assumption: 
%   x = A*x + B*u with input equal to the final output of G1
xg2_initial = Ag2 \ (-Bg2 * yg1(end)); 

% Combine initial states of G1 and G2
initial_state_G = [xg1(end), xg2_initial'];

% Combined State-Space Matrices for Simulation
AG = Ag_combined;                           % AG matrix
BG = Bg_combined;                           % BG matrix
CG = [Cg_combined; [Cg1*1.0611,0,0,0]];     % CG matrix with scaling 
% This scaling = 1.0611 is a fitting adjustment for static gain obtained 
% empirically 
DG = [Dg_combined; Dg1];                    % DG matrix 
CGp = [Cg1*1.0611,0,0,0];                   % CG matrix for plotting
DGp = Dg1;                                  % DG matrix for plotting

disp('Saved Open-Loop Model G');
save('ModelOL_G.mat','AG','BG','CG','DG','initial_state_G');