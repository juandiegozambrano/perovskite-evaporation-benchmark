%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: ComputeGainObs.m
%
% Description: Calcules Kalman Filter (KF) observer gain L for the 
% closed-loop augmented model H = (AH, BH, CH, DH). The observer is designed 
% according to:
%       x_hat(k+1) = AH * x_hat(k) + BH * w(k) 
%                    - L * ( y(k) - CH * x_hat(k) - DH * w(k) )
%
%   The pair (AH, CH) is assumed to be observable.
%
% Outputs:
%   - L : Kalman gain (8x2)
%   - (A_k, B_k, C_k, D_k) : State-space matrices of the observer model
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

clear all; clc;

%%% Load closed-loop model
load('ModelCL_H.mat');
AH = modelH_ss.A;  % State matrix 
BH = modelH_ss.B;  % Input matrix
CH = modelH_ss.C;  % Output-state matrix
DH = modelH_ss.D;  % Output-input matrix

% Dimensions
no = 8;     % Number of states
mo = 2;     % Number of control inputs
po = 1;     % Number of measured outputs

% Observability check
Obs_matrix = obsv(AH, CH);
if rank(Obs_matrix) < no
    error('System is not fully observable');
end
fprintf('System is fully observable\n');

% poles = linspace(0, 1, size(AH,1));  % Poles 10 times faster
% L = place(AH', CH', poles)';
% A_obs = AH - L * CH;
% B_obs = [BH, L];  
% C_obs = eye(no);     
% D_obs = zeros(no, mo + po);
% clc; close all;

% Noise statistics (Gaussian, zero mean)
sigma_noise = 0.00001378;       % Standard deviation from data
variance_noise = sigma_noise^2; % Variance

% Covariance matrices
Q = variance_noise * eye(8);    % Process noise covariance matrix Φ_x = σ²·I
R = diag([variance_noise, 0]);  % Measurement noise covariance Φ_y
N = zeros(8,2);                 % Cross-covariance correlation (zero)

% Solve discrete-time Riccati equation and extract KF gain
[P, ~, L] = dare(AH', CH', Q, R, N);
L = L'  % Kalman gain

% Observer state-space matrices with DH ≠ 0 
A_k = AH - L * CH;  
B_k = [BH - L * DH, L]; 
C_k = eye(8);          
D_k = zeros(8, 3);    

% Save 
disp('Saved Kalman Filter Gain and SS Matrices');
save('KalmanFilter.mat','L','A_k','B_k','C_k','D_k');
