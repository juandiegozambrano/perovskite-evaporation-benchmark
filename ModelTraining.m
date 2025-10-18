%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: ModelTraining.m
%
% Description: This script performs the training of the model by diving the
% dataset into training and validation subsets. By running this script, you 
% obtain:
%   - The estimated model for G1 (Output: Temperature)
%   - The estimated model for G2 (Output: Deposition rate)
%
% Requirements:
%   - System Identification Toolbox (iddata, tfest, compare, tfdata, merge)
%   - Control System Toolbox (tf, ss, tfdata)
%   - Signal Processing Toolbox (seconds, timeseries)
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

%%% Dataset 1
d1u1 = d1s1au;                              % Input: Power supply signal for dataset 1
d1y1 = d1s2au - d1s2au(1);                  % Output: Temperature deviation from initial value
data = iddata(d1y1,d1u1,Ts1);               % Create iddata object for system identification (input-output data)
np = 3;                                     % Number of poles for transfer function estimation
nz = 2;                                     % Number of zeros for transfer function estimation
nk = 0;                                     % Input-output delay
d1systf1 = tfest(data,np,nz,nk,'Ts',Ts1);   % Estimate discrete-time transfer function
figure;
compare(data,d1systf1);                     % Compare measured output with estimated system
[d1num, d1den] = tfdata(d1systf1, 'v');     % Extract numerator and denominator of transfer function

t = time1 - time1(1);               % Time vector starting at zero
ta = seconds(t);                    % Convert time difference to seconds
d1s1sim = timeseries(d1u1,ta);      % Create time series for input
d1s2sim = timeseries(d1s2au,ta);    % Create time series for measured output

d1y2 = d1s3au;                      % Secondary output: Deposition rate
d1s3sim = timeseries(d1y2,ta);      % Time series for deposition rate

d1y1f = d1y1(775:3069);         % Filtered temperature data (selected range)
d1y2f = d1y2(775:3069);         % Filtered deposition rate data (selected range)
data = iddata(d1y2f,d1y1f,Ts1); % Create iddata object for filtered data

plot(data)                                  % Plot input-output data
np = 3; nz = 0; nk = 0;                     % Define model structure
d1systf2 = tfest(data,np,nz,nk,'Ts',Ts1);   % Estimate transfer function for filtered data
figure;
compare(data,d1systf2);                     % Compare filtered data with estimated system
[d1num2, d1den2] = tfdata(d1systf2, 'v');   % Extract numerator and denominator

%%% Dataset 2
d2u1 = d2s1au;                              % Input: Power supply signal
d2y1 = d2s2au;                              % Output: Temperature
data = iddata(d2y1,d2u1,Ts2);               % Create iddata object
np = 3; nz = 2; nk = 0;                     % Model structure
d2systf1 = tfest(data,np,nz,nk,'Ts',Ts2);   % Transfer function estimation
figure;
compare(data,d2systf1);                  % Compare measured and modeled output
[d2num, d2den] = tfdata(d2systf1, 'v');  % Extract numerator and denominator of transfer function

t = time2 - time2(1);
ta = seconds(t);
d2s1sim = timeseries(d2u1,ta); % Time series for input
d2s2sim = timeseries(d2y1,ta); % Time series for measured output

d2y2 = d2s3au;                  % Secondary output: Deposition rate
d2s3sim = timeseries(d2y2,ta);  % Time series for deposition rate

d2y1f = d2y1(759:6063);                     % Filtered temperature
d2y2f = d2y2(759:6063);                     % Filtered deposition rate
data = iddata(d2y2f,d2y1f,Ts2);             % Filtered input-output data
plot(data)                                  % Plot input-output data
np = 3; nz = 1; nk = 0;                     % Model structure for filtered data
d2systf2 = tfest(data,np,nz,nk,'Ts',Ts2); 
figure;
compare(data,d2systf2);
[d2num2, d2den2] = tfdata(d2systf2, 'v');

%%% Dataset 3

d3u1 = d3s1au; 
d3y1 = d3s2au; 
data = iddata(d3y1,d3u1,Ts3);
np = 3; nz = 2; nk = 0;
d3systf1 = tfest(data,np,nz,nk,'Ts',Ts3);
figure;
compare(data,d3systf1);
[d3num, d3den] = tfdata(d3systf1, 'v');

t = time3 - time3(1);
ta = seconds(t);
d3s1sim = timeseries(d3u1,ta);
d3s2sim = timeseries(d3y1,ta);

d3y2 = d3s3au; 
d3s3sim = timeseries(d3y2,ta);

d3y1f = d3y1(736:6063); 
d3y2f = d3y2(736:6063); 
data = iddata(d3y2f,d3y1f,Ts3);
plot(data)
np = 3; nz = 0; nk = 0;
d3systf2 = tfest(data,np,nz,nk,'Ts',Ts3);
figure;
compare(data,d3systf2);
[d3num2, d3den2] = tfdata(d3systf2, 'v');

%%% Dataset 4

d4u1 = d4s1au; 
d4y1 = d4s2au; 
data = iddata(d4y1,d4u1,Ts4);
np = 3; nz = 2; nk = 0;
d4systf1 = tfest(data,np,nz,nk,'Ts',Ts4);
figure;
compare(data,d4systf1);
[d4num, d4den] = tfdata(d4systf1, 'v');

t = time4 - time4(1);
ta = seconds(t);
d4s1sim = timeseries(d4u1,ta);
d4s2sim = timeseries(d4y1,ta);

d4y2 = d4s3au;
d4s3sim = timeseries(d4y2,ta);

d4y1f = d4y1(736:6063); 
d4y2f = d4y2(736:6063); 
data = iddata(d4y2f,d4y1f,Ts4);
plot(data)

np = 3; nz = 0; nk = 0;
d4systf2 = tfest(data,np,nz,nk,'Ts',Ts4);
figure;
compare(data,d4systf2);
[d4num2, d4den2] = tfdata(d4systf2, 'v');

%%% Dataset 5

d5u1 = d5s1au; 
d5y1 = d5s2au; 
data = iddata(d5y1,d5u1,Ts5);
np = 3; nz = 2; nk = 0;
d5systf1 = tfest(data,np,nz,nk,'Ts',Ts5);
figure;
compare(data,d5systf1);
[d5num, d5den] = tfdata(d5systf1, 'v');

t = time5 - time5(1);
ta = seconds(t);
d5s1sim = timeseries(d5u1,ta);
d5s2sim = timeseries(d5y1,ta);

d5y2 = d5s3au;
d5s3sim = timeseries(d5y2,ta);

d5y1f = d5y1(26:10270); 
d5y2f = d5y2(26:10270); 
data = iddata(d5y2f,d5y1f,Ts5);
plot(data)

np = 3; nz = 1; nk = 0;
d5systf2 = tfest(data,np,nz,nk,'Ts',Ts5);
figure;
compare(data,d5systf2);
[d5num2, d5den2] = tfdata(d5systf2, 'v');

%%% Dataset 6

d6u1 = d6s1au; 
d6y1 = d6s2au; 
data = iddata(d6y1,d6u1,Ts6);
np = 3; nz = 2; nk = 0;
d6systf1 = tfest(data,np,nz,nk,'Ts',Ts6);
figure;
compare(data,d6systf1);
[d6num, d6den] = tfdata(d6systf1, 'v');

t = time6 - time6(1);
ta = seconds(t);
d6s1sim = timeseries(d6u1,ta);
d6s2sim = timeseries(d6y1,ta);

d6y2 = d6s3au;
d6s3sim = timeseries(d6y2,ta);

d6y1f = d6y1(753:10708); 
d6y2f = d6y2(753:10708); 
data = iddata(d6y2f,d6y1f,Ts6);
plot(data)
np = 3; nz = 1; nk = 0;
d6systf2 = tfest(data,np,nz,nk,'Ts',Ts6);
figure;
compare(data,d6systf2);
[d6num2, d6den2] = tfdata(d6systf2, 'v');
