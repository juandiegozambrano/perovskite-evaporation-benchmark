%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: ReadingDataPbI2.m
%
% Description: This script extracts from EXCEL files the six experimental 
% datasets of PbI2 with rate reference at 0.5 Å/s. By running this script,
% you obtain:
%   - Power input,
%   - Measured temperature and its setpoint,
%   - Deposition rate,
%   - Thickness and cumulated thickness,
% for each dataset.
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

%% Reading PbI2 Dataset 1

% Detect import options, skipping first 4 header lines
opts = detectImportOptions('Src2_PbI2_Master_0_5As_v201111 2023.12.12-15.37.44.CSV', 'NumHeaderLines', 4);
opts.VariableTypes(1:2) = {'char', 'char'};  % Ensure first two columns are treated as strings

% Read the whole CSV into a table
source2a_aux = readtable('Src2_PbI2_Master_0_5As_v201111 2023.12.12-15.37.44.CSV', opts);

% Extract date and time columns separately and merge
date_text = source2a_aux{:, 1};
time_text = source2a_aux{:, 2};
datetime_text = strcat(date_text, {' '}, time_text);

% Convert to datetime format
src2a.t = datetime(datetime_text, 'InputFormat', 'MMM-dd-yyyy HH:mm:ss.SSS', 'Locale', 'en_US');
src2a.t.Format = 'MMM-dd-yyyy HH:mm:ss.SSS';

% Sampling time calculation
dtt = diff(src2a.t);  
dt = seconds(dtt);
Ts1 = mode(dt); % Sampling time

% The third column contains concatenated values separated by commas
concatenated_data = source2a_aux{:, 3};
sp_data = split(concatenated_data, ',');

% Time vector
time1 = src2a.t;

% Extract variables
d1s1au = str2double(sp_data(:,4));   % Power Supply
d1s2au = str2double(sp_data(:,15));  % Temperature
d1s3au = str2double(sp_data(:,14));  % Deposition Rate
d1s4au = str2double(sp_data(:,16));  % Temperature Setpoint
d1s5au = str2double(sp_data(:,28));  % Thickness
d1s6au = str2double(sp_data(:,12));  % Cummulated thickness

% Plot of Dataset 1
f1 = figure;
subplot(4,1,1)
plot(time1,d1s1au)
ylabel('Power Supply [%]')

subplot(4,1,2)
plot(time1,d1s2au)
hold on
plot(time1,d1s4au)
ylabel('Temperature [°C]')
legend('Measured','Setpoint')

subplot(4,1,3)
plot(time1,d1s3au)
ylabel('Deposition Rate [Å/s]')

subplot(4,1,4)
plot(time1,d1s5au)
hold on
plot(time1,d1s6au)
ylabel('Thickness [Å]')
legend('Measured','Cumulative')

sgtitle('PbI2 - Dataset 1')
f1.WindowState = 'maximized';

%% Reading PbI2 Dataset 2

% Detect import options and skip header rows
opts = detectImportOptions('Src2_PbI2_Master_0_5As_v201111 2024.01.11-09.17.24.CSV', 'NumHeaderLines', 4);
opts.VariableTypes(1:2) = {'char', 'char'};  % Date and time as strings

% Load CSV data
source2a_aux = readtable('Src2_PbI2_Master_0_5As_v201111 2024.01.11-09.17.24.CSV', opts);

% Extract date and time fields and merge
date_text = source2a_aux{:, 1};
time_text = source2a_aux{:, 2};
datetime_text = strcat(date_text, {' '}, time_text);

% Convert to datetime format
src2a.t = datetime(datetime_text, 'InputFormat', 'MMM-dd-yyyy HH:mm:ss.SSS', 'Locale', 'en_US');
src2a.t.Format = 'MMM-dd-yyyy HH:mm:ss.SSS';

% Sampling time calculation
dtt = diff(src2a.t);  
dt = seconds(dtt);
Ts2 = mode(dt); % Sampling time

% The third column contains concatenated values separated by commas
concatenated_data = source2a_aux{:, 3};
sp_data = split(concatenated_data, ',');

% Time vector
time2 = src2a.t;

% Extract variables
d2s1au = str2double(sp_data(:,4));   % Power Supply
d2s2au = str2double(sp_data(:,15));  % Temperature
d2s3au = str2double(sp_data(:,14));  % Deposition Rate
d2s4au = str2double(sp_data(:,16));  % Temperature Setpoint
d2s5au = str2double(sp_data(:,28));  % Thickness
d2s6au = str2double(sp_data(:,12));  % Cumulative Thickness

% Plot of Dataset 2
f2 = figure;
subplot(4,1,1)
plot(time2,d2s1au)
ylabel('Power Supply [%]')

subplot(4,1,2)
plot(time2,d2s2au)
hold on
plot(time2,d2s4au)
ylabel('Temperature [°C]')
legend('Measured','Setpoint')

subplot(4,1,3)
plot(time2,d2s3au)
ylabel('Deposition Rate [Å/s]')

subplot(4,1,4)
plot(time2,d2s5au)
hold on
plot(time2,d2s6au)
ylabel('Thickness [Å]')
legend('Measured','Cumulative')

sgtitle('PbI2 - Dataset 2')
f2.WindowState = 'maximized';

km = 756;               % Sample index where the substrate shutter opens
tpid = 0:1:1999;        % Time vector of 2000 samples (assuming 1 unit per sample)
rpid = d2s3au(km:2755); % Extract 2000-point segment of deposition rate data starting from shutter opening
upid = d2s1au(km:2755); % Extract 2000-point segment of power data starting from shutter opening

% Save results
disp('Saved dataset: PID results ');
save('PIDresults.mat','tpid','rpid','upid');

% ------------ Code to generate a dataset comparative image --------------
% subplot(2,1,1)
% plot(time2, d2s2au, 'b', 'LineWidth', 2); 
% grid on;
% xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 12);
% ylabel('Temperature [°C]', 'Interpreter', 'latex', 'FontSize', 12);
% title('\textbf{Data response - PbI$_2$}', 'Interpreter', 'latex', 'FontSize', 14);
% set(gca, 'FontSize', 12, 'FontName', 'Times');
% 
% subplot(2,1,2)
% plot(time2, d2s3au, 'r', 'LineWidth', 2);
% grid on;
% xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 12);
% ylabel('Rate [Å/s]', 'Interpreter', 'latex', 'FontSize', 12);
% set(gca, 'FontSize', 12, 'FontName', 'Times');
% sgtitle('\textbf{Data response - PbI$_2$}', 'Interpreter', 'latex', 'FontSize', 16);
%-------------------------------------------------------------------------------------

%% Reading PbI2 Dataset 3

% Detect import options, skipping first 4 header lines
opts = detectImportOptions('Src2_PbI2_Master_0_5As_v201111 2024.01.31-09.22.37.CSV', 'NumHeaderLines', 4);
opts.VariableTypes(1:2) = {'char', 'char'};  % Ensure first two columns are treated as strings

% Read the whole CSV into a table
source2a_aux = readtable('Src2_PbI2_Master_0_5As_v201111 2024.01.31-09.22.37.CSV', opts);

% Extract date and time columns separately and merge
date_text = source2a_aux{:, 1};
time_text = source2a_aux{:, 2};
datetime_text = strcat(date_text, {' '}, time_text);

% Convert to datetime format
src2a.t = datetime(datetime_text, 'InputFormat', 'MMM-dd-yyyy HH:mm:ss.SSS', 'Locale', 'en_US');
src2a.t.Format = 'MMM-dd-yyyy HH:mm:ss.SSS';

% Sampling time calculation
dtt = diff(src2a.t);
dt = seconds(dtt);
Ts3 = mode(dt); % Sampling time

% The third column contains concatenated values separated by commas
concatenated_data = source2a_aux{:, 3};
sp_data = split(concatenated_data, ',');

% Time vector
time3 = src2a.t;

% Extract variables
d3s1au = str2double(sp_data(:,4));   % Power Supply
d3s2au = str2double(sp_data(:,15));  % Temperature
d3s3au = str2double(sp_data(:,14));  % Deposition Rate
d3s4au = str2double(sp_data(:,16));  % Temperature Setpoint
d3s5au = str2double(sp_data(:,28));  % Thickness
d3s6au = str2double(sp_data(:,12));  % Cummulated thickness

% Plot of Dataset 3
f3 = figure;
subplot(4,1,1)
plot(time3,d3s1au)
ylabel('Power Supply [%]')

subplot(4,1,2)
plot(time3,d3s2au)
hold on
plot(time3,d3s4au)
ylabel('Temperature [°C]')
legend('Measured','Setpoint')

subplot(4,1,3)
plot(time3,d3s3au)
ylabel('Deposition Rate [Å/s]')

subplot(4,1,4)
plot(time3,d3s5au)
hold on
plot(time3,d3s6au)
ylabel('Thickness [Å]')
legend('Measured','Cumulative')

sgtitle('PbI2 - Dataset 3')
f3.WindowState = 'maximized';

%% Reading PbI2 Dataset 4

% Detect import options, skipping first 4 header lines
opts = detectImportOptions('Src2_PbI2_Master_0_5As_v201111 2024.02.15-11.55.35.CSV', 'NumHeaderLines', 4);
opts.VariableTypes(1:2) = {'char', 'char'};  % Ensure first two columns are treated as strings

% Read the whole CSV into a table
source2a_aux = readtable('Src2_PbI2_Master_0_5As_v201111 2024.02.15-11.55.35.CSV', opts);

% Extract date and time columns
date_text = source2a_aux{:, 1};
time_text = source2a_aux{:, 2};

% Split third column values (comma-separated)
concatenated_data = source2a_aux{:, 3};
sp_data = split(concatenated_data, ',');

% Extract AM/PM info (first entry of concatenated data)
ap = sp_data(:, 1);

% Combine date, time, and AM/PM into datetime strings
datetime_text = strcat(date_text, {' '}, time_text, {' '}, ap);

% Convert to datetime with AM/PM formatting
src2a.t = datetime(datetime_text, 'InputFormat', 'MMM-dd-yyyy hh:mm:ss.SSS a', 'Locale', 'en_US');
src2a.t.Format = 'MMM-dd-yyyy hh:mm:ss.SSS a';

% Sampling time calculation
dtt = diff(src2a.t);
dt = seconds(dtt);
Ts4 = mode(dt);

% Time vector
time4 = src2a.t;

% Extract variables
d4s1au = str2double(sp_data(:,4));   % Power Supply
d4s2au = str2double(sp_data(:,15));  % Temperature
d4s3au = str2double(sp_data(:,14));  % Deposition Rate
d4s4au = str2double(sp_data(:,16));  % Temperature Setpoint
d4s5au = str2double(sp_data(:,28));  % Thickness
d4s6au = str2double(sp_data(:,12));  % Cummulated thickness

% Plot of Dataset 4
f4 = figure;
subplot(4,1,1)
plot(time4,d4s1au)
ylabel('Power Supply [%]')

subplot(4,1,2)
plot(time4,d4s2au)
hold on
plot(time4,d4s4au)
ylabel('Temperature [°C]')
legend('Measured','Setpoint')

subplot(4,1,3)
plot(time4,d4s3au)
ylabel('Deposition Rate [Å/s]')

subplot(4,1,4)
plot(time4,d4s5au)
hold on
plot(time4,d4s6au)
ylabel('Thickness [Å]')
legend('Measured','Cumulative')

sgtitle('PbI2 - Dataset 4')
f4.WindowState = 'maximized';

%% Reading PbI2 Dataset 5

% Detect import options, skipping first 4 header lines
opts = detectImportOptions('Src2_PbI2_Master_0_5As_v201111 2024.05.07-09.33.32.CSV', 'NumHeaderLines', 4);
opts.VariableTypes(1:2) = {'char', 'char'};  % Ensure first two columns are treated as strings

% Read the whole CSV into a table
source2a_aux = readtable('Src2_PbI2_Master_0_5As_v201111 2024.05.07-09.33.32.CSV', opts);

% Extract date and time columns separately and merge
date_text = source2a_aux{:, 1};
time_text = source2a_aux{:, 2};
datetime_text = strcat(date_text, {' '}, time_text);

% Convert to datetime format
src2a.t = datetime(datetime_text, 'InputFormat', 'MMM-dd-yyyy HH:mm:ss.SSS', 'Locale', 'en_US');
src2a.t.Format = 'MMM-dd-yyyy HH:mm:ss.SSS';

% Sampling time calculation
dtt = diff(src2a.t);
dt = seconds(dtt);
Ts5 = mode(dt); % Sampling time

% The third column contains concatenated values separated by commas
concatenated_data = source2a_aux{:, 3};
sp_data = split(concatenated_data, ',');

% Time vector
time5 = src2a.t;

% Extract variables
d5s1au = str2double(sp_data(:,4));   % Power Supply
d5s2au = str2double(sp_data(:,15));  % Temperature
d5s3au = str2double(sp_data(:,14));  % Deposition Rate
d5s4au = str2double(sp_data(:,16));  % Temperature Setpoint
d5s5au = str2double(sp_data(:,28));  % Thickness
d5s6au = str2double(sp_data(:,12));  % Cummulated thickness

% Plot of Dataset 5
f5 = figure;
subplot(4,1,1)
plot(time5,d5s1au)
ylabel('Power Supply [%]')

subplot(4,1,2)
plot(time5,d5s2au)
hold on
plot(time5,d5s4au)
ylabel('Temperature [°C]')
legend('Measured','Setpoint')

subplot(4,1,3)
plot(time5,d5s3au)
ylabel('Deposition Rate [Å/s]')

subplot(4,1,4)
plot(time5,d5s5au)
hold on
plot(time5,d5s6au)
ylabel('Thickness [Å]')
legend('Measured','Cumulative')

sgtitle('PbI2 - Dataset 5')
f5.WindowState = 'maximized';


%% Reading PbI2 Dataset 6

% Detect import options, skipping first 4 header lines
opts = detectImportOptions('Src2_PbI2_Master_0_5As_v201111 2024.05.08-09.46.37.CSV', 'NumHeaderLines', 4);
opts.VariableTypes(1:2) = {'char', 'char'};  % Ensure first two columns are treated as strings

% Read the whole CSV into a table
source2a_aux = readtable('Src2_PbI2_Master_0_5As_v201111 2024.05.08-09.46.37.CSV', opts);

% Extract date and time columns separately and merge
date_text = source2a_aux{:, 1};
time_text = source2a_aux{:, 2};
datetime_text = strcat(date_text, {' '}, time_text);

% Convert to datetime format
src2a.t = datetime(datetime_text, 'InputFormat', 'MMM-dd-yyyy HH:mm:ss.SSS', 'Locale', 'en_US');
src2a.t.Format = 'MMM-dd-yyyy HH:mm:ss.SSS';

% Sampling time calculation
dtt = diff(src2a.t);
dt = seconds(dtt);
Ts6 = mode(dt); % Sampling time

% The third column contains concatenated values separated by commas
concatenated_data = source2a_aux{:, 3};
sp_data = split(concatenated_data, ',');

% Time vector
time6 = src2a.t;

% Extract variables
d6s1au = str2double(sp_data(:,4));   % Power Supply
d6s2au = str2double(sp_data(:,15));  % Temperature
d6s3au = str2double(sp_data(:,14));  % Deposition Rate
d6s4au = str2double(sp_data(:,16));  % Temperature Setpoint
d6s5au = str2double(sp_data(:,28));  % Thickness
d6s6au = str2double(sp_data(:,12));  % Cummulated thickness

% Plot of Dataset 6
f6 = figure;
subplot(4,1,1)
plot(time6,d6s1au)
ylabel('Power Supply [%]')

subplot(4,1,2)
plot(time6,d6s2au)
hold on
plot(time6,d6s4au)
ylabel('Temperature [°C]')
legend('Measured','Setpoint')

subplot(4,1,3)
plot(time6,d6s3au)
ylabel('Deposition Rate [Å/s]')

subplot(4,1,4)
plot(time6,d6s5au)
hold on
plot(time6,d6s6au)
ylabel('Thickness [Å]')
legend('Measured','Cumulative')

sgtitle('PbI2 - Dataset 6')
f6.WindowState = 'maximized';


