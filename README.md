# Perovskite Evaporation Benchmark
This repository contains all MATLAB and Simulink scripts developed for the data-driven modeling and predictive control of perovskite deposition. In addition to the implementation files, the repository provides the Excel datasets with six experimental samples of PbI₂ deposition at a rate of 0.5 Å/s:

Scripts must be runned according to the following workflow:

<img width="1804" height="782" alt="image" src="https://github.com/user-attachments/assets/38dbe589-9f63-498a-84cc-9528c27dfb5e" />

To run the benchmark, it is necessary to install MATLAB® R2021b (or later versions) with standard toolboxes (System Identification Toolbox, Control system Toolbox, Signal Processing Toolbox, and Optimization Toolbox are required) and ensure that Simulink is available. 


A step-by-step summary of the files to be executed for modeling is shown hereunder:
1. (Data preprocessing) The model (1) for thermal evaporation dynamics of a material is identified based on experimental data using the MATLAB System Identification Toolbox. The script readingDataPbI2.m  is provided in the repository to first load and prepare the PbI2 datasets. Next, you have to run ModelTraining.m to proceed with the identification process, where the datasets are divided into training and validation subsets.
    
 2. (Model identification) Run script ModelIdentification.m to estimate subsystems 𝐺1 and 𝐺2 based on previously processed data, and obtain the unified grey-box model 𝐺 in (1) of the PbI2 vaporation process. At the end, the model matrices are saved in ModelOL_G.mat.
    
 3. (Closed-loop model) The PRG method requires the closed-loop model. Running ModelClosedLoop.m generates the PID state-space model and the final closed-loop model, which is saved in ModelCL_H.mat. The repository also includes OverallSysSim.slx, a block-based Simulink file for the PID control using the open-loop model 𝐺, employed to validate the identified 8th-order closed-loop model that requires the PRG design. Once the process models are obtained, the files that need to be executed for implementing the proposed control methods are the following:

 1. (MPC method): The MPC problem is solved in MATLAB using the quadprog solver. You must runMPCalgorithm.m,inwhichthesimulationparameters(simulation length, prediction horizon, sampling time, etc) and the weights of the cost function 𝐽 are defined, to solve the optimal control problem. It also generates a plot with the power input and rate trajectories obtained over the simulation.
    
 2. (PRG method): The PRG approach is implemented in MATLAB using the CasADI toolbox with the ipopt solver. This approach considers the closed loop dynamics stored in ModelCL_H.mat and the full closed-loop state of the process, whichcanbeestimatedviathestateobserver in ComputeGainObs.m. The PRG problem is implemented and solved in PRGalgorithm.m, which also generates graphical results for the rate, temperature, and virtual rate reference throughout the simulation.
    
 Finally, a graphical comparison of the three control methods (PID, MPC, and PRG)  is produced by running ResultsComparison.m, which also calculates the performance metrics ISE and ISCO.
