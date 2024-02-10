%% 2
%% A
clc 
clear 


% a0 = 1 - deltat/R1C1
% b0 = Rs*Deltat/R1C1 - Rs + deltat/C1
% B1 = Rs
%% B

load 'polynomial.mat';
% Initialize variables 
T0a = 25;
Delta_t  = 1;
Tspan = 600;
Q = 5;
Current_Input = -2.5;
SOC_0 = 0.5;
mdl = "HW3_virtual_testbed"
set_param(mdl, 'SimulationCommand', 'update');
simIn = Simulink.SimulationInput(mdl);
simIn = setModelParameter(simIn,"StopTime",'Tspan');
out = sim(simIn);
Current = -2.5 * ones(601, 1);

figure;

subplot(2,2,1);
plot(out.tout, Current);
title('Input Current');
xlabel('Time');
ylabel('Current');

subplot(2,2,2);
plot(out.tout, out.V_out);
title('Terminal Voltage');
xlabel('Time');
ylabel('Voltage');

subplot(2,2,3);
plot(out.tout, out.SOC_out);
title('Calculated SOC');
xlabel('Time');
ylabel('SOC');

subplot(2,2,4);
plot(out.tout, out.Temp);
title('Temp');
xlabel('Time');
ylabel('y');

fprintf("The value of SOC at T=10mins is %d", out.SOC_out(end))

image_data = imread('HW3_2b_socestimator.PNG'); 
figure;
imshow(image_data);



%% C
clc
clear


load 'polynomial.mat';
Step_Time = [0, 100, 400, 1000, 1300, 1500];
Step_Current = [0, 5, 0, -5, 0];
T = 0:1500;
%Preallocatie input current 
Current = zeros(size(T));

for i = 1:5
    Current(T >= Step_Time(i) & T < Step_Time(i+1)) = Step_Current(i);
end

Current_Input = [0:1500;Current]';
plot(T, Current);

% Initialize variables
T0a = 25;
Delta_t  = 1;
Tspan = 1500;
Q = 5;
SOC_0 = 0.5;
alpha = flipud(alpha);

mdl = "HW3_virtual_testbed_2";
set_param(mdl, 'SimulationCommand', 'update');
simIn = Simulink.SimulationInput(mdl);
simIn = setModelParameter(simIn,"StopTime",'Tspan');
out = sim(simIn);

% Plot of all graphs 
figure;

subplot(2,2,1);
plot(out.tout, out.I_out);
title('Input Current');
xlabel('Time');
ylabel('Current');

subplot(2,2,2);
plot(out.tout, out.V_out);
title('Terminal Voltage');
xlabel('Time');
ylabel('Voltage');

subplot(2,2,3);
plot(out.tout(1:1501), out.SOC_out);
title('Calculated SOC');
xlabel('Time');
ylabel('SOC');

subplot(2,2,4);
plot(out.tout, out.Y_k);
title('Calculated Y_k');
xlabel('Time');
ylabel('y');

%% D 

% 0.25 SOC
T0a = 25;
Delta_t  = 1;
Tspan = 1500;
Q = 5;
SOC_0 = 0.25;


mdl = "HW3_virtual_testbed_2";
set_param(mdl, 'SimulationCommand', 'update');
simIn = Simulink.SimulationInput(mdl);
simIn = setModelParameter(simIn,"StopTime",'Tspan');
out_25 = sim(simIn);

%Preallocate big PHI_25
PHI_25 = zeros(1501, 3);
PHI_25(:, 1) = out_25.I_out;
PHI_25(2:end, 2) = out_25.I_out(1:end-1);
PHI_25(2:end, 3) = out_25.Y_k(1:end-1);

%solve least squares 
Y_25 = out_25.Y_k;
Theta_25 = PHI_25\Y_25
b1_25 = Theta_25(1)
b0_25 = Theta_25(2)
a0_25 = Theta_25(3)

% 0.5 SOC
T0a = 25;
Delta_t  = 1;
Tspan = 1500;
Q = 5;
SOC_0 = 0.5;


mdl = "HW3_virtual_testbed_2";
set_param(mdl, 'SimulationCommand', 'update');
simIn = Simulink.SimulationInput(mdl);
simIn = setModelParameter(simIn,"StopTime",'Tspan');
out_50 = sim(simIn);

%Preallocate big PHI_50
PHI_50 = zeros(1501, 3);
PHI_50(:, 1) = out_50.I_out;
PHI_50(2:end, 2) = out_50.I_out(1:end-1);
PHI_50(2:end, 3) = out_50.Y_k(1:end-1);

%solve least squares 
Y_50 = out_50.Y_k;
Theta_50 = PHI_50\Y_50
b1_50 = Theta_50(1)
b0_50 = Theta_50(2)
a0_50 = Theta_50(3)

% 0.75 SOC
T0a = 25;
Delta_t  = 1;
Tspan = 1500;
Q = 5;
SOC_0 = 0.75;


mdl = "HW3_virtual_testbed_2";
set_param(mdl, 'SimulationCommand', 'update');
simIn = Simulink.SimulationInput(mdl);
simIn = setModelParameter(simIn,"StopTime",'Tspan');
out_75 = sim(simIn);

%Preallocate big PHI_75
PHI_75 = zeros(1501, 3);
PHI_75(:, 1) = out_75.I_out;
PHI_75(2:end, 2) = out_75.I_out(1:end-1);
PHI_75(2:end, 3) = out_75.Y_k(1:end-1);

%solve least squares 
Y_75 = out_75.Y_k;
Theta_75 = PHI_75\Y_75
b1_75 = Theta_75(1)
b0_75 = Theta_75(2)
a0_75 = Theta_75(3)



%% E
%Espressions:
% Rs = b1
% C1 = deltat/(b0+b1a0)
% R1 = -deltat/c1(a0-1)

% 25% SOC 
Rs_25 = b1_25
C1_25 = Delta_t/(b0_25+b1_25*a0_25)
R1_25 = -Delta_t/(C1_25*(a0_25-1))
% 50% SOC 
Rs_50 = b1_50
C1_50 = Delta_t / (b0_50 + b1_50 * a0_50)
R1_50 = -Delta_t / (C1_50 * (a0_50 - 1))
% 75% SOC 
Rs_75 = b1_75
C1_75 = Delta_t / (b0_75 + b1_75 * a0_75)
R1_75 = -Delta_t / (C1_75 * (a0_75 - 1))


%% P3
%% A
OCV_Poly = [3.2152,-5.2313,3.0532,3.1264];
Delta_t  = 1;
Tspan = 1500;
Rs = 0.04;
R1 = 0.1;
C1 = 300;
Q = 5;
Soc_Init = 0.5;

mdl = "HW3_virtual_testbed_3";
set_param(mdl, 'SimulationCommand', 'update');
simIn = Simulink.SimulationInput(mdl);
simIn = setModelParameter(simIn,"StopTime",'Tspan');
out = sim(simIn);

image_data = imread('HW3_3a.PNG');

figure;
imshow(image_data);
title('3A Simulink'); 


%% B

%Preallocate big phi
PHI = zeros(1501, 3);
PHI(:, 1) = out.I_out;
PHI(2:end, 2) = out.I_out(1:end-1);
PHI(2:end, 3) = out.Y_k(1:end-1);

%solve least squares 
Y = out.Y_k;
Theta = PHI\Y
b1 = Theta(1)
b0 = Theta(2)
a0 = Theta(3)


Rs = b1;
C1 = Delta_t / (b0 + b1 * a0);
R1 = -Delta_t / (C1 * (a0 - 1));

figure;

subplot(2,2,1);
plot(out.tout, out.I_out);
title('Input Current');
xlabel('Time');
ylabel('Current');

subplot(2,2,2);
plot(out.tout, out.Vt);
title('Terminal Voltage');
xlabel('Time');
ylabel('Voltage');

subplot(2,2,3);
plot(out.tout, out.Z_soc);
title('Calculated SOC');
xlabel('Time');
ylabel('SOC');

subplot(2,2,4);
plot(out.tout, out.Y_k);
title('Calculated Y_k');
xlabel('Time');
ylabel('y');


%% C
% make table
col = {'Actual', 'Estimated'};
titles = {'0.04', '0.01', '300'}';
Nums = double([Rs,R1,C1]);
vals =  Nums';
configTable = cell2table([titles, num2cell(vals)], 'VariableNames', col);
configTable.(1) = categorical(configTable.(1));
disp(configTable)

%Parameters calculated and estimated are the same. 
