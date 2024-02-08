clc 
clear 

load 'polynomial.mat';

%% B
% Initialize variables 
T0a = 25;
Delta_t  = 1;
Tspan = 600;
Q = 5;
I = -2.5;
SOC_0 = 0.5;

mdl = "HW3_virtual_testbed"
set_param(mdl, 'SimulationCommand', 'update');
simIn = Simulink.SimulationInput(mdl);
simIn = setModelParameter(simIn,"StopTime",'Tspan');
out = sim(simIn);

fprintf("The value of SOC at T=10mins is %d", out.SOC_out(end))

%put image of soc estimator in

%% C

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
I = -2.5;
SOC_0 = 0.5;
alpha = flipud(alpha)

mdl = "HW3_virtual_testbed"
set_param(mdl, 'SimulationCommand', 'update');
simIn = Simulink.SimulationInput(mdl);
simIn = setModelParameter(simIn,"StopTime",'Tspan');
out = sim(simIn);




