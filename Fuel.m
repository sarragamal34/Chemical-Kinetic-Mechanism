clc; clear;
%% Chemical Kinematic Mechanism
% K function 
Kinetic_Mechanism_data=readtable('Combustion Mechanism.xlsx');
T=1000;
%%
K_data=K_Arrhenius(Kinetic_Mechanism_data,T);