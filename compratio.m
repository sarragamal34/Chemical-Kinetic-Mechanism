%% Syngas Compression Ratio Sensitivity Study
% Variable Cp(T), Residuals, and Composition Sweep
clc; clear; close all;

%% UNIVERSAL CONSTANTS
R_u = 8.314462618; % J/mol-K

%% USER INPUTS
T1 = 298;             % Intake temperature [K]
p1 = 101325;          % Intake pressure [Pa]
residual_frac = 0.05; % Residual gas fraction (by mass)
CR_range = 6:0.5:18;  % Compression ratios
q_LHV = 1e8;          % J/kg fuel (syngas average, updated later per composition)
phi = 1.0;            % Equivalence ratio (stoichiometric)
T_res = 800;          % Residual gas temperature [K]

%% SPECIES DATA
species = {'H2','CO','CH4','CO2','N2'};
M = [2.016 28.01 16.04 44.01 28.014]*1e-3; % kg/mol

% NASA Cp coefficients (valid 300–1500K)
NASA.H2  = [3.33727920E+00, -4.94024731E-05, 4.99456778E-07, -1.79566394E-10, 2.00255376E-14];
NASA.CO  = [3.57953347E+00, -6.10353680E-04, 1.01681433E-06, 9.07005884E-10, -9.04424499E-13];
NASA.CH4 = [1.70227300E+00, 9.08170100E-03, -2.16447300E-06, 2.18663400E-10, -7.95260700E-15];
NASA.CO2 = [2.35677352E+00, 8.98459677E-03, -7.12356269E-06, 2.45919022E-09, -1.43699548E-13];
NASA.N2  = [3.29812400E+00, 1.40824040E-03, -3.96322200E-06, 5.64151500E-09, -2.44485400E-12];

%% DEFINE FUNCTIONS
Cp_T = @(T,coeff) R_u*(coeff(1)+coeff(2)*T+coeff(3)*T.^2+coeff(4)*T.^3+coeff(5)./T.^2); % J/molK

%% H2 fraction sweep (0.1 → 0.9)
H2_fracs = 0.1:0.1:0.9;
nH2 = numel(H2_fracs);
nCR = numel(CR_range);

eta = zeros(nH2, nCR);

%% MAIN LOOP
for i = 1:nH2
    x_H2 = H2_fracs(i);
    x_CO = 1 - x_H2 - 0.1;   % keep CH4 fixed at 0.1 fraction of fuel mixture
    x_CH4 = 0.1;
    x_CO2 = 0.05;
    x_N2  = 0.05;
    
    % normalize mole fractions
    x = [x_H2 x_CO x_CH4 x_CO2 x_N2];
    x = x / sum(x);
    
    % Molar mass of mixture
    M_mix = sum(x .* M);
    
    % Define Cp(T) for mixture
    Cp_mix = @(T) sum(x .* arrayfun(@(k) Cp_T(T, NASA.(species{k})), 1:length(species))) / M_mix;
    R_spec = R_u / M_mix;
    gamma_mix = @(T) Cp_mix(T) / (Cp_mix(T) - R_spec);
    
    % Include residual gas effect
    Cp_charge = @(T) (1-residual_frac)*Cp_mix(T) + residual_frac*Cp_mix(T_res);
    gamma_eff = @(T) Cp_charge(T) / (Cp_charge(T) - R_spec);
    
    % Estimate effective LHV (weighted)
    LHV = [120e6 10.1e6 50e6 0 0]; % J/kg for H2, CO, CH4, CO2, N2
    q_in = sum(x .* LHV) * 0.95; % approximate, include burn efficiency
    
    %% Compression ratio sweep
    for j = 1:nCR
        cr = CR_range(j);
        
        % variable gamma compression (numerical)
        steps = 200;
        V_ratio = linspace(1,1/cr,steps);
        Ttemp = zeros(1,steps);
        Ttemp(1) = T1;
        for k = 2:steps
            gamma = gamma_eff(Ttemp(k-1));
            Ttemp(k) = Ttemp(k-1)*(V_ratio(k-1)/V_ratio(k))^(gamma-1);
        end
        T2 = Ttemp(end);
        
        % Add heat (constant volume)
        T3 = T2 + q_in/Cp_charge(T2);
        
        % Mean gamma
        gamma_avg = mean([gamma_eff(T1), gamma_eff(T3)]);
        
        % Thermal efficiency (Otto)
        eta(i,j) = 1 - 1/(cr^(gamma_avg - 1));
    end
end

%% PLOT RESULTS
[CRg,H2g] = meshgrid(CR_range, H2_fracs);

figure('Position',[200 150 1000 600]);
surf(CRg, H2g, eta*100);
xlabel('Compression Ratio');
ylabel('H_2 Mole Fraction');
zlabel('Thermal Efficiency (%)');
title('Syngas Efficiency Map (Variable Cp, with Residuals)');
grid on; shading interp; colorbar;
view(45,30);

%% PRINT SAMPLE OUTPUT
fprintf('Syngas Sensitivity Study Completed.\n');
fprintf('Sample Results:\n');
for i=[1 round(nH2/2) nH2]
    fprintf('H2 = %.1f → η @ CR=10: %.2f%% | η @ CR=16: %.2f%%\n', ...
        H2_fracs(i), eta(i,8)*100, eta(i,end)*100);
end
%% Syngas Compression Ratio Sensitivity and Target Study
% Variable Cp(T), Residuals, and Composition Sweep
clc; clear; close all;

%% CONSTANTS
R_u = 8.314462618; % J/mol-K

%% INPUTS
T1 = 298;             % K, intake temperature
p1 = 101325;          % Pa, intake pressure
residual_frac = 0.05; % by mass
CR_range = 6:0.5:18;  % compression ratios
T_res = 800;          % K, residual gas temp
target_eff = 0.35;    % desired efficiency (35%)

%% SPECIES DATA
species = {'H2','CO','CH4','CO2','N2'};
M = [2.016 28.01 16.04 44.01 28.014]*1e-3; % kg/mol

% NASA Cp coefficients (300–1500 K)
NASA.H2  = [3.33727920E+00,-4.94024731E-05,4.99456778E-07,-1.79566394E-10,2.00255376E-14];
NASA.CO  = [3.57953347E+00,-6.10353680E-04,1.01681433E-06,9.07005884E-10,-9.04424499E-13];
NASA.CH4 = [1.70227300E+00,9.08170100E-03,-2.16447300E-06,2.18663400E-10,-7.95260700E-15];
NASA.CO2 = [2.35677352E+00,8.98459677E-03,-7.12356269E-06,2.45919022E-09,-1.43699548E-13];
NASA.N2  = [3.29812400E+00,1.40824040E-03,-3.96322200E-06,5.64151500E-09,-2.44485400E-12];

%% FUNCTIONS
Cp_T = @(T,coeff) R_u*(coeff(1)+coeff(2)*T+coeff(3)*T.^2+coeff(4)*T.^3+coeff(5)./T.^2); % J/mol K

%% COMPOSITION SWEEP
H2_fracs = 0.1:0.1:0.9;  % vary H2
nH2 = numel(H2_fracs);
CRs = CR_range;
nCR = numel(CRs);
eta = zeros(nH2,nCR);
req_CR = zeros(1,nH2); % to reach target η

%% MAIN LOOP
for i = 1:nH2
    x_H2  = H2_fracs(i);
    x_CO  = 1 - x_H2 - 0.1;
    x_CH4 = 0.1;
    x_CO2 = 0.05;
    x_N2  = 0.05;
    x = [x_H2 x_CO x_CH4 x_CO2 x_N2];
    x = x / sum(x);

    M_mix = sum(x.*M);
    R_spec = R_u/M_mix;
    
    Cp_mix = @(T) sum(x .* arrayfun(@(k) Cp_T(T,NASA.(species{k})),1:length(species))) / M_mix;
    Cp_charge = @(T) (1-residual_frac)*Cp_mix(T) + residual_frac*Cp_mix(T_res);
    gamma_eff = @(T) Cp_charge(T) / (Cp_charge(T) - R_spec);
    
    LHV = [120e6 10.1e6 50e6 0 0]; % J/kg
    q_in = sum(x .* LHV) * 0.95;
    
    for j = 1:nCR
        cr = CRs(j);
        steps = 200;
        V_ratio = linspace(1,1/cr,steps);
        Ttemp = zeros(1,steps); Ttemp(1)=T1;
        for k = 2:steps
            gamma = gamma_eff(Ttemp(k-1));
            Ttemp(k)=Ttemp(k-1)*(V_ratio(k-1)/V_ratio(k))^(gamma-1);
        end
        T2 = Ttemp(end);
        T3 = T2 + q_in/Cp_charge(T2);
        gamma_avg = mean([gamma_eff(T1),gamma_eff(T3)]);
        eta(i,j) = 1 - 1/(cr^(gamma_avg-1));
    end

    % Find required CR for target efficiency
    if any(eta(i,:) >= target_eff)
        req_CR(i) = interp1(eta(i,:), CRs, target_eff, 'linear');
    else
        req_CR(i) = NaN;
    end
end

%% PLOTS
[CRg,H2g] = meshgrid(CRs,H2_fracs);
figure('Position',[200 150 1000 600]);
surf(CRg,H2g,eta*100);
xlabel('Compression Ratio'); ylabel('H_2 Mole Fraction');
zlabel('Thermal Efficiency (%)');
title('Syngas Efficiency Map (Variable Cp, Residuals)');
grid on; shading interp; colorbar;

% Required CR plot
figure('Position',[250 200 700 400]);
plot(H2_fracs, req_CR,'-o','LineWidth',1.8);
xlabel('H_2 Mole Fraction'); ylabel('Required CR for Target Efficiency');
title(sprintf('Compression Ratio Needed for %.0f%% Efficiency', target_eff*100));
grid on;

%% SAMPLE OUTPUT
fprintf('Syngas Compression Ratio Sensitivity Study Complete.\n');
fprintf('Target efficiency: %.1f%%\n',target_eff*100);
for i=[1 round(nH2/2) nH2]
    fprintf('H2 = %.1f → Required CR = %.2f\n', H2_fracs(i), req_CR(i));
end
