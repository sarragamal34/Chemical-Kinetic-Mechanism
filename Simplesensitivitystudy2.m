%% Advanced Syngas Sensitivity Study 

clc; clear; close all;

%% --- Constants ---
R = 8.314462618;        % J/(mol·K)
T0 = 298.15;            % K
P0 = 101325;            % Pa
M_air = 28.9647e-3;     % kg/mol
O2_frac_air = 0.21;
rho_air_ref = 1.184;    % kg/m^3 at 25°C (for Wobbe reference)

%% --- Component database (row vectors) ---
species = {'H2','CO','CH4','CO2','N2'};

M = [2.016, 28.01, 16.04, 44.01, 28.014] * 1e-3;   % kg/mol (row)
LHV = [120.0, 10.1, 50.0, 0.0, 0.0];              % MJ/kg (row)
O2_permol = [0.5, 0.5, 2.0, 0.0, 0.0];            % mol O2 per mol fuel
RON_proxy = [120, 0, 100, 0, 0];                  % placeholder knock index (row)
Cp = [14.3, 29.1, 35.7, 37.1, 29.1];              % J/mol·K (row)

%% --- Example syngas compositions (mole fractions) ---
syngas_examples = [
    0.20, 0.30, 0.05, 0.10, 0.35;
    0.40, 0.40, 0.05, 0.05, 0.10;
    0.10, 0.20, 0.02, 0.08, 0.60;
    0.60, 0.30, 0.05, 0.00, 0.05;
];
names = {'A','B','C','D'};

fprintf('Syngas property summary (25°C, 1 atm):\n');
fmt = '%s: Mmix=%.4f kg/mol | rho=%.3f kg/m3 | LHV=%.2f MJ/kg | AFR=%.2f | Cp=%.3f kJ/kgK | RON_proxy=%.1f\n';

for i = 1:size(syngas_examples,1)
    x = syngas_examples(i,:);
    [Mmix, rho, LHVm, AFRm, RONm, Cp_mass] = compute_props_row(x, M, LHV, O2_permol, RON_proxy, Cp, P0, R, T0, M_air, O2_frac_air);

    % crude adiabatic flame temperature estimate (LHV -> sensible heat)
    T_ad = T0 + (LHVm * 1e6) / (Cp_mass * 1000); % K (very crude)

    % Wobbe index: LHV_vol / sqrt(relative density)
    LHV_vol = LHVm * rho; % MJ/m^3
    WI = LHV_vol / sqrt(rho / rho_air_ref);

    fprintf(fmt, names{i}, Mmix, rho, LHVm, AFRm, Cp_mass, RONm);
    fprintf('    -> T_ad (est) = %.0f K,  Wobbe Index = %.3f MJ/m^3\n', T_ad, WI);
end

%% --- Sensitivity sweep: vary H2 mole fraction 0..0.8 ---
H2_grid = linspace(0,0.8,41);
CH4_fixed = 0.02;
CO2_fixed = 0.03;
results = nan(numel(H2_grid),8); % columns: Mmix,rho,LHV,AFR,RON,Cp,Tad,WI

for k = 1:numel(H2_grid)
    h2 = H2_grid(k);
    remaining = 1 - h2 - CH4_fixed - CO2_fixed;
    if remaining <= 0
        continue; % leave NaNs
    end
    co = remaining * 0.6;
    n2 = remaining - co;
    x = [h2, co, CH4_fixed, CO2_fixed, n2];
    [Mmix, rho, LHVm, AFRm, RONm, Cp_mass] = compute_props_row(x, M, LHV, O2_permol, RON_proxy, Cp, P0, R, T0, M_air, O2_frac_air);
    T_ad = T0 + (LHVm * 1e6) / (Cp_mass * 1000);
    LHV_vol = LHVm * rho;
    WI = LHV_vol / sqrt(rho / rho_air_ref);
    results(k,:) = [Mmix, rho, LHVm, AFRm, RONm, Cp_mass, T_ad, WI];
end

%% --- Plot results (simple) ---
figure('Name','Syngas Sensitivity Analysis','NumberTitle','off','Position',[100 100 1000 600]);

subplot(2,3,1);
plot(H2_grid, results(:,1)*1e3,'LineWidth',1.2);
xlabel('H_2 mole fraction'); ylabel('M_{mix} (g/mol)'); grid on;

subplot(2,3,2);
plot(H2_grid, results(:,2),'LineWidth',1.2);
xlabel('H_2 mole fraction'); ylabel('Density (kg/m^3)'); grid on;

subplot(2,3,3);
plot(H2_grid, results(:,3),'LineWidth',1.2);
xlabel('H_2 mole fraction'); ylabel('LHV (MJ/kg)'); grid on;

subplot(2,3,4);
plot(H2_grid, results(:,4),'LineWidth',1.2);
xlabel('H_2 mole fraction'); ylabel('AFR (mass)'); grid on;

subplot(2,3,5);
plot(H2_grid, results(:,7),'LineWidth',1.2);
xlabel('H_2 mole fraction'); ylabel('T_{ad} (K)'); grid on;

subplot(2,3,6);
plot(H2_grid, results(:,8),'LineWidth',1.2);
xlabel('H_2 mole fraction'); ylabel('Wobbe Index (MJ/m^3)'); grid on;

sgtitle('Syngas Thermochemical Property Sensitivity vs H_2 fraction');

%% --- Local function definition ---
function [Mmix, rho, LHV_mass, AFR_mass, RONp, Cp_mass_kJkg] = compute_props_row(x, M, LHV, O2_permol, RON_proxy, Cp, P0, R, T0, M_air, O2_frac_air)
    % Ensure row vector
    x = reshape(x,1,[]);
    s = sum(x);
    if s <= 0
        error('Composition sum must be positive.');
    end
    x = x / s;

    % Mmix (kg/mol)
    Mmix = sum(x .* M);

    % density (ideal gas) kg/m^3
    rho = P0 * Mmix / (R * T0);

    % LHV (mass-based)
    mass_frac = (x .* M) / Mmix;
    LHV_mass = sum(mass_frac .* LHV);

    % Stoichiometric AFR (mass basis)
    O2_needed = sum(x .* O2_permol);
    mol_air = O2_needed / O2_frac_air;
    AFR_mass = (mol_air * M_air) / Mmix;

    % RON proxy
    RONp = sum(x .* RON_proxy);

    % Cp (kJ/kg·K)
    Cp_mix_molar = sum(x .* Cp);
    Cp_mass_kJkg = (Cp_mix_molar / Mmix) / 1000;
end

