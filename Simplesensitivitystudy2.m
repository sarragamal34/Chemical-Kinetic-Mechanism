%% Advanced Syngas Sensitivity Study (NASA Cp + Residuals)
clc; clear; close all;

%% --- Constants ---
R = 8.314462618;        % J/mol·K
T0 = 298.15;            % K
P0 = 101325;            % Pa
M_air = 28.9647e-3;     % kg/mol
O2_frac_air = 0.21;
rho_air_ref = 1.184;    % kg/m3 at 25°C (for Wobbe reference)
residual_frac = 0.05;   % 5% residual gas fraction (approx)

%% --- NASA polynomial coefficients (200–1000 K range) ---
% a1 a2 a3 a4 a5 a6 a7  (Low temperature coefficients)
nasa.H2  = [3.298124, 8.249442e-04, -8.143015e-07, -9.475434e-11, 4.134872e-13, -1012.521, 3.508409];
nasa.CO  = [3.579533, -6.103537e-04, 1.016814e-06, 9.070059e-10, -9.044245e-13, -14344.086, 3.508409];
nasa.CH4 = [5.149876, -1.367097e-02, 4.918005e-05, -4.847431e-08, 1.666939e-11, -10246.647, -4.641303];
nasa.CO2 = [2.356774, 8.984597e-03, -7.123562e-06, 2.45919e-09, -1.436995e-13, -48372.0, 9.901052];
nasa.N2  = [3.531005, -1.236609e-04, -5.029994e-07, 2.435306e-09, -1.408812e-12, -1046.976, 2.967474];

%% --- Species database ---
species = {'H2','CO','CH4','CO2','N2'};
M = [2.016, 28.01, 16.04, 44.01, 28.014]*1e-3;   % kg/mol
LHV = [120.0, 10.1, 50.0, 0.0, 0.0];             % MJ/kg
O2_permol = [0.5, 0.5, 2.0, 0.0, 0.0];           % mol O2/mol fuel
RON_proxy = [120, 0, 100, 0, 0];                 % placeholder knock index

%% --- Example syngas compositions ---
syngas_examples = [
    0.20, 0.30, 0.05, 0.10, 0.35;
    0.40, 0.40, 0.05, 0.05, 0.10;
    0.10, 0.20, 0.02, 0.08, 0.60;
    0.60, 0.30, 0.05, 0.00, 0.05;
];
names = {'A','B','C','D'};

fprintf('Syngas property summary (25°C, 1 atm):\n');
fmt = '%s: Mmix=%.4f kg/mol | rho=%.3f kg/m3 | LHV=%.2f MJ/kg | AFR=%.2f | RON_proxy=%.1f\n';

for i = 1:size(syngas_examples,1)
    x = syngas_examples(i,:);
    [Mmix, rho, LHVm, AFRm, RONm, Tad, WI] = ...
        compute_props_NASA_residuals(x, M, LHV, O2_permol, RON_proxy, nasa, ...
        P0, R, T0, M_air, O2_frac_air, residual_frac, rho_air_ref);
    fprintf(fmt, names{i}, Mmix, rho, LHVm, AFRm, RONm);
    fprintf('   -> T_ad = %.0f K, Wobbe Index = %.3f MJ/m^3\n', Tad, WI);
end

%% --- Sensitivity sweep: vary H2 mole fraction 0..0.8 ---
H2_grid = linspace(0,0.8,41);
CH4_fixed = 0.02; CO2_fixed = 0.03;
results = nan(numel(H2_grid),7);

for k = 1:numel(H2_grid)
    h2 = H2_grid(k);
    remain = 1 - h2 - CH4_fixed - CO2_fixed;
    if remain <= 0, continue; end
    co = remain * 0.6;
    n2 = remain - co;
    x = [h2, co, CH4_fixed, CO2_fixed, n2];
    [Mmix, rho, LHVm, AFRm, RONm, Tad, WI] = ...
        compute_props_NASA_residuals(x, M, LHV, O2_permol, RON_proxy, nasa, ...
        P0, R, T0, M_air, O2_frac_air, residual_frac, rho_air_ref);
    results(k,:) = [Mmix, rho, LHVm, AFRm, RONm, Tad, WI];
end

%% --- Plot results ---
figure('Name','Syngas Sensitivity (NASA Cp + Residuals)','Position',[100 100 1000 600]);
subplot(2,3,1); plot(H2_grid, results(:,1)*1e3,'LineWidth',1.2);
xlabel('H_2 mole fraction'); ylabel('M_{mix} (g/mol)'); grid on;
subplot(2,3,2); plot(H2_grid, results(:,2),'LineWidth',1.2);
xlabel('H_2 mole fraction'); ylabel('\rho (kg/m^3)'); grid on;
subplot(2,3,3); plot(H2_grid, results(:,3),'LineWidth',1.2);
xlabel('H_2 mole fraction'); ylabel('LHV (MJ/kg)'); grid on;
subplot(2,3,4); plot(H2_grid, results(:,4),'LineWidth',1.2);
xlabel('H_2 mole fraction'); ylabel('AFR (mass)'); grid on;
subplot(2,3,5); plot(H2_grid, results(:,6),'LineWidth',1.2);
xlabel('H_2 mole fraction'); ylabel('T_{ad} (K)'); grid on;
subplot(2,3,6); plot(H2_grid, results(:,7),'LineWidth',1.2);
xlabel('H_2 mole fraction'); ylabel('Wobbe Index (MJ/m^3)'); grid on;
sgtitle('Syngas Property Sensitivity using NASA Cp and Residuals');

%% ====== Local Functions ======
function [Mmix, rho, LHV_mass, AFR_mass, RONp, Tad, WI] = ...
    compute_props_NASA_residuals(x, M, LHV, O2_permol, RON_proxy, nasa, ...
    P0, R, T0, M_air, O2_frac_air, residual_frac, rho_air_ref)

    % Normalize composition
    x = x/sum(x);
    if residual_frac>0
        x(5) = x(5) + residual_frac; % Add to N2
        x = x/sum(x);
    end

    % Mmix, density, and stoichiometric AFR
    Mmix = sum(x.*M);
    rho = P0*Mmix/(R*T0);
    mass_frac = (x.*M)/Mmix;
    LHV_mass = sum(mass_frac.*LHV);
    O2_needed = sum(x.*O2_permol);
    mol_air = O2_needed/O2_frac_air;
    AFR_mass = mol_air*M_air/Mmix;
    RONp = sum(x.*RON_proxy);

    % Cp function (NASA)
    Cp_mass_func = @(T) mix_Cp(T,x,nasa)/Mmix;

    % Energy balance for adiabatic flame temperature
    q_release = LHV_mass*1e6;
    Cp_int = @(T) integral(@(t)Cp_mass_func(t),T0,T,'ArrayValued',true);
    fun = @(T) Cp_int(T) - q_release;
    try
        Tad = fzero(fun, [1000 3500]);
    catch
        Tad = NaN;
    end

    % Wobbe index
    LHV_vol = LHV_mass * rho;
    WI = LHV_vol / sqrt(rho/rho_air_ref);
end

function Cp_molar = mix_Cp(T,x,nasa)
    species = {'H2','CO','CH4','CO2','N2'};
    Cp_species = zeros(numel(T),numel(species));
    for j = 1:numel(species)
        coeff = nasa.(species{j});
        Cp_species(:,j) = nasa_Cp(coeff,T);
    end
    Cp_molar = sum(Cp_species.*x,2);
end

function Cp = nasa_Cp(coeff,T)
    a1=coeff(1); a2=coeff(2); a3=coeff(3);
    a4=coeff(4); a5=coeff(5);
    Cp_R = a1 + a2.*T + a3.*T.^2 + a4.*T.^3 + a5.*T.^4;
    Cp = Cp_R * 8.314462618; % J/mol·K
end


