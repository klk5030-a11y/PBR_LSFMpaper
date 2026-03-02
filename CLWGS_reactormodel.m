% Chemical Engineering Journal 461 (2023) 141896
% Intensification of the reverse water_gas shift process using a
% countercurrent chemical looping regenerative reactor
% Reproduced by Yikyeom Kim (finally modified at 23/04/07)

% Matched w/ Metcalfe version
% Keep in mind that superficial velocity definition is different w/
% Metcalfe version.
% Aldo ver -> convection term - STP basis
% Metcalfe ver -> convection term - high temp. basis
% Whether flowrate determined at STP is accelerated at high T?

function packed_bed

%%% Process parameters %%%
R = 8.314;                              % Universial gas const : J/mol/K
T = 850 + 273.15;                       % Temperature : K
P = 1;                                  % Pressure : bar -> Pa
L = 0.0051;                              % Length : m
r = 0.0035;                              % reactor radius in m
A = pi()*r^2;                           % reactor cross section in m2 
Coc = 6/(355.09e-30*6.02e23);           % Molar conc : mol/m3
Cgas = P*1e5/R/T;                       % Gas phase conc : mol/m3 (c = n/V = p/RT)
eps = 0.6;                              % Bed void fraction
porosity = 0;                        % Particle internal inaccessible volume 

flowrate = 50;                          % sccm total gas flow 
flowrate2 = flowrate/1000/1000/60;      % m3/sec total gas flow
u_s = flowrate2/A * T/273.15;           % Superficial gas velocity at working T: m/s
Deff_red = 1e-5;                           % Diffusion coefficient : m2/s
Deff_ox = 1e-5;
k0 = 1;                                % Source term rate const : 1/sec
Duration_red = 58;                      % Duration : sec
Duration_ox = 60;                      % Duration : sec

cycle_number = 5;                      % Number of redox cycling : #
yC = 0.1; yH = 0.1;                     % sum of carbon and hydrogen species : vol. composition
xC = yC*Cgas; xH = yH*Cgas;             % sum conc. of carbon and hydrogen species : mol/m3
fprintf('mass of OC used in this simulation is %f g.\n', L*A*(1-eps)*Coc*228)      % mass in g
c1 = u_s/eps; c2 = Deff_red/eps; c3 = Coc*(1-porosity)*(1-eps)/eps;    

%%% Calculation accuracy parameters %%%
nt = 600; dt_red = Duration_red/nt; dt_ox = Duration_ox/nt;
nz = 200; dz = L / nz; dz2 = dz^2;
nsolid = 1; ngas = 5;
nvar = nsolid+ngas; 


%%% Thermodynamic properties of gas %%%
% 2H2 + O2 = 2H2O
g_rxn_H = -484236.17 + 86.4262.*T + 0.01253.*T.^2;
Kwat = exp(-g_rxn_H/R/T); 
% 2CO + O2 = 2CO2
g_rxn_C = -567258.15 + 177.219.*T -0.00138.*T.^2;
Kcar = exp(-g_rxn_C/R/T);

Kw = Kwat^(-0.5);                   % modification for the paper equation convention
Kc = Kcar^(-0.5);
%Kw = 1.107e-9;                      % From NIST H2O = H2+0.5O2
%Kc = 1.094e-9;                      % From NIST CO+0.5O2 = CO2
% 
delta_LSF311 = load("LSFM3131.mat").delta_LSFM3131;
pO2 = logspace(-23, 0, 50);
lpO2 = log10(pO2)/10;
coeff = polyfit(lpO2, delta_LSF311, 5);

fdelta_LSF311 = polyf(pO2, coeff); 

semilogx(pO2, delta_LSF311, 'o'); hold on
semilogx(pO2, fdelta_LSF311);


x_red = zeros(nt, nz*nvar,cycle_number);
x_ox = zeros(nt, nz*nvar,cycle_number);



tspan_red = linspace(0, Duration_red, nt);
tspan_ox = linspace(0, Duration_ox, nt);
zspan = linspace(0, L, nz);

% ODE solver Setups
% 1. Sparse Jacobian pattern matrix should be provided for speed up
% 2. Strong source term (far deviation from equilibrium) can lead to numerical instabilities 
% 2-1. Check species balance when ODE solver fails => If so, tighten calculation tolerances and parameters
% 3. Increase mesh size when numerical instabilities are observed or fails at t = 0 and z = 0

Jred = JPat_red(); Jox = JPat_ox();
    ODEoptions_red = odeset('RelTol',1e-8, 'AbsTol', 1e-12, 'JPattern', Jred, 'NonNegative', 1:(nz*nvar));
    ODEoptions_ox = odeset('RelTol',1e-8, 'AbsTol', 1e-12, 'JPattern', Jox, 'NonNegative', 1:(nz*nvar));

P_eq = 9.38e-8*ones(nz, 1);
IC_delta = polyf(P_eq, coeff); 
IC_CO = carbon_yCO(P_eq).*Cgas; %IC for oxidation half cycle
dummy = ones(nz, 1)*Cgas*1e-14;

%%% Solving cycles of chemical looping ODE system %%%
for cycle = 1:cycle_number  
    IC_red = [IC_delta; IC_CO; yC*Cgas-IC_CO; dummy; dummy; dummy];
    BC_red = [0, yC*Cgas*1, 0, 0, 0, 0];    
    tic
    [t, x_red(:,:,cycle)] = ode15s(@ode_red, tspan_red, IC_red, ODEoptions_red, BC_red); % reduction step
    toc
    fprintf("*** %dth CO Reduction cycle done ***\n", cycle);    
    
    IC_delta = x_red(nt, 1:nz, cycle)';    
    IC_H2O = water_yH2O(carbon_pO2(x_red(nt, 1+nz:2*nz,cycle)'/Cgas, x_red(nt, 1+2*nz:3*nz,cycle)'/Cgas))*Cgas;      
    IC_ox = [IC_delta; dummy; dummy; yH*Cgas-IC_H2O; IC_H2O; dummy];
    BC_ox = [0, 0, 0, 0, yH*Cgas*1, 0];
    tic
    [t, x_ox(:,:,cycle)] = ode15s(@ode_ox, tspan_ox, IC_ox, ODEoptions_ox, BC_ox); % H2O splitting step
    toc   
    fprintf("*** %dth Steam oxidation cycle done ***\n", cycle);
    
    IC_delta = x_ox(nt, 1:nz, cycle)';
    IC_CO = carbon_yCO(water_pO2(x_ox(nt, 1+4*nz:5*nz,cycle)'/Cgas, x_ox(nt, 1+3*nz:4*nz,cycle)'/Cgas))*Cgas;     
end
fname=sprintf('reactor_cycle_variables.mat');     % save the process variables
save(fname, 'x_red', 'x_ox');

%%% Post processing of ODE solution %%%
load('reactor_cycle_variables.mat');        % load the saved process variables
for cycle = 1:cycle_number
    plot_gas_profile(x_red(:,:,cycle), x_ox(:,:,cycle));
    plot_nonstoichiometry(x_red(:,:,cycle), x_ox(:,:,cycle));
    plot_partialpressure(x_red(:,:,cycle), x_ox(:,:,cycle));    
    [sol_red, sol_ox, gas_red, gas_ox] = recovery_ratio(x_red(:,:,cycle), x_ox(:,:,cycle));
    fprintf("%d cycle - solid : %.3f, %.3f (mmol), gas : %.3f, %.3f (mmol) \n", cycle, sol_ox, sol_red, gas_ox, gas_red);
end


function [delta, xCO, xCO2, xH2, xH2O, xO2] = comp_split(x)
    delta = x(:,1:nz);
    xCO = x(:, 1+nz:2*nz);
    xCO2 = x(:, 1+2*nz:3*nz);
    xH2 = x(:, 1+3*nz:4*nz);
    xH2O = x(:, 1+4*nz:5*nz);
    xO2 = x(:, 1+5*nz:6*nz);
end

function [delta, yCO, yCO2, yH2, yH2O, yO2] = comp_split2(x)
    delta = x(:,1:nz);
    yCO = x(:, 1+nz:2*nz)/Cgas;
    yCO2 = x(:, 1+2*nz:3*nz)/Cgas;
    yH2 = x(:, 1+3*nz:4*nz)/Cgas;
    yH2O = x(:, 1+4*nz:5*nz)/Cgas;
    yO2 = x(:, 1+5*nz:6*nz)/Cgas;
end
function [delta, yCO, yCO2, yH2, yH2O, yO2] = comp_split3(x)
    delta = x(1:nz);
    yCO = x(1+nz:2*nz)/Cgas;
    yCO2 = x(1+2*nz:3*nz)/Cgas;
    yH2 = x(1+3*nz:4*nz)/Cgas;
    yH2O = x(1+4*nz:5*nz)/Cgas;
    yO2 = x(1+5*nz:6*nz)/Cgas;
end


function plot_gas_profile(x_red, x_ox)  % Realtime gas-profile, nonstoichiometry
    [delta_red, yCO_red, yCO2_red, yH2_red, yH2O_red, yO2_red] = comp_split2(x_red(:,:));
    [delta_ox, yCO_ox, yCO2_ox, yH2_ox, yH2O_ox, yO2_ox] = comp_split2(x_ox(:,:));

    t = tiledlayout(2, 2, 'TileSpacing', 'Compact');
    set(gcf, 'position', [50, 50, 1200, 720]);
    title(t, sprintf('Cycle %d', cycle), 'FontSize', 14, 'FontWeight', 'Bold');
    xlabel(t, 'time span (hr)');
    ylabel(t, 'axial position (m)');

    nexttile
    pO2=CeO2_pO2(delta_red);
    plot_data = delta_red';
    imagesc(tspan_red, zspan, plot_data);
    title_OC = sprintf('Reduction step');
    title(title_OC, 'FontSize', 10);        
    grid on
    colormap jet
    colorbar
    caxis([0.0 0.5]);        

    nexttile
    pO2 = CeO2_pO2(delta_ox);
    plot_data = delta_ox';
    imagesc(tspan_ox, zspan, plot_data);
    title_OC = sprintf('Oxidation step');
    title(title_OC, 'FontSize', 10);
    grid on
    colormap jet
    colorbar
    caxis([0.0 0.5]);   

    nexttile
    plot(tspan_red, yCO_red(:,end-1), tspan_red, yCO2_red(:,end-1));
    legend('CO', 'CO2');    
    hold on
    red_title = sprintf('H2 Reduction step');
    title(red_title, 'FontSize', 10);
    xlim([0 Duration_red]);       
    ylim([0 yC]);

    nexttile     
    plot(tspan_ox, yH2_ox(:,2), tspan_ox, yH2O_ox(:,2));
    legend('H2', 'H2O');     
    hold on
    ox_title = sprintf('Oxidation step');
    title(ox_title, 'FontSize', 10);
    xlim([0 Duration_ox]);
    ylim([0 yH]);
    ylabel('composition');
    
  


    fpath = './data/';
    fname1 = sprintf('Gas_profile_cycle%d', cycle);
    saveas(gcf, fullfile(fpath, fname1), 'jpeg');
end
function plot_nonstoichiometry(x_red, x_ox) % Nonstoichiometry variation
    [delta_red, yCO_red, yCO2_red, yH2_red, yH2O_red, yO2_red] = comp_split2(x_red(:,:));
    [delta_ox, yCO_ox, yCO2_ox, yH2_ox, yH2O_ox, yO2_ox] = comp_split2(x_ox(:,:));
    
    t = tiledlayout(1, 2, 'TileSpacing', 'Compact');
    set(gcf, 'position', [50, 50, 1200, 400]);
    title(t, sprintf('Oxygen nonstoichiometry (Cycle %d)', cycle), 'FontSize', 14, 'FontWeight', 'Bold');

    nexttile    
    time_pick = [1/Duration_red; 1/6; 2/6; 3/6; 4/6; 5/6; 6/6]*Duration_red;
    for i = 1:7
        plot(zspan, delta_red(round(time_pick(i)/dt_red), :));
        hold on
    end
    xlim([0 L]);
    ylim([0 0.5]);
    legend(num2str(time_pick)); 
    xlabel('position x [m]');
    ylabel('Non-stoichiometry delta');    
    
    nexttile    
    time_pick = [1/Duration_ox; 1/6; 2/6; 3/6; 4/6; 5/6; 6/6]*Duration_ox;
    for i = 1:7
        plot(zspan, delta_ox(round(time_pick(i)/dt_ox), :));
        hold on
    end
    xlim([0 L]);
    ylim([0 0.5]);
    legend(num2str(time_pick)); 
    xlabel('position x [m]');
    ylabel('Non-stoichiometry delta');    
    
    fpath = './data/';
    fname1 = sprintf('Nonstoichiometry_cycle%d', cycle);
    saveas(gcf, fullfile(fpath, fname1), 'jpeg');
end
function plot_partialpressure(x_red, x_ox) % Oxygen partial pressure variation
    [delta_red, yCO_red, yCO2_red, yH2_red, yH2O_red, yO2_red] = comp_split2(x_red(:,:));
    [delta_ox, yCO_ox, yCO2_ox, yH2_ox, yH2O_ox, yO2_ox] = comp_split2(x_ox(:,:));
    
    pO2_ox = water_pO2(yH2O_ox, yH2_ox);
    pO2_red = carbon_pO2(yCO_red, yCO2_red);
    t = tiledlayout(1, 2, 'TileSpacing', 'Compact');
    set(gcf, 'position', [50, 50, 1200, 400]);
    title(t, sprintf('Oxygen partial pressure (cycle %d)', cycle), 'FontSize', 14, 'FontWeight', 'Bold');

    nexttile    
    time_pick = [1/Duration_red; 1/6; 2/6; 3/6; 4/6; 5/6; 6/6]*Duration_red;
    for i = 1:7
        semilogy(zspan, pO2_red(round(time_pick(i)/dt_red), :));
        hold on
    end
    xlim([0 L]);
    ylim([1e-25 1]);
    legend(num2str(time_pick)); 
    xlabel('position x [m]');
    ylabel('pO2');    
    
    nexttile    
    time_pick = [1/Duration_ox; 1/6; 2/6; 3/6; 4/6; 5/6; 6/6]*Duration_ox;
    for i = 1:7
        semilogy(zspan, pO2_ox(round(time_pick(i)/dt_ox), :));
        hold on
    end
    xlim([0 L]);
    ylim([1e-25 1]);    
    legend(num2str(time_pick)); 
    xlabel('position x [m]');
    ylabel('pO2');    
    
    fpath = './data/';
    fname1 = sprintf('oxygen partial pressure_cycle%d', cycle);
    saveas(gcf, fullfile(fpath, fname1), 'jpeg');
end

%%% Systems of ordinary differential equation %%%
function dydt = ode_red(t, x, BC)
    dydt = zeros(nz*nvar, 1);  
    delta_eq = zeros(nz, 1);    
    rxn_coeff = [1, -1, 1, 0, 0, 0];
    [delta, yCO, yCO2, yH2, yH2O, yO2] = comp_split3(x);
       
    pO2 = carbon_pO2(yCO, yCO2);
    yCsum = yCO+yCO2;
    
    delta_eq = polyf(pO2, coeff);
    delta_diff = delta_eq - delta;
    
    for k = 1:nvar                                                                              % k : index of species (CO, CO2, H2, H2O, O2)
        if k >= 4
            continue;
        end
        idx = (k-1)*nz;
        if k == 1                                                                               % oxygen carrier source term.
            reactant_limiter = sigf(yCO(1:nz)-1e-10, 1e10);
            dydt(1:nz) = rxn_coeff(k)*k0.*sigf(delta_diff(1:nz), 1e3).*delta_diff(1:nz).*reactant_limiter;
        else
            dxdz = [(x(idx+1)-BC(k)); diff(x(idx+1:idx+nz))]/dz;                                % Convection term
            d2xdz2 = [0; diff(x(idx+1:idx+nz), 2); 0]/dz2;                                      % Diffusion term - No diffusional flux at boundary         
            dydt(idx+1:idx+nz) = -c1.*dxdz(1:nz)+c2.*d2xdz2(1:nz)+rxn_coeff(k)*c3.*dydt(1:nz);                          % Conservation eqn - zone1.
        end            
    end
end
function dydt = ode_ox(t, x, BC)
    dydt = zeros(nz*nvar, 1);
    delta_eq = zeros(nz, 1);    
    rxn_coeff = [-1, 0, 0, -1, 1, 0];      
    [delta, yCO, yCO2, yH2, yH2O, yO2] = comp_split3(x);    
    yHsum = yH2+yH2O;
    pO2 = water_pO2(yH2O, yH2);    
    delta_eq = polyf(pO2, coeff);
    delta_diff = delta - delta_eq;
    
    for k = 1:nvar
        if k == 2 || k == 3 || k == 6
            continue;
        end
        idx = (k-1)*nz;
        if k == 1                                                                                % k : index of species (delta, CO, CO2, H2, H2O, O2)
            reactant_limiter = sigf(yH2O(1:nz)-1e-10, 1e10);
            dydt(idx+1:idx+nz) = rxn_coeff(k)*k0.*sigf(delta_diff(1:nz), 1e3).*delta_diff(1:nz).*reactant_limiter;
        else
            dxdz = [diff(x(idx+1:idx+nz)); (BC(k)-x(idx+nz))]/dz;                                % Convection term
            d2xdz2 = [0; diff(x(idx+1:idx+nz), 2); 0]/dz2;                                       % Diffusion term - No diffusional flux at boundary         
            dydt(idx+1:idx+nz) = c1.*dxdz(1:nz)+c2.*d2xdz2(1:nz)+rxn_coeff(k)*c3.*dydt(1:nz);                          % Conservation eqn - zone1.
        end
    end    
end

%%% Helper functions %%%
function [sum_solid_red, sum_solid_ox, sum_gas_red, sum_gas_ox] = recovery_ratio(x_red, x_ox)
    [delta_red, yCOred, yCO2red, yH2red, yH2Ored] = comp_split2(x_red);
    [delta_ox, yCOox, yCO2ox, yH2ox, yH2Oox] = comp_split2(x_ox);
    
    solid_ox_st = trapz(delta_ox(1, :));
    solid_ox_ed = trapz(delta_ox(nt, :)); 
    sum_solid_ox = trapz(delta_ox(1,:) - delta_ox(nt,:))*dz*Coc*(1-eps)*(1-porosity)*A*1000;
    sum_gas_ox = trapz(yH2ox(:, 1))*flowrate*dt_ox/22.4/60;              % mmol
    recovery_ox = sum_gas_ox/sum_solid_ox;    
    
    solid_red_st = trapz(delta_red(1, :));           
    solid_red_ed = trapz(delta_red(nt, :)); % mol/m2
    sum_solid_red = trapz(delta_red(nt,:) - delta_red(1,:))*dz*Coc*(1-eps)*(1-porosity)*A*1000;
    sum_gas_red = trapz(yCO2red(:, nz))*flowrate*dt_red/22.4/60;           % mmol
    recovery_red = sum_gas_red/sum_solid_red;
    
end

function pO2 = water_pO2(yH2O, yH2)
    y_floor = 1e-16;
    pO2 = (P).*(Kw.*max(yH2O, y_floor)./max(yH2, y_floor)).^2;
end
function yH2O = water_yH2O(pO2)
    lhs = (pO2./P).^0.5./Kw;
    yH2O = yH.*lhs./(lhs+1);
end
function pO2 = carbon_pO2(yCO, yCO2)
    y_floor = 1e-16;
    pO2 = (P).*(Kc.*max(yCO2, y_floor)./max(yCO, y_floor)).^2;
%     for i = 1:nz
%         if yCO(i) < 1e-5
%             pO2(i) = 1e-7;
%         end
%     end
end
function yCO = carbon_yCO(pO2)
    lhs = (pO2./P).^0.5./Kc;
    yCO2 = yC.*lhs./(lhs+1);
    yCO = yC-yCO2;
end

function delta_eq = CeO2_delta(pO2)
    % Thermodynamic properties of CeO2 OC
    delta_s = 165;                          % CeO2 Entropy : J/mol/K
    delta_h = 430 * 1000;                   % J/mol
    
    rhs = ((pO2./P).^(-0.5).*exp(delta_s/R).*exp(-delta_h/R/T)).^(-1/2.31);
    delta_eq = 0.35./(rhs + 1);
end
function pO2 = CeO2_pO2(delta_eq)
    % Thermodynamic properties of CeO2 OC
    delta_s = 165;                          % CeO2 Entropy : J/mol/K
    delta_h = 430 * 1000;                   % J/mol
    
    lhs = (delta_eq./(0.35-delta_eq)).^2.31;
    gibbs = exp(delta_s/R)*exp(-delta_h/R/T);
    
    pO2 = P.*(lhs./gibbs).^(-2);
end

function delta_eq = polyf(pO2, coeff)
    pO2 = min(max(pO2, 1e-23), 1);
    lpO2 = log10(pO2)/10;
    delta_eq = coeff(1).*lpO2.^5 +coeff(2).*lpO2.^4 +coeff(3).*lpO2.^3 +coeff(4).*lpO2.^2 +coeff(5).*lpO2 +coeff(6); 
end

function y = sigf(x, eta)
    % sigmoid delta function, eta can determine the slope of ramping
    y = 1./(1+exp(-eta.*(x)));
%     if x > 0
%         y = x;
%     else
%         y = 0;
%     end
end
function out = JPat_red()
    S0 = spdiags(ones(nz, 1), 0:0, nz, nz); 
    S1 = spdiags(ones(nz,3),-1:1,nz,nz);
    ZZ = zeros(nz,nz);
    out =[S0 S0 S0 ZZ ZZ ZZ;
          S0 S1 S0 ZZ ZZ ZZ;
          S0 S0 S1 ZZ ZZ ZZ;
          ZZ ZZ ZZ ZZ ZZ ZZ;
          ZZ ZZ ZZ ZZ ZZ ZZ;
          ZZ ZZ ZZ ZZ ZZ ZZ];
end
function out = JPat_ox()
    S0 = spdiags(ones(nz, 1), 0:0, nz, nz); 
    S1 = spdiags(ones(nz,3),-1:1,nz,nz);
    ZZ = zeros(nz,nz);
    out =[S0 ZZ ZZ S0 S0 ZZ;
          ZZ ZZ ZZ ZZ ZZ ZZ;
          ZZ ZZ ZZ ZZ ZZ ZZ;
          S0 ZZ ZZ S1 S0 ZZ;
          S0 ZZ ZZ S0 S1 ZZ;
          ZZ ZZ ZZ ZZ ZZ ZZ];
end
function out = JPat_ox2()
    S0 = spdiags(ones(nz, 1), 0:0, nz, nz); 
    S1 = spdiags(ones(nz,3),-1:1,nz,nz);
    ZZ = zeros(nz,nz);
    out =[S0 ZZ ZZ ZZ ZZ S0;
          ZZ ZZ ZZ ZZ ZZ ZZ;
          ZZ ZZ ZZ ZZ ZZ ZZ;
          ZZ ZZ ZZ ZZ ZZ ZZ;
          ZZ ZZ ZZ ZZ ZZ ZZ;
          S0 ZZ ZZ ZZ ZZ S1];
end

function coeff = PBFO_delta_coeff(T)
R = 8.314459;
n = 500;

s0_i = 15.6; h0_i = 133.8*1000; Gi = h0_i - T*s0_i; Ki = exp(-Gi/R/T);
s0_af = 58.6; h0_af = 102.7*1000; Gaf = h0_af - T*s0_af; Kaf = exp(-Gaf/R/T);
s0_ox = -97.7; h0_ox = -46.4*1000; Gox = h0_ox - T*s0_ox; Kox = exp(-Gox/R/T);
Kred = Ki^4*Kaf^2/Kox;

% Ki = e*h
% Kaf = Oi" + Vo
% Kox = Oi^2*h^4/pO2
% Kred = Vo^2*e^4*pO2

% Kox*Kred = Ki^4 * Kaf^2
% 2*Kaf/Vo + Ki/h = Pr_Ba(1, fixed) + h + 2*Vo
% => h^2 + (2Vo + 1 - 2Kaf/Vo)h - Ki = 0
% Solving routine : Vo -> Oi" -> h -> i -> pO2
Vo = logspace(-1.8 , -0.05, n);
Oi = Kaf./Vo;
% solving quadratic eqn. ax^2+bx+c = 0  
B = 2.*Vo+1-2.*Kaf./Vo;
C = -Ki;
h = (-B + sqrt(B.^2 - 4*C))/2;

i = Ki./h;
pO2 = Kaf^2/Kox./Vo.^2.*h.^4;
delta_eq = Oi - Vo;
delta_eq = (1-delta_eq)/2;

lpO2 = log10(pO2)/10;
coeff = polyfit(lpO2, delta_eq, 5);


end

end
