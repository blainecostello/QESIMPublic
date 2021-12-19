% circuit model code
N = 100;    % number of datapoints in sweep

% Set geometric parameters (these will be used for discrete QE simulation)
ds = 250E-9
A  = (25E-9)^2

% set material parameters 
%% NOTE::d1 IS THE DIELECTRIC... d2 IS NC / FERROELECTRIC.
conds = [7.9E-8, -7.9E-8];
perms = [  3.6 ,   -3];

% Set signal parameters
w = 1000*2*pi;
vapp = 1;

% Constants
e0 = 8.8541878128E-12;

% Initialize arrays
E1cpx = zeros(1,N);
E2cpx = zeros(1,N);
Ecpx  = zeros(2,N);
effRe = zeros(1,N);

gam1 = (1i*w*e0*perms(1)+conds(1))
gam2 = (1i*w*e0*perms(2)+conds(2))

%% Sweep this (loop, store values in array of size N) (Sweeping ferroelectric thickness)
for i = 1:N
    d1 = ds * (i/N);
%   d1 = ds * 0.5;
    d2 = ds - d1;
    % Calculate Complex E-field
    E1cpx(i) = 1./(d1+d2*(gam1/gam2));
    E2cpx(i) = 1./(d2+d1*(gam2/gam1));
    Eeez = [E1cpx(i), E2cpx(i)];
    % Calculate total energy in volume...
    deez = [d1,d2];
    EinEa = ((1/2)*e0.*(perms).*abs(Eeez).^2) .* deez*A; % energy in each voxel
    EtotRept = sum(EinEa);
    
    % Calculate permittivity with parallel plate assumption for capacitance
    effRe(i) = (2*EtotRept*ds)/((vapp^2)*e0*A);    

end

%% Calculat expeected permittivity peak (analytical)
% calculate optimal ferroelectric thickness for polarization catastrophe (PC)
%d1_opt = - ds / (((gam2))/((gam1) - 1))




%% Plot analytical model data, 
% Data array for permittivity values
figure
hold on
plot([1:N], effRe,'Linewidth', 2)
% vertical dashed line for expected PC
% plot(abs([real(d1_opt)*N, real(d1_opt)*N]), [min(effRe), max(effRe)], '-.k', 'Linewidth', 2)
title("Effective Permittivity of NC Laminate (Z-direction)")
xlabel("% Dielectric")
ylabel("Effective Permittivity")
legend("Simulation Model", "Analytical Model", "Location", "northeast")

%% Run Simulator


%% Plot Simulated data, note peak 




%% Quantify error

