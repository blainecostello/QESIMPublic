% clc 
% clear all
% close all
function [error] = AnaValid3(N,wExp)
% Updated version with equations as derived in the LaTeX document...

% circuit model code
%N = 200;    % number of datapoints in sweep

% Set geometric parameters (these will be used for discrete QE simulation)
ds = 250E-9
x  = (7.5E-9);
y  = (7.5E-9);
A = x*y

GD = [3 3 N];
e0 = 8.8541878128E-12;

dVol = [x, y, ds] ./ GD


% set material parameters 
%% NOTE::d1 IS THE DIELECTRIC... d2 IS NC / FERROELECTRIC.
conds = [7.9E-8, 7.9E-8];
% conds = [0 0];
perms = [  3.6 ,   1];



% Set signal parameters
w = (10^(wExp))*2*pi;
vapp = 1;

volts = [-1 -1; -1 -1; 0 1];

cpxPm = perms - 1i .* (conds/(e0 * w));


% Constants
e0 = 8.8541878128E-12;

% Initialize arrays
E1cpx = zeros(1,N-1);
E2cpx = zeros(1,N-1);
Ecpx  = zeros(2,N-1);
effRe = zeros(1,N-1);

gam2 = (1i*w*e0*perms(1)-conds(1))
gam1 = (1i*w*e0*perms(2)-conds(2))

  ds = ((ds/(GD(3)))*(GD(3)-1));

    
%% Sweep this (loop, store values in array of size N) (Sweeping ferroelectric thickness)
for i = 1:N-2
    
    %% Analytical
    d2 = ds * (i/(N-1));
%   d1 = ds * 0.5;
    d1 = ds - d2;
    % Calculate Complex E-field
    E1cpx(i) = 1./(d1+d2*(gam1/gam2));
    E2cpx(i) = 1./(d2+d1*(gam2/gam1));
    Eeez = [E1cpx(i), E2cpx(i)];
    % Calculate total energy in volume...
    deez = [d1,d2];
%     EinEa = ((1/2)*e0.*abs(cpxPm).*abs(Eeez).^2) .* deez*A; % energy in each voxel
%     EtotRept = sum(EinEa);
    
    % Calculate permittivity with parallel plate assumption for capacitance
%     effReANA(i) = (2*EtotRept*ds)/((vapp^2)*e0*A); 
    
    effReANA(i) = ds*(((real(cpxPm(1))/d1)*(d1/cpxPm(1))^2 + (real(cpxPm(2))/d2)*(d2/cpxPm(2))^2)/((d1/cpxPm(1)+d2/cpxPm(2))^2));
    
    
    
    %% Simulation
        % Reset Geometry
    Omega = zeros(GD);
    Omega(:,:,1:round((GD(3)*((i/(N)))))) = 1;
    
    % Compute Stuff
    [Evec, E2p, Potential] = electrifyQE(Omega, volts, cpxPm, dVol);

    
    
    % {Permittivity}
    EinEa = ((1/2) * e0 .* real(cpxPm(Omega(:,:,2:GD(3))+1)) .* abs(E2p).^2) .* prod(dVol); % energy in each voxel
    Ed = sum(EinEa, 3);
    Ed = sum(Ed, 2);
    EtotRept = sum(Ed, 1);
    effReSIM(i) = (2*EtotRept*dVol(3)*(GD(3)-1))/((volts(3,2)^2)*e0*dVol(1)*GD(1)*dVol(2)*GD(2));
    
    

end

%% Calculat expeected permittivity peak (analytical)
% calculate optimal ferroelectric thickness for polarization catastrophe (PC)
%d1_opt = - ds / (((gam2))/((gam1) - 1))




%% Plot analytical model data, 
% Data array for permittivity values
figure
hold on
plot([1:N-2], abs(effReANA), '-x',  'Linewidth', 2)
plot([1:N-2], abs(effReSIM), '-o','Linewidth', 2)
% vertical dashed line for expected PC
% plot(abs([real(d1_opt)*N, real(d1_opt)*N]), [min(effRe), max(effRe)], '-.k', 'Linewidth', 2)
title("Effective Permittivity of NC Laminate (Z-direction)")
xlabel("% Ferroelectric")
ylabel("Effective Permittivity")
legend("Analytical Model", "Simulation Model", "Location", "northwest")

set(gcf, 'Position', [100 100 500 300])
print(gcf, '-djpeg', "anaValid_" + wExp + "_"+ N +".jpeg");

%% Run Simulator
probeInd = round((50E-9/ds)*N)
error = abs(real(effReANA(probeInd)) - real(effReSIM(probeInd)))./real(effReANA(probeInd))
%% Plot Simulated data, note peak 




%% Quantify error



end