clc
clear all
close all

% circuit model code
N = 200;    % number of datapoints in sweep

% Set geometric parameters (these will be used for discrete QE simulation)
ds = 250E-9
A  = (25E-9)^2

% set material parameters 
%% NOTE::d1 IS THE DIELECTRIC... d2 IS NC / FERROELECTRIC.
conds = [7.9E-8, -7.9E-8];
perms = [  3.6 ,   -3];

% Set signal parameters
w = 1E-2*2*pi;
vapp = 1;

% Constants
e0 = 8.8541878128E-12;

% Initialize arrays
E1cpx = zeros(1,N);
E2cpx = zeros(1,N);
Ecpx  = zeros(2,N);
effRe = zeros(1,N);

% Initialize Video
X1 = VideoWriter("NC_Surf"+wexp+".mp4");
X1.FrameRate = 20;
open(X1);


gam1 = (1i*w*e0*perms(1)+conds(1))

condsNC = -1*linspace(1E-8,5E-7,N);
permsNC = -1*linspace(1,10,N);

%% Sweep this (loop, store values in array of size N) (Sweeping ferroelectric thickness)
for i = 1:N % This i index will denote the frame within the animation.
    for j = 1:N  % This j index will denote the real NC permittivity.
        for k = 1:N  % This k index will denote the NC conductivity.
            d1 = ds * (j/N);
            %   d1 = ds * 0.5;
            d2 = ds - d1;
            gam2 = (1i*w*e0*permsNC(i)+condsNC(k));
            % Calculate Complex E-field
            E1cpx(i,j,k) = 1./(d1+d2*(gam1/gam2));
            E2cpx(i,j,k) = 1./(d2+d1*(gam2/gam1));
            Eeez = [E1cpx(i,j,k), E2cpx(i,j,k)];
            % Calculate total energy in volume...
            deez = [d1,d2];
            perms(2) = permsNC(i);
            EinEa = ((1/2)*e0.*(perms).*abs(Eeez).^2) .* deez*A; % energy in each voxel
            EtotRept = sum(EinEa);

            % Calculate permittivity with parallel plate assumption for capacitance
            effRe(j,k) = (2*EtotRept*ds)/((vapp^2)*e0*A);
        end
    end
    subplot(1,2,1)
    hold on
    surf(condsNC, [1:N]/2, effRe);
    ylabel("% of FE")
    xlabel("Conductivity of FE")
    colorbar
    caxis([-1E2, 1E2])
    title("Top View (Permittivity = "+permsNC(i)+")")
    
    subplot(1,2,2)
    surf(condsNC,[1:N]/2, effRe);
    zlim([-1E2,1E2]);
    % vertical dashed line for expected PC
%     plot(abs([real(d1_opt)*N, real(d1_opt)*N]), [min(effRe), max(effRe)], '-.k', 'Linewidth', 2)
    title("Effective Permittivity of NC Laminate")
    ylabel("% of FE")
    xlabel("Conductivity of FE")
    zlabel("Effective Composite Permittivity")
    view([60,30]);
    colorbar
    caxis([-2E2, 2E2])

    set(gcf, "Position", [000,100,1400,500])
    
    % Save frame
    e2p(i) = getframe(gcf);
    writeVideo(X1, e2p(i).cdata);
    
    
    close
end

%% Calculat expeected permittivity peak (analytical)
% calculate optimal ferroelectric thickness for polarization catastrophe (PC)
% d1_opt = - ds / (((gam2))/((gam1) - 1))


close(X1)

%% Plot analytical model data, 
% Data array for permittivity values


%% Run Simulator


%% Plot Simulated data, note peak 




%% Quantify error

