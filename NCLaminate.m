% Simulation of stacked NCFET type capacitive devices

% Constants
e0 = 8.8541878128E-12;


% Geometric Simulation Parameters
GD = [10 10 100];
tSamp = 100E-9
Omega = zeros(GD);
Omega(:,:,(GD(3)*0.94):GD(3)) = 2;
Omega(:,:,1:(GD(3)*0.06)) = 2;
Omega(:,:,(GD(3)*0.39):(GD(3)*0.62)) = 2;
Omega(:,:,(GD(3)*0.45):(GD(3)*0.55)) = 1;

dVol = tSamp ./ GD;

% Electronic Simulation Parameters
freq = 1E3;
w = freq * (2*pi);
volts = [-1 -1; -1 -1; 0 1]

% Material Simulation Parameters
perms = [3.9 ,  1 ,  -1  ];
conds = [1E-8, 1E3, -1E-8];
cpxPm = perms - 1i .* (conds/(e0 * w))

% Run E-Field Approximation
[Evec, E2p, Potential] = electrifyQE(Omega, volts, cpxPm, dVol)

% Calculate Effective Parameters
% {Permittivity, Charge, Loss Tangent, etc...}
EinEa = ((1/2)*e0 .* abs(cpxPm(Omega(:,:,1:GD(3)-1)+1)).* abs(E2p).^2) .* prod(dVol); % energy in each voxel
Ed = sum(EinEa, 3);
Ed = sum(Ed, 2);
EtotRept = sum(Ed, 1);
effRe = (2*EtotRept*dVol(3)*(GD(3)-1))/((volts(3,2)^2)*e0*dVol(1)*GD(1)*dVol(2)*GD(2))  


% Question:  Would a voltage sweep show non-linear response or is this
%            outside of the scope of this simulator?   (Test this.)

subplot(4, 1, 1); 
imagesc(squeeze(rot90(real(Omega(5,:,:))))); h = colorbar; set(get(h,'title'),'string','Phase')
title("Material Arrangement")


subplot(4, 1, 2); 
imagesc(squeeze(rot90(real(Potential(5,:,:))))); h = colorbar; set(get(h,'title'),'string','Volts');
title("Scalar Potential Field")


subplot(4, 1, 3); 
plot(1:GD(3),squeeze(Potential(1,1,:)), "Linewidth", 2)
title("Electric Potential")
xlabel("Position (z-axis)")
ylabel("Voltage")

subplot(4, 1, 4); 
imagesc(squeeze(rot90(real(E2p(5,:,:))))); h = colorbar; set(get(h,'title'),'string','V/m');
title("Internal Electric Field Magnitude")



set(gcf, "Position", [100 100 600 700])




