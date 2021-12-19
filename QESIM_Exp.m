function QESIM_Exp(gd, seed, volfracs, npDiams, tSamp, ovlp, Perms, Conds, bdfs, freq, runID)

% Simulate cubic region of nanoparticle composite material subject to
% parameters defined in the function argument

warning on verbose
warning on backtrace


% Parameters for Experiment:
%npDiams = 7E-9
% 20 um
%tSamp   = 100E-9       % Half sample thickness to capture adequate resolution of NPs...

% Verify below params for specific papers...
% Perms   = [4.35,   1]
% Conds   = [1E-8, 2.3E7]
% freq    = 1E6 % 1MHZ

%runID = 2000+gd



% Sample Run: QESIM(25, 6.5E-9, 50E-9, 27, 1155)

% GD     - length of single edge of cubic simulation window
% npDiam - diameter of a single nanoparticle
% tSamp  - Sample Thickness
% seed   - seed for pseudorandomization of microstructure
% runID  - determines save folder location 

% Initialize Wt% Vector:
% WtFrac = [0.05, 0.1];

% for Gold (Au) inclusions in an SU-8) host
% Pinc = 19300;   % Au
% Phst = 1200;    % SU-8

% Convert to Volume Fraction
% VolFrac = WtFrac./(WtFrac+(1-WtFrac).*((Pinc)/(Phst)));

% MolFrac = [0.005, 0.02, 0.1];
% 
% % VolFrac = MolFrac.*0.40208;
% 
% % for Gold (Au) inclusions in an Alumina (Al2O3) host
% Minc = 196.96655;
% Mhst = 101.96;
% Pinc = 19300;
% Phst = 3950;
% 
% VolFrac = MolFrac./(MolFrac+(1-MolFrac).*((Mhst*Pinc)/(Minc*Phst)));

wtPrc = [0.05,  0.1];

rhoAu = 19300;
rhoSu = 1200;

e0       = 8.8541878128E-12;

VolFrac = volfracs;% wtPrc./(wtPrc+(1-wtPrc).*(rhoAu/rhoSu));

% Initialize Material Parameters for Host & Inclusion materials
    sigmaH   = Conds(1);   
    relPermH = Perms(1);   
    sigmaI   = Conds(2);         % 1/Ohm-Meters
    relPermI = Perms(2);
    if(length(Conds) > 2)
        sigmaInt = Conds(3);          % Conductivity of Interphase
        rPermInt = Perms(3);             % Interphase Real permittivity
    end
    
    
    % Calculate Complex Permittivity of all constituent phases.
if(length(Conds) > 2)
    Mats     = [relPermH - 1i*(sigmaH/(freq*e0)), relPermI - 1i*(sigmaI/(freq*e0)), rPermInt - 1i*(sigmaInt/(freq*e0))];
    %bdfs     = [395.5    1E-12     1   ];
else
    
    Mats     = [relPermH - 1i*(sigmaH/(freq*e0)), relPermI - 1i*(sigmaI/(freq*e0))];
    %bdfs     = [395.5    1E-12  ];
end

% Define grid & cell dimensions
gd = [gd   gd   gd];
Delta = [tSamp/gd(1) tSamp/gd(2) tSamp/gd(3)];

% Initialize geometric nanoparticle parameters
rads    = [  0  ; npDiams'./(2*Delta(1))];  % Cubic Voxels Necessary for this to work...
shells  = [  0:length(npDiams)];



if(length(npDiams) == 1)
    % Host, Core, Shell1, Shell2,...
    bdfs    = bdfs(1:2);
    Mats    = Mats(1:2);
elseif(length(npDiams) == 2)
    bdfs    = bdfs(1:3);
    Mats    = Mats(1:3);
else
    disp("Incomplete Material Description: Please declare the breakdown of constituent phases")
    % Here....
end





%ovlp    = 2; 

% Initialize Applied Voltage in Z-Direction
%Volt = [-1,-1;-1,-1;-0.5,0.5]; 
Volt = [-1,-1;-1,-1;0,1]; 
%Volt = [-1,-1;-1,-1;1,2]; 



%==========================={ Run-Da-Traaap }=============================%


% Start Parallel Pool 
%ppool = parpool;

% Run all Mat-Gens (Time Each)
for(i = 1:length(VolFrac))
    tic
    MatGenTQE(gd, rads, shells, ovlp, VolFrac(i), runID, seed);
    toc
end

% Display Volume Fractions
disp(VolFrac);
%close all;

% Run all EM-Scans (Time Each)
for(i = 1:length(VolFrac))
    tic
    EMScanTQE(Volt, gd, Mats, bdfs, VolFrac(i), Delta, freq, runID, seed);
    toc
end

%delete(ppool)

%MeasDatCompare2(runID);




end