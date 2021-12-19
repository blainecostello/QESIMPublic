function MeasDatCompPerc(runID)

% loop through all files in  ::RunID/Solved::

% Initialize Measured Data for Comparison:

%MolFrac = [0, 0.005, 0.02, 0.1];

VolFrac = [0.0,0.1,0.2,0.3];%MolFrac.*0.40208;

freq=10000
e0 = 8.854187812813E-12

actConds = [1.5E-9, 5.87E7]
actMats  = [2.5 - 1i*actConds(1)/(freq*e0), 1 - 1i*actConds(2)/(freq*e0)]


% for Gold (Au) inclusions in an Alumina (Al2O3) host
Minc = 196.96655;
Mhst = 101.96;
Pinc = 19300;
Phst = 3950;

%VolFrac = MolFrac./(MolFrac+(1-MolFrac).*((Mhst*Pinc)/(Minc*Phst)));

% Placeholder until bdfs data is obtained...
MBDFS   = [395.5, 335, 327, 266];
% CU in PMMA @ 10kHz
Meff    = [2.5, 5.7, 6.8, 7.9];
MTanD   = [0.013, 0.02, 0.025, 0.03];

% initialize arrays for plotting simulated data
Seff = zeros(size(Meff));
Seff(1) = Meff(1);

StanD = zeros(size(MTanD));
StanD(1) = MTanD(1);

SBDFS = zeros(size(MBDFS));
SBDFS(1) = MBDFS(1);

SVF = zeros(size(VolFrac));

% Iterate through files in /Solved/ folder
%runID = 2234;

dirname = "Run_" + runID + "/Solved/";
files = dir(dirname);
[numF,~] = size(files);

% Initialize scatterplot storage structures Nx3 
btVFr = zeros(numF,3);
btVFr = zeros(numF,3);
effSct = btVFr;
bdfsSct = btVFr;
lossSct = btVFr;
condSct = btVFr;


vfs = [0.03:0.03:0.3];
% initialize average line storage structures and counters
%(do this later...)
vfAve = zeros(length(vfs),4);
counter = vfAve



inv = 0;
alltheVF = zeros(1,numF);
for cuou = 3:numF
    Fn = cuou - inv ;
    fname = files(Fn).name;
    fname = [dirname+fname];
    %importfile(fname)
    load(fname)
    
    % Check if geometry is valid... (assume true unless we run into issues)
        if(btVF > 0) % && other stuff?
            
            % Get prox & categorize into perc (1), nearPerc (2), and nonPerc (3)
            if prox == 0
                p = 1;
            elseif prox == 1
                p = 2;
            elseif prox > 1 || prox < 0
                p = 3;
            else
                p = 4;
            end
            if p~=4      % Determine percolating configuration (1=P,2=NP,3=NoP)
                
                % get index associated with volume fraction
                volfracs = round(vfs*1000) == round(VF*1000);

                
                
                % Recalculate Permittivity 
                EinEa = ((1/2)*e0 .* real(actMats(Omega(:,:,1:GD(3)-1)+1)).* abs(E2p).^2) .* (Delta(1)*Delta(2)*Delta(3)); % energy in each voxel
                Ed = sum(EinEa, 3);
                Ed = sum(Ed, 2);
                EtotRept = sum(Ed, 1);
                effRe = (2*EtotRept*Delta(3)*(GD(3)-1))/((1^2)*e0*Delta(1)*GD(1)*Delta(2)*GD(2));    


                Seff(cuou-1)  = real(effRe);


                % Use Total Power Dissipated through the simulated volume
                powTot = sum(sum(sum(abs(E2p(:,:,:)).^2 .* (actConds(Omega(:,:,2:GD(3))+1)),3),2),1) * Delta(1)*Delta(2)*(Delta(3))
                % To calculate the effective conductivity of the simulated composite
                cond_fromPow = powTot / (Eapp^2*Delta(3)*GD(3))     % [S/m]

                % add in freq-dependent component...
                eff_recalc = real(effRe) - 1i*(real(cond_fromPow))/(freq * e0);
                
                % Calculate loss tangent 
                tanD_recalc = abs(imag(eff_recalc)/real(eff_recalc))


                % Assign everything...
                Seff(cuou-1) = real(eff_recalc)
                StanD(cuou-1) = tanD_recalc%tanD_recalc;
                SBDFS(cuou-1) = effComBdfs;
                SVF(cuou-1)   = btVF;
                
                
                % increment counter
                counter(volfracs,p) = counter(volfracs,p) + 1;

                % Assign New Stuff
                vfsSct(sum(counter(:,p)),p) = VF;
                effSct(sum(counter(:,p)),p) = real(eff_recalc);
                

                
                % 
                
                
                
                if length(btVF) < 2
                    alltheVF(Fn-1) = btVF;
                else
                    files(Fn) = [];
                    inv = inv + 1;
                end
            end
    else 
        disp("eff = " + real(effSp))
        disp("tanD = " + tanD_Sp)
        disp("bdfs = " + effComBdfs)
        disp("-----------------------")
    end
    
end
alltheVF = alltheVF(1:(numF-inv-1));
[~,VFsort]=sort(alltheVF); %Get the order of B
files = files(VFsort);




disp(cond_fromPow)
disp(eff_recalc)

figure
hold on
plot(VolFrac, Meff, '-ok', 'linewidth', 2)
plot(SVF(), effSct(:,1), 'o', 'linewidth', 2)
plot(SVF, effSct(:,1), 'o', 'linewidth', 2)
plot(SVF, effSct(:,1), 'o', 'linewidth', 2)

title("Effective Composite Permittivity")
xlabel("Volume Fraction of Inclusions")
ylabel("Permittivity")
legend("Measured Data", "Simulated Data", "Location", "Northwest")
set(gcf, 'Position', [0 300 300 200])


figure
hold on
plot(VolFrac, MTanD, '-ok', 'linewidth', 2)
plot(SVF, StanD, 'o', 'linewidth', 2)
title("Loss Tangent")
xlabel("Volume Fraction of Inclusions")
ylabel("Loss Tangent")
%ylim([0, 0.1])

legend("Measured Data", "Simulated Data", "Location", "Northwest")
set(gcf, 'Position', [300 300 300 200])



figure
hold on
plot(VolFrac, MBDFS, '-ok', 'linewidth', 2)
plot(SVF, SBDFS, 'o', 'linewidth', 2)
title("Breakdown Field Strength")
xlabel("Volume Fraction of Inclusions")
ylabel("Breakdown Strength")
%ylim([100, 450])
legend("Measured Data", "Simulated Data", "Location","Northeast")
set(gcf, 'Position', [600 300 300 200])
ylim([0,500])



end


