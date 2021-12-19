function MeasDatCompare(runID)

% loop through all files in  ::RunID/Solved::

% Initialize Measured Data for Comparison:

%MolFrac = [0, 0.005, 0.02, 0.1];

VolFrac = [0.0,0.1,0.2,0.3];%MolFrac.*0.40208;

actConds = [1.5E-9, 5.87E7]


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

inv = 0;
alltheVF = zeros(1,numF);
for cuou = 3:numF
    Fn = cuou - inv ;
    fname = files(Fn).name;
    fname = [dirname+fname];
    %importfile(fname)
    load(fname)
    
    if  1% real(effSp) < 1E6  && effComBdfs <  2000% breakdown is below vacuum and permittivity is reasonably small
        Seff(cuou-1)  = real(effSp);
        
        % recalculate total current [A/m^2]
        %jtot_recalc = sum(sum(sum(abs(Exyz(:,:,:,3)).*imag(Mats(Omega(:,:,1:GD(3)-1)+1))*e0*freq,3),2),1) ;
        jtot_recalc = sum(sum(sum(abs(E2p(:,:,:)).*imag(Mats(Omega(:,:,2:GD(3))+1))*e0*freq,3),2),1) ;% * ((GD(1)*Delta(1))*(GD(2)*(Delta(2))))/((GD(3)-1)*20)
        
        % recalculate conductivity from total current
        cond_recalc = (((jtot_recalc * ((GD(1)*Delta(1))*(GD(2)*(Delta(2))))/((GD(3)-1)))/Eapp))
        
        powTot = sum(sum(sum(abs(E2p(:,:,:)).^2 .* imag(Mats(Omega(:,:,2:GD(3))+1))*e0*freq,3),2),1) * Delta(1)*Delta(2)*Delta(3)
        
        cond_fromPow = powTot / (Eapp^2*Delta(1)*Delta(2)*GD(1)*GD(2)*GD(3)*1000)
        
        % recalculate effective permittivity
        %eff_recalc = real(effSp) - 1i*(real(cond_recalc))/(freq * e0)
        eff_recalc = real(effSp) - 1i*(real(cond_fromPow))/(freq * e0)
        
        % reccalculate loss tangent 
        tanD_recalc = abs(imag(eff_recalc)/real(eff_recalc))
        
        %Seff(cuou-1) = real(eff_recalc)
        StanD(cuou-1) = tanD_recalc%tanD_recalc;
        SBDFS(cuou-1) = effComBdfs;
        SVF(cuou-1)   = btVF;

        if length(btVF) < 2
            alltheVF(Fn-1) = btVF;
        else
            files(Fn) = [];
            inv = inv + 1;
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







figure
hold on
plot(VolFrac, Meff, '-ok', 'linewidth', 2)
plot(SVF, Seff, 'o', 'linewidth', 2)
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
set(gcf, 'Position', [450 300 300 200])



figure
hold on
plot(VolFrac, MBDFS, '-ok', 'linewidth', 2)
plot(SVF, SBDFS, 'o', 'linewidth', 2)
title("Breakdown Field Strength")
xlabel("Volume Fraction of Inclusions")
ylabel("Breakdown Strength")
%ylim([100, 450])
legend("Measured Data", "Simulated Data", "Location","Northeast")
set(gcf, 'Position', [1000 300 300 200])




end

