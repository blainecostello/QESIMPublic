function MeasDatCompare3(runID)

% loop through all files in  ::RunID/Solved::

% Initialize Measured Data for Comparison:

%MolFrac = [0, 0.005, 0.02, 0.1];

VolFrac = [0.0,0.01,0.02,0.03,0.04];%MolFrac.*0.40208;

% for Gold (Au) inclusions in an Alumina (Al2O3) host
Minc = 196.96655;
Mhst = 101.96;
Pinc = 19300;
Phst = 3950;

%VolFrac = MolFrac./(MolFrac+(1-MolFrac).*((Mhst*Pinc)/(Minc*Phst)));

% Placeholder until bdfs data is obtained...
MBDFS   = [395.5, 335, 327, 266];
% nano CU in PMMA @ 20Hz
%Meff    = [5.7, 7.0, 7.65, 8.1, 8.8];
% micro CU in PMMA @ 20Hz
%Meff    = [5.7, 6, 6.3, 6.65, 7.1];
% micro CU in PMMA @ 464Hz
Meff    = [5.4, 5.7, 5.95, 6.3, 6.7];

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
    
    if  real(effSp) < 50 % && tanD_Sp < 1  && effComBdfs < 1E12 % breakdown is below vacuum and permittivity is reasonably small
        Seff(cuou-1)  = real(effSp);
        StanD(cuou-1) = tanD_Sp;
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


% figure
% hold on
% plot(VolFrac, MTanD, '-ok', 'linewidth', 2)
% plot(SVF, StanD, 'o', 'linewidth', 2)
% title("Loss Tangent")
% xlabel("Volume Fraction of Inclusions")
% ylabel("Loss Tangent")
% set(gcf, 'Position', [450 300 400 300])


% figure
% hold on
% plot(VolFrac, MBDFS, '-ok', 'linewidth', 2)
% plot(SVF, SBDFS, 'o', 'linewidth', 2)
% title("Breakdown Field Strength")
% xlabel("Volume Fraction of Inclusions")
% ylabel("Breakdown Strength")
% set(gcf, 'Position', [1000 300 400 300])



end

