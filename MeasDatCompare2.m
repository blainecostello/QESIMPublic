function MeasDatCompare2(runID)

% loop through all files in  ::RunID/Solved::

% Initialize Measured Data for Comparison:

wtPrc = [0,   0.05,  0.1  ];

rhoAu = 19300;
rhoSu = 1200;

e0       = 8.8541878128E-12;

VolFrac = wtPrc./(wtPrc+(1-wtPrc).*(rhoAu/rhoSu));

tanD1 = [0.033, 0.01,  0.095];
tanD2 = [0.014, 0.022, 0.070];
tanD3 = [0.015, 0.029, 0.055];

effP1 = [4.35,  5.55 , 13.75];
effP2 = [4.45,  5.20,  13.30];
effP3 = [4.35,  5.15,  12.75];


freq1 = 1E3;
freq2 = 2E3;
freq3 = 1E4;

w1 = 2*pi*freq1;
w2 = 2*pi*freq2;
w3 = 2*pi*freq3;

cond1 = tanD1(1) * effP1(1) * e0 * w1;
cond2 = tanD2(1) * effP2(1) * e0 * w2;
cond3 = tanD3(1) * effP3(1) * e0 * w3;

% initialize arrays for plotting simulated data
Seff = zeros(size(effP1));
Seff(1) = effP1(1);

StanD = zeros(size(tanD1));
StanD(1) = tanD1(1);


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
    
    if  real(effSp) < 100 && tanD_Sp < 10  && effComBdfs < 1E12 % breakdown is below vacuum and permittivity is reasonably small
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



if length(Mats) == 3
    ShellCond = -imag(Mats(3)) * e0*freq;
    ShellPerm = real(Mats(3));
else
    ShellCond = -imag(Mats(2)) * e0*freq;
    ShellPerm = real(Mats(2));
end

figure
hold on
plot(VolFrac, effP1, '-ok', 'linewidth', 2)
plot(SVF, Seff, 'o', 'linewidth', 2)
sc = sprintf('%2.0e|', ShellCond);
sp = sprintf('%2.0e|', ShellPerm);

title("Shell Cond = " + sc + "   |  Perm = " + sp)
xlabel("Volume Fraction of Inclusions")
ylabel("Permittivity")
legend("Measured Data", "Simulated Data", "Location", "Northwest")
set(gcf, 'Position', [0 300 300 200])


figure
hold on
plot(VolFrac, tanD1, '-ok', 'linewidth', 2)
plot(SVF, StanD, 'o', 'linewidth', 2)
title("Shell Cond = " + sc + "   |  Perm = " + sp)
xlabel("Volume Fraction of Inclusions")
ylabel("Loss Tangent")
set(gcf, 'Position', [1000 300 300 200])




end

