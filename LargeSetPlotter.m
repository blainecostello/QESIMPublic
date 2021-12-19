%% Unpack & Plot Large Datasets
clc
clear all
PSeeds = [3,31,40,47];
NPseeds = [12,23,24,37,53,59,62,67,68,71,77,79];
therest = [13,14,18,41,48,55,66,73];
%BTbd = 1e6;
dirnm = 'Run_12702/Solved/';

% perrydat = [0,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6; 11.9, 12.5, 15.2, 16.4, 18.6, 22.3, 26.2, 35.1, 37.5]
% perrybd = [0,0.05,0.1,0.15,0.2,0.3,0.4,0.5;380,344,345,300,235,230,218,213];
% perryed = [0,0.1,0.2,0.3,0.4,0.5;1.6E6,1.6E6,2E6,2.3E6,2.35E6,3.2E6];
perm = figure;
hold on
efld = figure;
hold on
eden = figure;
hold on
bdfs = figure;
hold on
efldBT = figure;
hold on
countbars = figure;
hold on
perm2 = figure;
hold on
efld2 = figure;
hold on
eden2 = figure;
hold on
bdfs2 = figure;
hold on
efldBT2 = figure;
hold on
volts = figure;
hold on
volts2 = figure;
hold on


Em2p = (1)/(200e-9)
vfs = [0.005:0.005:0.05]

vl = length(vfs)
files = dir(dirnm);
effchart = zeros(vl,100,4);
avEchart = effchart;
avENPchart = effchart;
bdchart = effchart;
brugchart = effchart;
maxEdchart = effchart;
[numF,~] = size(files);
numF = numF -2;
files = files(3:numF+2);
vfAve = zeros(1,length(vfs),4);
effAve = vfAve;
aveMaxED = vfAve;
counter = vfAve+1;
countvf = vfAve;
aveE = vfAve;
aveENP = vfAve;
avBD = vfAve;


clear VFrac
clear eff
clear btVFr
clear bdv2pBT
clear BTstdE
clear maxVolt2pNP
clear MaxEd2pNP
inv = 0;
alltheVF = zeros(1,numF);
for j = 2:numF
    Fn = j - inv ;
    fnm = files(Fn).name;
    fnm = [dirnm,fnm];
    importfile(fnm)
    if length(btVF) < 2
        alltheVF(Fn-1) = btVF;
    else
        files(Fn) = [];
        inv = inv + 1;
    end
end
alltheVF = alltheVF(1:(numF-inv-1))
[~,VFsort]=sort(alltheVF); %Get the order of B
files = files(VFsort);

brg = zeros(1,length(vfs))
btVFr = zeros(numF,3);
btVFsct = btVFr;
effSct = btVFr;
avESct = btVFr;
avENPSct = btVFr;
bdSct = btVFr;
maxEdSct = btVFr;
effSctS = btVFr;
avESctS = btVFr;
avENPSctS = btVFr;
bdSctS = btVFr;
maxEdSctS = btVFr;
              
              
for j = 2:numF-1
    i = j;
    fnm = files(i).name;
    fnm = [dirnm,fnm];
    if length(fnm) > 30
        importfile(fnm)
        
        if(btVF > 0)
            
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
            if p~=4
                % Scale Things...
                eff2p = real(effSp);
%                 Emean2p = Emean2p;
%                 Emean2pNP = Emean2pNP;
                PVDFbd = PVDFbd*1e-6;
                BTbd = BTbd*1e-6;
                % ^^ fix later in Tester3 so that variation (stdev) is 
                % correctly calculated
                
                
                eff(j,1,p) = eff2p;
                eff(j,2,p) = eff3p;
                eff(j,3,p) = eff3pH;
                eff(j,4,p) = eff3pS;
                % use path to compute average electric field along breakdown path
                stE2pNP = std(std(E2p(OmegaBT(:,:,1:gd-1)==2))); %.*(Delta(1)*Delta(2)*Delta(3));
                bdv2pBT(j,p) = (gd*BTbd)./((Emean2pNP+stE2p));
                maxVolt2pNP(j,p) =   BTbd ./ ( Emean2pNP + stE2pNP );
                MaxEd2pNP(j,p) = Edens2p * maxVolt2pNP(j,p).^2;
                MaxVolt2p(j,p) =  PVDFbd ./ ( Emean2p + stE2p );
                avE(j,1,p) = avE2p;
                
                Ed = (E2p.^2 .*(1/2)*e0 .* Mats(OmegaBT(:,:,1:Zsz-1)+1)); % Delta(1)*Delta(2)*Delta(3)
                Ed = sum(Ed,3);
                Ed = sum(Ed,2);
                Ed2p = sum(Ed,1);
                Edens2p(vfi) = Ed2p/(gd*gd*gd);
                
                % mcE(vfs == volfrac) = mcE
                % maxEd(j) = MaxEd2p;
                if ((maxVolt2pNP(j,p) < MaxVolt2p(j,p)) && p~=3)
                    maxEd(j,p) = (maxVolt2pNP(j,p)*1E6)^2*Edens2p;
                    % FIX THIS --->> maxEdS(j,p) = (maxVolt3pNPS(j,p)*1E6)^2*Edens3pS;
                    bd = maxVolt2pNP(j,p)/(Delta(3)*gd);%BTbd/((Emean2pNP+stE2pNP)/(3.1E6)); % Approximate y intercept
%                     bdS = maxVolt2pNPS(j,p)/(Delta(3)*gd);%BTbd/((Emean2pNP+stE2pNP)/(3.1E6)); % Approximate y intercept
                     % Account for breakdown voltage inside nanoparticle
                     % HERE
                else
                    maxEd(j,p) = (MaxVolt2p(j,p)*1E6)^2*Edens2p;
%                     maxEdS(j,p) = (maxVolt3pS(j,p)*1E6)^2*Edens3pS;
                    bd = MaxVolt2p(j,p)/(Delta(3)*gd);%PVDFbd/((Emean2p+stE2p)/(7.8E6)); % divide by approximate y intercept...
%                     bdS = MaxVolt3pS(j,p)/(Delta(3)*gd);%PVDFbd/((Emean2p+stE2p)/(7.8E6)); % divide by approximate y intercept...

                    % AND HERE
                end
                volfracs = round(vfs*1000) == round(volfrac*1000);
                disp(volfracs)
                eff2p = eff2p/2; % Adjust permittivity
                aveE(1,volfracs,p) = aveE(1,volfracs,p) + Emean2p;
                aveENP(1,volfracs,p) = aveENP(1,volfracs,p) + Emean2pNP;
                counter(1,volfracs,p) = counter(1,volfracs,p) + 1;
                effAve(1,volfracs,p) = effAve(1,volfracs,p) + eff2p;
                avBD(1,volfracs,p) = avBD(1,volfracs,p) + bd;
                aveMaxED(1,volfracs,p) = aveMaxED(1,volfracs,p) + maxEd(j,p);
                effchart(volfracs, counter(1,volfracs,p),p) = eff2p;
                avEchart(volfracs, counter(1,volfracs,p),p) = Emean2p;
                avENPchart(volfracs, counter(1,volfracs,p),p) = Emean2pNP;
                bdchart(volfracs, counter(1,volfracs,p),p) = bd;
                maxEdchart(volfracs, counter(1,volfracs,p),p) = maxEd(j,p);
%                 
%                 eff3pS = eff3pS/2; % Adjust permittivity
%                 aveES(1,volfracs,p) = aveES(1,volfracs,p) + Emean3pS;
%                 aveENPS(1,volfracs,p) = aveENPS(1,volfracs,p) + Emean3pNPS;
%                 counterS(1,volfracs,p) = counterS(1,volfracs,p) + 1;
%                 effAveS(1,volfracs,p) = effAveS(1,volfracs,p) + eff3pS;
%                 avBDS(1,volfracs,p) = avBDS(1,volfracs,p) + bd;
%                 aveMaxEDS(1,volfracs,p) = aveMaxEDS(1,volfracs,p) + maxEdS(j,p);
%                 effchartS(volfracs, counterS(1,volfracs,p),p) = eff3pS;
%                 avEchartS(volfracs, counterS(1,volfracs,p),p) = Emean3pS;
%                 avENPchartS(volfracs, counterS(1,volfracs,p),p) = Emean3pNPS;
%                 bdchartS(volfracs, counterS(1,volfracs,p),p) = bdS;
%                 maxEdchartS(volfracs, counterS(1,volfracs,p),p) = maxEdS(j,p);
                
                vfchart(volfracs, counter(1,volfracs,p),p) = btVF;
                
                cSct = sum(counter-1,2);
                btVFr(   cSct(p),p) = btVF;
                btVFsct( cSct(p),p) = btVF;
                effSct(  cSct(p),p) = eff2p;
                avESct(  cSct(p),p) = Emean2p;
                avENPSct(cSct(p),p) = Emean2pNP;
                bdSct(   cSct(p),p) = bd;
                maxEdSct(j,p) = maxEd(j,p);
                                
%                 cSctS = sum(counterS-1,2)
%                 btVFrS(   cSct(p),p) = btVF;
%                 btVFsctS( cSct(p),p) = btVF;
%                 effSctS(  cSctS(p),p) = eff3pS;
%                 avESctS(  cSctS(p),p) = Emean3pS;
%                 avENPSctS(cSctS(p),p) = Emean3pNPS;
%                 bdSctS(   cSctS(p),p) = bdS;
%                 maxEdSctS(j,p) = maxEdS(j,p);
                
                
                vfAve(1,volfracs,p) = vfAve(1,volfracs,p) + btVF;
                vfAve(1,volfracs,4) = vfAve(1,volfracs,4) + btVF;
%                 vfAveS(1,volfracs,p) = vfAveS(1,volfracs,p) + btVF;
%                 vfAveS(1,volfracs,4) = vfAveS(1,volfracs,4) + btVF;
                countvf(1,volfracs,p) = countvf(1,volfracs,p) + 1;
                
                countvf(1,volfracs,4) = countvf(1,volfracs,4) + 1;
                counter(1,volfracs,4) = counter(1,volfracs,4) + 1;
%                 effchartS(volfracs, counter(1,volfracs,4),4) = eff2p;
%                 avEchartS(volfracs, counter(1,volfracs,4),4) = Emean2p;
%                 avENPchartS(volfracs, counter(1,volfracs,4),4) = Emean2pNP;
%                 bdchartS(volfracs, counter(1,volfracs,4),4) = bd;
%                 maxEdchartS(volfracs, counter(1,volfracs,4),4) = maxEd(j,p);
                
                                syms x % effective permittivity of solid
            if((phases==2) && (brg(mci,volfracs)==0))
                brgEq = (1-volfrac)*((CxPvty(1)-x)/(CxPvty(1)+x)) + (volfrac)*((CxPvty(3)-x)/(CxPvty(3)+x)) == 0;
            elseif((phases == 3)  && (brg(mci,volfracs)==0))
                brgEq = (1-volfrac)*((CxPvty(1)-x)/(CxPvty(1)+x)) + (volfrac)*((CxPvty(2)-x)/(CxPvty(2)+x)) + (volfrac2)*((CxPvty(3)-x)/(CxPvty(3)+x))==0;
            end
            soln = solve(brgEq,x);
            brg(mci,volfracs) = max(soln);


            end
        end
    end
end

for bcnt = 1:length(volfracs)
    brgEq = (1-vfs(bcnt))*((Mats(1)-x)/(Mats(1)+x)) + (vfs(bcnt))*((Mats(3)-x)/(Mats(3)+x)) == 0;
    soln = solve(brgEq,x);
    brugg = max(soln);
    brug(bcnt) = brugg;
end
%     figure(perm);
%     plot(vfs, brug, '--', 'linewidth',2);
    
% line identifiers
linetps = ['-o';'-x';'-*';'-v'];
dottps = ['o';'x';'*';'v'];
for i = 1:3
    ctemp = nonzeros(counter(1,:,i)-1)
    l = length(ctemp)
    
    effAve(1,1:l,i) = nonzeros(effAve(1,:,i))./ctemp;
    effAve(1,l+1:vl,i) = 0;
    aveE(1,1:l,i) = nonzeros(aveE(1,:,i))./ctemp;
    aveE(1,l+1:vl,i) = 0;
    aveENP(1,1:l,i) = nonzeros(aveENP(1,:,i))./ctemp;
    aveENP(1,l+1:vl,i) = 0;
    avBD(1,1:l,i) = nonzeros(avBD(1,:,i))./ctemp;
    avBD(1,l+1:vl,i) = 0;
    vfAve(1,1:l,i)  = nonzeros(vfAve(1,:,i))./ctemp;
    vfAve(1,l+1:vl,i) = 0;
    aveMaxED(1,1:l,i) = nonzeros(aveMaxED(1,:,i))./ctemp;
    aveMaxED(1,l+1:vl,i) = 0;
    for num = 1:length(vfs)
        stdev(num,i) = std(nonzeros(effchart(num,1:nonzeros(counter(1,num,i)-1),i)));
        stdevE(num,i) = std(nonzeros(avEchart(num,1:nonzeros(counter(1,num,i)-1),i)));
        stdevENP(num,i) = std(nonzeros(avENPchart(num,1:nonzeros(counter(1,num,i)-1),i)));
        stdevBD(num,i) = std(nonzeros(bdchart(num,1:nonzeros(counter(1,num,i)-1),i)))
        stdevMED(num,i) = std(nonzeros(maxEdchart(num,1:nonzeros(counter(1,num,i)-1),i)));
        
        stdev(num,4) = std(nonzeros(effchart(num,1:nonzeros(counter(1,num,4)-1),4)));
        stdevE(num,4) = std(nonzeros(avEchart(num,1:nonzeros(counter(1,num,4)-1),4)));
        stdevENP(num,4) = std(nonzeros(avENPchart(num,1:nonzeros(counter(1,num,4)-1),4)));
        stdevBD(num,4) = std(nonzeros(bdchart(num,1:nonzeros(counter(1,num,4)-1),4)))
        stdevMED(num,4) = std(nonzeros(maxEdchart(num,1:nonzeros(counter(1,num,4)-1),4)));
        
%         stdevS(num,i) = std(nonzeros(effchartS(num,1:nonzeros(counterS(1,num,i)-1),i)));
%         stdevES(num,i) = std(nonzeros(avEchartS(num,1:nonzeros(counterS(1,num,i)-1),i)));
%         stdevENPS(num,i) = std(nonzeros(avENPchartS(num,1:nonzeros(counterS(1,num,i)-1),i)));
%         stdevBDS(num,i) = std(nonzeros(bdchartS(num,1:nonzeros(counterS(1,num,i)-1),i)))
%         stdevMEDS(num,i) = std(nonzeros(maxEdchartS(num,1:nonzeros(counterS(1,num,i)-1),i)));
%         
%         stdevS(num,4) = std(nonzeros(effchartS(num,1:nonzeros(counterS(1,num,4)-1),4)));
%         stdevES(num,4) = std(nonzeros(avEchartS(num,1:nonzeros(counterS(1,num,4)-1),4)));
%         stdevENPS(num,4) = std(nonzeros(avENPchartS(num,1:nonzeros(counterS(1,num,4)-1),4)));
%         stdevBDS(num,4) = std(nonzeros(bdchartS(num,1:nonzeros(counterS(1,num,4)-1),4)))
%         stdevMEDS(num,4) = std(nonzeros(maxEdchartS(num,1:nonzeros(counterS(1,num,4)-1),4)));
        

    end
    
%     figure(perm);
%     errorbar(nonzeros(vfAve(:,:,i)), nonzeros(effAve(1,:,i))', 3*stdev(counter(1,:,i)-1~=0,i), dottps(i), 'linewidth',1.5);
%     
%     figure(efld)
%     errorbar(nonzeros(vfAve(:,:,i)), nonzeros(aveE(1,:,i))', 3*stdevE(counter(1,:,i)-1~=0,i), linetps(i), 'linewidth',1.5);
% 
%     figure(efld2)
%     errorbar(nonzeros(vfAve(:,:,i)), nonzeros(aveE(1,:,i))', 3*stdevE(counter(1,:,i)-1~=0,i), dottps(i), 'linewidth',1.5);
%     
%     figure(efldBT)
%     errorbar(nonzeros(vfAve(:,:,i)), nonzeros(aveENP(1,:,i))', 3*stdevENP(counter(1,:,i)-1~=0,i), dottps(i), 'linewidth',1.5);
%     
%     figure(eden);
%     errorbar(nonzeros(vfAve(:,:,i)), nonzeros(aveMaxED(1,:,i))', 3*stdevMED(counter(1,:,i)-1~=0,i), dottps(i),'linewidth',1.5);
%     
%     figure(bdfs);
%     errorbar(nonzeros(vfAve(:,:,i)), nonzeros(avBD(1,:,i))', stdevBD(counter(1,:,i)-1~=0,i), dottps(i),'linewidth',1.5);
%     
end
ctemp = nonzeros(counter(1,:,4)-1);
l = length(ctemp)

ctemp = nonzeros(counter(1,:,4)-1);
effAve(1,1:l,4) = nonzeros(sum(sum(effchart(:,:,1:3),3),2))./ctemp;
effAve(1,l+1:vl,4) = 0;
aveE(1,1:l,4) = nonzeros(sum(sum(avEchart(:,:,1:3),2),3))./ctemp;
aveE(1,l+1:vl,4) = 0;
aveENP(1,1:l,4) = nonzeros(sum(sum(avENPchart(:,:,1:3),2),3))./ctemp;
aveENP(1,l+1:vl,4) = 0;
avBD(1,1:l,4) = nonzeros(sum(sum(bdchart(:,:,1:3),2),3))./ctemp;
avBD(1,l+1:vl,4) = 0;
aveMaxED(1,1:l,4) = nonzeros(sum(sum(maxEdchart(:,:,1:3),2),3))./ctemp;
aveMaxED(1,l+1:vl,4) = 0;

vfAve(1,1:l,4)  = nonzeros(vfAve(1,:,4))./ctemp;
vfAve(1,l+1:vl,4) = 0;

counter(1,:,4) = sum(counter(1,:,:)-1,3);

%------------------for loop ended here ------------%
%% EFFECTIVE PERMITTIVITY PLOTS
figure(perm)
c = get(gca,'colororder') % variation of permittivity as scatterplot
scatter(btVFsct(:,1),effSct(:,1),[],c(1,:),'.');
scatter(btVFsct(:,2),effSct(:,2),[],c(1,:),'o');
scatter(btVFsct(:,3),effSct(:,3),[],c(1,:),'x');
errorbar(nonzeros(vfAve(:,:,4)), nonzeros(effAve(1,:,4))', 3*stdev(counter(1,:,4)~=0,4), '-k', 'linewidth',1);
ylabel("Effective Relative Permittivity [F/M]");
xlabel("Volume Fraction of Nanoparticles [%]");
legend('Percolating Geometries','Near Percolating Geometries', 'Non-Percolating Geometries', 'All Geometries', 'Location', 'southeast')

figure(perm2)% mean of permittivity with perry and bruggeman
plot(nonzeros(vfAve(:,:,4)), nonzeros(effAve(1,:,4))', '-k', 'linewidth',1.5);
plot(perrydat(1,:), perrydat(2,:),'ok','linewidth',1.5)
plot(vfs, brug,'--k','linewidth',1.5)
ylabel("Effective Relative Permittivity [F/M]");
xlabel("Volume Fraction of Nanoparticles [%]");
legend('Mean of tested geometries','Experimental Data (Perry 2009)', 'Bruggeman Model', 'All Geometries', 'Location', 'southeast')



%% ENERGY DENSITY  PLOTS

figure(eden)
errorbar(nonzeros(vfAve(:,:,4)), nonzeros(aveMaxED(1,:,4))', 3*stdevMED(counter(1,:,4)~=0,4), dottps(4), 'linewidth',2);

ylabel("Maximum Energy Density [J/cm^3]");
xlabel("Volume Fraction of Nanoparticles [%]");
legend('Percolating Geometries','Near Percolating Geometries', 'Non-Percolating Geometries', 'All Geometries', 'Location', 'southeast')

figure(eden2)
%plot(nonzeros(vfAve(:,:,4)), nonzeros(aveMaxED(1,:,4))', linetps(4), 'linewidth',1.5);
c = get(gca,'colororder')
scatter(nonzeros(btVFsct(:,1)),nonzeros(maxEdSct(:,1)),[],c(1,:),'.');
scatter(nonzeros(btVFsct(:,2)),nonzeros(maxEdSct(:,2)),[],c(1,:),'o');
scatter(nonzeros(btVFsct(:,3)),nonzeros(maxEdSct(:,3)),[],c(1,:),'x');
plot(vfAve(:,:,4),aveMaxED(:,:,4),'-k','linewidth',1.5) % Plot mean line
%plot(perryed(1,:), perryed(2,:),'.-.','linewidth',1.5)
ylabel("Maximum Energy Density [J/cm^3]");
xlabel("Volume Fraction of Nanoparticles [%]");
legend('Percolating Geometries','Near Percolating Geometries', 'Non-Percolating Geometries', 'Mean', 'Location', 'southeast')
 
%% BREAKDOWN FIELD STRENGTH PLOTS
figure(bdfs)
c = get(gca,'colororder') % variation of permittivity as scatterplot
scatter(nonzeros(btVFsct(:,3)),nonzeros(bdSct(:,3)).*1E-3,[],c(1,:),'.');
scatter(nonzeros(btVFsct(:,2)),nonzeros(bdSct(:,2)).*1E-3,[],c(1,:),'o');
scatter(nonzeros(btVFsct(:,1)),nonzeros(bdSct(:,1)).*1E-3,[],c(1,:),'x');

errorbar(nonzeros(vfAve(:,:,4)), nonzeros(avBD(1,:,4))'.*1E-3, 3*stdevBD(counter(1,:,4)~=0,4).*1E-3, '-k', 'linewidth',1);
%plot(perrybd(1,:), perrybd(2,:),'ok','linewidth',1.5)
ylabel("Breakdown Field Strength [MV/M]");
xlabel("Volume Fraction of Nanoparticles [%]");
legend('Percolating Geometries','Near Percolating Geometries', 'Non-Percolating Geometries', 'All Geometries','Measured Data (Perry)', 'Location', 'northeast')

figure(bdfs2)
plot(nonzeros(vfAve(:,:,4)), nonzeros(avBD(1,:,4))', '-k', 'linewidth',1.5);
plot(perrybd(1,:), perrybd(2,:),'ok','linewidth',1.5)
ylabel("Breakdown Field Strength [MV/M]");
xlabel("Volume Fraction of Nanoparticles[%]");
legend('Simulated Geometries','Measured Data (Perry)', 'Non-Percolating Geometries', 'All Geometries', 'Location', 'northeast')


%% MEAN E-FIELD PLOTS
figure(efld)
scatter(btVFsct(:,1),avESct(:,1)/Em2p,[],c(1,:),'.');
scatter(btVFsct(:,2),avESct(:,2)/Em2p,[],c(1,:),'o');
scatter(btVFsct(:,3),avESct(:,3)/Em2p,[],c(1,:),'x');
scatter(btVFsct(:,1),avENPSct(:,1)/Em2p,[],c(1,:),'.');
scatter(btVFsct(:,2),avENPSct(:,2)/Em2p,[],c(1,:),'o');
scatter(btVFsct(:,3),avENPSct(:,3)/Em2p,[],c(1,:),'x');
errorbar(nonzeros(vfAve(:,:,4)), nonzeros(aveE(1,:,4)/Em2p)', 3*stdevE(counter(1,:,4)~=0,4)/Em2p, '-k', 'linewidth',1);
errorbar(nonzeros(vfAve(:,:,4)), nonzeros(aveENP(1,:,4)/Em2p)', 3*stdevENP(counter(1,:,4)~=0,4)/Em2p, '-k', 'linewidth',1);
ylabel("Average Normalized Electric Field");
xlabel("Volume Fraction of Nanoparticles [%]");
%legend('Inside the Host Material','Inside the Nanoparticles', 'Location', 'southeast')

% figure(efldBT)
% scatter(btVFsct(:,1),avENPSct(:,1),[],c(1,:),'o');
% scatter(btVFsct(:,2),avENPSct(:,2),[],c(2,:),'o');
% scatter(btVFsct(:,3),avENPSct(:,3),[],c(3,:),'o');
% ylabel("Average Electric Field within Nanoparticles");
% xlabel("Volume Fraction of BT Nanoparticles");
% errorbar(nonzeros(vfAve(:,:,4)), nonzeros(aveENP(1,:,4))', stdevENP(counter(1,:,4)~=0,4), linetps(4), 'linewidth',2);

% legend('Percolating Geometries','Near Percolating Geometries', 'Non-Percolating Geometries', 'All Geometries', 'Location', 'southeast')
%% OTHER MISC. PLOTS
% max V vs. permittivity
figure(volts)
%errorbar(nonzeros(vfAve(:,:,4)), nonzeros(avBD(1,:,4))', 3*stdevBD(counter(1,:,4)~=0,4), linetps(4), 'linewidth',1);
c = get(gca,'colororder') % variation of permittivity as scatterplot
scatter(nonzeros(btVFsct(:,1)),nonzeros(MaxVolt2p(:,1)),[],c(1,:),'.');
scatter(nonzeros(btVFsct(:,2)),nonzeros(MaxVolt2p(:,2)),[],c(1,:),'o');
scatter(nonzeros(btVFsct(:,3)),nonzeros(MaxVolt2p(:,3)),[],c(1,:),'x');

scatter(nonzeros(btVFsct(:,1)),nonzeros(maxVolt2pNP(:,1)),[],c(1,:),'.');
scatter(nonzeros(btVFsct(:,2)),nonzeros(maxVolt2pNP(:,2)),[],c(1,:),'o');
scatter(nonzeros(btVFsct(:,3)),nonzeros(maxVolt2pNP(:,3)),[],c(1,:),'x');
ylabel("Maximum Voltage [V]");
xlabel("Volume Fraction of Nanoparticles [%]");
legend('Percolating in Host','Near Percolating in Host', 'Non-Percolating in Host', 'Percolating in NP','Near Percolating in NP', 'Non-Percolating in NP', 'Location', 'southeast')

figure(volts2)
%errorbar(nonzeros(vfAve(:,:,4)), nonzeros(avBD(1,:,4))', 3*stdevBD(counter(1,:,4)~=0,4), linetps(4), 'linewidth',1);
c = get(gca,'colororder') % variation of permittivity as scatterplot
scatter(nonzeros(effSct(:,1)),nonzeros(MaxVolt2p(:,1)),[],c(1,:),'o');
scatter(nonzeros(effSct(:,2)),nonzeros(MaxVolt2p(:,2)),[],c(1,:),'o');
scatter(nonzeros(effSct(:,3)),nonzeros(MaxVolt2p(:,3)),[],c(1,:),'o');

scatter(nonzeros(effSct(:,1)),nonzeros(maxVolt2pNP(:,1)),[],c(1,:),'.');
scatter(nonzeros(effSct(:,2)),nonzeros(maxVolt2pNP(:,2)),[],c(1,:),'.');
scatter(nonzeros(effSct(:,3)),nonzeros(maxVolt2pNP(:,3)),[],c(1,:),'.');
ylabel("Maximum Voltage [V]");
xlabel("Effective Permittivity [F/M]");
legend('Percolating in Host','Near Percolating in Host', 'Non-Percolating in Host','Percolating in NP','Near Percolating in NP', 'Non-Percolating in NP', 'Location', 'northeast')


%% COUNTER BAR PLOT
figure(countbars)
allcounts = [counter(:,:,1)'-1, counter(:,:,2)'-1, counter(:,:,3)'-1];% counter(:,:,4)'];
%hb = bar(vfAve(:,:,4)', allcounts, 'stacked');
ylabel("Number of Microstructures");
xlabel("Volume Fraction of Nanoparticles [%]");
ylim([0,80]);
legend('Percolaging Geometries', 'Near Percolating Geometries', 'Non-Percolating Geometries')
saveas(countbars,'bars2','fig');
saveas(countbars,'bars2','jpeg');
saveas(perm,'permValid2','fig');
saveas(perm,'permValid2','jpeg');
saveas(eden,'eden2','fig');
saveas(eden,'eden2','jpeg');
saveas(efldBT,'efldBT2','fig');
saveas(efldBT,'efldBT2','jpeg');
saveas(bdfs,'bdfs2','fig');
saveas(bdfs,'bdfs2','jpeg');
saveas(efld,'efld_hst2','fig');
saveas(efld,'efld_hst2','jpeg');
saveas(perm2,'perm3','fig');
saveas(perm2,'perm3','jpeg');
saveas(eden2,'eden3','fig');
saveas(eden2,'eden3','jpeg');
saveas(efldBT2,'efldBT3','fig');
saveas(efldBT2,'efldBT3','jpeg');
saveas(bdfs2,'bdfs3','fig');
saveas(bdfs2,'bdfs3','jpeg');
saveas(efld2,'efld_hst3','fig');
saveas(efld2,'efld_hst3','jpeg');
saveas(volts,'volts1','fig');
saveas(volts,'volts1','jpeg');
saveas(volts2,'volts2','fig');
saveas(volts2,'volts2','jpeg');

MCP = figure; 
hold on
%MCPdat = [3.9,3.9,6.3,6.7,7.1,7.1,8.5,8.5,21.0,21.0,21.0,19.0,19.0,23.4,26,26,26.8,29,31,60,95,95,95,95,50,98;1400,1310,630,870,760,920,620,720,470,400,300,250,600,350,320,420,420,380,400,320,300,250,200,100,330]
Mcpdat = [3.9,1400;3.9,1310;6.3,630;6.7,870;7.1,760;7.1,920;8.5,620;8.5,720;21.0,470;21.0,400;21.0,300;19.0,250;19.0,600;23.4,350;26,320;26,420;26.8,420;29,380;31,400;60,320;95,300;95,250;95,200;95,100;50,300;98,330;]
scatter(effSct(:,3),bdSct(:,3).^2,[],c(1,:),'.');
scatter(effSct(:,2),bdSct(:,2).^2,[],c(1,:),'.');
scatter(effSct(:,1),bdSct(:,1).^2,[],c(1,:),'.');
c = get(gca,'colororder') % variation of permittivity as scatterplot

plot(Mcpdat(:,1), Mcpdat(:,2).^2, 'ok', 'linewidth',1);
ylabel("Breakdown Field Strength^2 [MV/M]");
xlabel("Effective Relative Permittivity [F/M]");
legend('Percolating Geometries','Near Percolating Geometries', 'Non-Percolating Geometries', 'Data from ', 'Location', 'southeast')
set(gca,'xscale','log')
set(gca,'yscale','log')

% Get a list of all variables
allvars = whos;
% Identify the variables that ARE NOT graphics handles. This uses a regular
% expression on the class of each variable to check if it's a graphics object
tosave = cellfun(@isempty, regexp({allvars.class}, '^matlab\.(ui|graphics)\.'));
% Pass these variable names to save
save('dataset2_nf.mat', allvars(tosave).name)


%% Unpack & Plot Large Datasets
% clc
% clear all
PSeeds = [3,31,40,47];
NPseeds = [12,23,24,37,53,59,62,67,68,71,77,79];
therest = [13,14,18,41,48,55,66,73];
BTbd = 1e6;
dirnm = 'PercRuns/Run_1031/';

perrydat = [0,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6; 11.9, 12.5, 15.2, 16.4, 18.6, 22.3, 26.2, 35.1, 37.5]
perrybd = [0,0.05,0.1,0.15,0.2,0.3,0.4,0.5;380,344,345,300,235,230,218,213];
perryed = [0,0.1,0.2,0.3,0.4,0.5;1.6E6,1.6E6,2E6,2.3E6,2.35E6,3.2E6];
% perm = figure;
% hold on
% efld = figure;
% hold on
% eden = figure;
% hold on
% bdfs = figure;
% hold on
% efldBT = figure;
% hold on
% countbars = figure;
% hold on
% perm2 = figure;
% hold on
% efld2 = figure;
% hold on
% eden2 = figure;
% hold on
% bdfs2 = figure;
% hold on
% efldBT2 = figure;
% hold on
% volts = figure;
% hold on
% volts2 = figure;
% hold on


Em2p = (1)/(65*2e-9)
vfs = [0.05:0.05:0.6]

vl = length(vfs)
files = dir(dirnm);
effchart = zeros(vl,100,4);
avEchart = effchart;
avENPchart = effchart;
bdchart = effchart;
brugchart = effchart;
maxEdchart = effchart;
[numF,~] = size(files);
numF = numF -2;
files = files(3:numF+2);
vfAve = zeros(1,length(vfs),4);
effAve = vfAve;
aveMaxED = vfAve;
counter = vfAve+1;
countvf = vfAve;
aveE = vfAve;
aveENP = vfAve;
avBD = vfAve;


clear VFrac
clear eff
clear btVFr
clear bdv2pBT
clear BTstdE
clear MaxVolt2p
clear maxVolt2pNP
clear MaxEd2pNP
inv = 0;
alltheVF = zeros(1,numF);
for j = 2:numF
    Fn = j - inv ;
    fnm = files(Fn).name;
    fnm = [dirnm,fnm];
    importfile(fnm)
    if length(btVF) < 2
        alltheVF(Fn-1) = btVF;
    else
        files(Fn) = [];
        inv = inv + 1;
    end
end
alltheVF = alltheVF(1:(numF-inv-1))
[~,VFsort]=sort(alltheVF); %Get the order of B
files = files(VFsort);

brg = zeros(1,length(vfs))
btVFr = zeros(numF,3);
btVFsct = btVFr;
effSct = btVFr;
avESct = btVFr;
avENPSct = btVFr;
bdSct = btVFr;
maxEdSct = btVFr;
effSctS = btVFr;
avESctS = btVFr;
avENPSctS = btVFr;
bdSctS = btVFr;
maxEdSctS = btVFr;
              
              
for j = 2:numF-1
    i = j;
    fnm = files(i).name;
    fnm = [dirnm,fnm];
    if length(fnm) > 30
        importfile(fnm)
        
        if(btVF > 0)
            
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
            if p~=4
                % Scale Things...
                eff2p = eff2p*2;
%                 Emean2p = Emean2p;
%                 Emean2pNP = Emean2pNP;
                PVDFbd = PVDFbd*1e-6;
                BTbd = BTbd*1e-6;
                % ^^ fix later in Tester3 so that variation (stdev) is 
                % correctly calculated
                
                
                eff(j,1,p) = eff2p;
                eff(j,2,p) = eff3p;
                eff(j,3,p) = eff3pH;
                eff(j,4,p) = eff3pS;
                % use path to compute average electric field along breakdown path
                stE2pNP = std(std(std(E2p(OmegaBT(:,:,1:gd-1)==2)))); %.*(Delta(1)*Delta(2)*Delta(3));
                bdv2pBT(j,p) = (gd*BTbd)./((Emean2pNP+stE2p));
                maxVolt2pNP(j,p) =   BTbd ./ ( Emean2pNP + stE2pNP );
                MaxEd2pNP(j,p) = Edens2p * maxVolt2pNP(j,p).^2;
                MaxVolt2p(j,p) =  PVDFbd ./ ( Emean2p + stE2p );
                avE(j,1,p) = avE2p;
                
                Ed = (E2p.^2 .*(1/2)*e0 .* Mats(OmegaBT(:,:,1:Zsz-1)+1)); % Delta(1)*Delta(2)*Delta(3)
                Ed = sum(Ed,3);
                Ed = sum(Ed,2);
                Ed2p = sum(Ed,1);
                Edens2p(vfi) = Ed2p/(gd*gd*gd);
                
                % mcE(vfs == volfrac) = mcE
                % maxEd(j) = MaxEd2p;
                if ((maxVolt2pNP(j,p) < MaxVolt2p(j,p)) && p~=3)
                    maxEd(j,p) = (maxVolt2pNP(j,p)*1E6)^2*Edens2p;
                    % FIX THIS --->> maxEdS(j,p) = (maxVolt3pNPS(j,p)*1E6)^2*Edens3pS;
                    bd = maxVolt2pNP(j,p)/(Delta(3)*gd);%BTbd/((Emean2pNP+stE2pNP)/(3.1E6)); % Approximate y intercept
%                     bdS = maxVolt2pNPS(j,p)/(Delta(3)*gd);%BTbd/((Emean2pNP+stE2pNP)/(3.1E6)); % Approximate y intercept
                     % Account for breakdown voltage inside nanoparticle
                     % HERE
                else
                    maxEd(j,p) = (MaxVolt2p(j,p)*1E6)^2*Edens2p;
%                     maxEdS(j,p) = (maxVolt3pS(j,p)*1E6)^2*Edens3pS;
                    bd = MaxVolt2p(j,p)/(Delta(3)*gd);%PVDFbd/((Emean2p+stE2p)/(7.8E6)); % divide by approximate y intercept...
%                     bdS = MaxVolt3pS(j,p)/(Delta(3)*gd);%PVDFbd/((Emean2p+stE2p)/(7.8E6)); % divide by approximate y intercept...

                    % AND HERE
                end
                volfracs = round(vfs*1000) == round(volfrac*1000);
                disp(volfracs)
                eff2p = eff2p/2; % Adjust permittivity
                aveE(1,volfracs,p) = aveE(1,volfracs,p) + Emean2p;
                aveENP(1,volfracs,p) = aveENP(1,volfracs,p) + Emean2pNP;
                counter(1,volfracs,p) = counter(1,volfracs,p) + 1;
                effAve(1,volfracs,p) = effAve(1,volfracs,p) + eff2p;
                avBD(1,volfracs,p) = avBD(1,volfracs,p) + bd;
                aveMaxED(1,volfracs,p) = aveMaxED(1,volfracs,p) + maxEd(j,p);
                effchart(volfracs, counter(1,volfracs,p),p) = eff2p;
                avEchart(volfracs, counter(1,volfracs,p),p) = Emean2p;
                avENPchart(volfracs, counter(1,volfracs,p),p) = Emean2pNP;
                bdchart(volfracs, counter(1,volfracs,p),p) = bd;
                maxEdchart(volfracs, counter(1,volfracs,p),p) = maxEd(j,p);
%                 
%                 eff3pS = eff3pS/2; % Adjust permittivity
%                 aveES(1,volfracs,p) = aveES(1,volfracs,p) + Emean3pS;
%                 aveENPS(1,volfracs,p) = aveENPS(1,volfracs,p) + Emean3pNPS;
%                 counterS(1,volfracs,p) = counterS(1,volfracs,p) + 1;
%                 effAveS(1,volfracs,p) = effAveS(1,volfracs,p) + eff3pS;
%                 avBDS(1,volfracs,p) = avBDS(1,volfracs,p) + bd;
%                 aveMaxEDS(1,volfracs,p) = aveMaxEDS(1,volfracs,p) + maxEdS(j,p);
%                 effchartS(volfracs, counterS(1,volfracs,p),p) = eff3pS;
%                 avEchartS(volfracs, counterS(1,volfracs,p),p) = Emean3pS;
%                 avENPchartS(volfracs, counterS(1,volfracs,p),p) = Emean3pNPS;
%                 bdchartS(volfracs, counterS(1,volfracs,p),p) = bdS;
%                 maxEdchartS(volfracs, counterS(1,volfracs,p),p) = maxEdS(j,p);
                
                vfchart(volfracs, counter(1,volfracs,p),p) = btVF;
                
                cSct = sum(counter-1,2)
                btVFr(   cSct(p),p) = btVF;
                btVFsct( cSct(p),p) = btVF;
                effSct(  cSct(p),p) = eff2p;
                avESct(  cSct(p),p) = Emean2p;
                avENPSct(cSct(p),p) = Emean2pNP;
                bdSct(   cSct(p),p) = bd;
                maxEdSct(j,p) = maxEd(j,p);
                                
%                 cSctS = sum(counterS-1,2)
%                 btVFrS(   cSct(p),p) = btVF;
%                 btVFsctS( cSct(p),p) = btVF;
%                 effSctS(  cSctS(p),p) = eff3pS;
%                 avESctS(  cSctS(p),p) = Emean3pS;
%                 avENPSctS(cSctS(p),p) = Emean3pNPS;
%                 bdSctS(   cSctS(p),p) = bdS;
%                 maxEdSctS(j,p) = maxEdS(j,p);
                
                
                vfAve(1,volfracs,p) = vfAve(1,volfracs,p) + btVF;
                vfAve(1,volfracs,4) = vfAve(1,volfracs,4) + btVF;
%                 vfAveS(1,volfracs,p) = vfAveS(1,volfracs,p) + btVF;
%                 vfAveS(1,volfracs,4) = vfAveS(1,volfracs,4) + btVF;
                countvf(1,volfracs,p) = countvf(1,volfracs,p) + 1;
                
                countvf(1,volfracs,4) = countvf(1,volfracs,4) + 1;
                counter(1,volfracs,4) = counter(1,volfracs,4) + 1;
%                 effchartS(volfracs, counter(1,volfracs,4),4) = eff2p;
%                 avEchartS(volfracs, counter(1,volfracs,4),4) = Emean2p;
%                 avENPchartS(volfracs, counter(1,volfracs,4),4) = Emean2pNP;
%                 bdchartS(volfracs, counter(1,volfracs,4),4) = bd;
%                 maxEdchartS(volfracs, counter(1,volfracs,4),4) = maxEd(j,p);
                
                                syms x % effective permittivity of solid
            if((phases==2) && (brg(mci,volfracs)==0))
                brgEq = (1-volfrac)*((CxPvty(1)-x)/(CxPvty(1)+x)) + (volfrac)*((CxPvty(3)-x)/(CxPvty(3)+x)) == 0;
            elseif((phases == 3)  && (brg(mci,volfracs)==0))
                brgEq = (1-volfrac)*((CxPvty(1)-x)/(CxPvty(1)+x)) + (volfrac)*((CxPvty(2)-x)/(CxPvty(2)+x)) + (volfrac2)*((CxPvty(3)-x)/(CxPvty(3)+x))==0;
            end
            soln = solve(brgEq,x);
            brg(mci,volfracs) = max(soln);


            end
        end
    end
end

for bcnt = 1:length(volfracs)
    brgEq = (1-vfs(bcnt))*((Mats(1)-x)/(Mats(1)+x)) + (vfs(bcnt))*((Mats(3)-x)/(Mats(3)+x)) == 0;
    soln = solve(brgEq,x);
    brugg = max(soln);
    brug(bcnt) = brugg;
end
%     figure(perm);
%     plot(vfs, brug, '--', 'linewidth',2);
    
% line identifiers
linetps = ['-o';'-x';'-*';'-v'];
dottps = ['o';'x';'*';'v'];
for i = 1:3
    ctemp = nonzeros(counter(1,:,i)-1)
    l = length(ctemp)
    
    effAve(1,1:l,i) = nonzeros(effAve(1,:,i))./ctemp;
    effAve(1,l+1:vl,i) = 0;
    aveE(1,1:l,i) = nonzeros(aveE(1,:,i))./ctemp;
    aveE(1,l+1:vl,i) = 0;
    aveENP(1,1:l,i) = nonzeros(aveENP(1,:,i))./ctemp;
    aveENP(1,l+1:vl,i) = 0;
    avBD(1,1:l,i) = nonzeros(avBD(1,:,i))./ctemp;
    avBD(1,l+1:vl,i) = 0;
    vfAve(1,1:l,i)  = nonzeros(vfAve(1,:,i))./ctemp;
    vfAve(1,l+1:vl,i) = 0;
    aveMaxED(1,1:l,i) = nonzeros(aveMaxED(1,:,i))./ctemp;
    aveMaxED(1,l+1:vl,i) = 0;
    for num = 1:length(vfs)
        stdev(num,i) = std(nonzeros(effchart(num,1:nonzeros(counter(1,num,i)-1),i)));
        stdevE(num,i) = std(nonzeros(avEchart(num,1:nonzeros(counter(1,num,i)-1),i)));
        stdevENP(num,i) = std(nonzeros(avENPchart(num,1:nonzeros(counter(1,num,i)-1),i)));
        stdevBD(num,i) = std(nonzeros(bdchart(num,1:nonzeros(counter(1,num,i)-1),i)))
        stdevMED(num,i) = std(nonzeros(maxEdchart(num,1:nonzeros(counter(1,num,i)-1),i)));
        
        stdev(num,4) = std(nonzeros(effchart(num,1:nonzeros(counter(1,num,4)-1),4)));
        stdevE(num,4) = std(nonzeros(avEchart(num,1:nonzeros(counter(1,num,4)-1),4)));
        stdevENP(num,4) = std(nonzeros(avENPchart(num,1:nonzeros(counter(1,num,4)-1),4)));
        stdevBD(num,4) = std(nonzeros(bdchart(num,1:nonzeros(counter(1,num,4)-1),4)))
        stdevMED(num,4) = std(nonzeros(maxEdchart(num,1:nonzeros(counter(1,num,4)-1),4)));
        
%         stdevS(num,i) = std(nonzeros(effchartS(num,1:nonzeros(counterS(1,num,i)-1),i)));
%         stdevES(num,i) = std(nonzeros(avEchartS(num,1:nonzeros(counterS(1,num,i)-1),i)));
%         stdevENPS(num,i) = std(nonzeros(avENPchartS(num,1:nonzeros(counterS(1,num,i)-1),i)));
%         stdevBDS(num,i) = std(nonzeros(bdchartS(num,1:nonzeros(counterS(1,num,i)-1),i)))
%         stdevMEDS(num,i) = std(nonzeros(maxEdchartS(num,1:nonzeros(counterS(1,num,i)-1),i)));
%         
%         stdevS(num,4) = std(nonzeros(effchartS(num,1:nonzeros(counterS(1,num,4)-1),4)));
%         stdevES(num,4) = std(nonzeros(avEchartS(num,1:nonzeros(counterS(1,num,4)-1),4)));
%         stdevENPS(num,4) = std(nonzeros(avENPchartS(num,1:nonzeros(counterS(1,num,4)-1),4)));
%         stdevBDS(num,4) = std(nonzeros(bdchartS(num,1:nonzeros(counterS(1,num,4)-1),4)))
%         stdevMEDS(num,4) = std(nonzeros(maxEdchartS(num,1:nonzeros(counterS(1,num,4)-1),4)));
        

    end
    
%     figure(perm);
%     errorbar(nonzeros(vfAve(:,:,i)), nonzeros(effAve(1,:,i))', 3*stdev(counter(1,:,i)-1~=0,i), dottps(i), 'linewidth',1.5);
%     
%     figure(efld)
%     errorbar(nonzeros(vfAve(:,:,i)), nonzeros(aveE(1,:,i))', 3*stdevE(counter(1,:,i)-1~=0,i), linetps(i), 'linewidth',1.5);
% 
%     figure(efld2)
%     errorbar(nonzeros(vfAve(:,:,i)), nonzeros(aveE(1,:,i))', 3*stdevE(counter(1,:,i)-1~=0,i), dottps(i), 'linewidth',1.5);
%     
%     figure(efldBT)
%     errorbar(nonzeros(vfAve(:,:,i)), nonzeros(aveENP(1,:,i))', 3*stdevENP(counter(1,:,i)-1~=0,i), dottps(i), 'linewidth',1.5);
%     
%     figure(eden);
%     errorbar(nonzeros(vfAve(:,:,i)), nonzeros(aveMaxED(1,:,i))', 3*stdevMED(counter(1,:,i)-1~=0,i), dottps(i),'linewidth',1.5);
%     
%     figure(bdfs);
%     errorbar(nonzeros(vfAve(:,:,i)), nonzeros(avBD(1,:,i))', stdevBD(counter(1,:,i)-1~=0,i), dottps(i),'linewidth',1.5);
%     
end
ctemp = nonzeros(counter(1,:,4)-1);
l = length(ctemp)

ctemp = nonzeros(counter(1,:,4)-1);
effAve(1,1:l,4) = nonzeros(sum(sum(effchart(:,:,1:3),3),2))./ctemp;
effAve(1,l+1:vl,4) = 0;
aveE(1,1:l,4) = nonzeros(sum(sum(avEchart(:,:,1:3),2),3))./ctemp;
aveE(1,l+1:vl,4) = 0;
aveENP(1,1:l,4) = nonzeros(sum(sum(avENPchart(:,:,1:3),2),3))./ctemp;
aveENP(1,l+1:vl,4) = 0;
avBD(1,1:l,4) = nonzeros(sum(sum(bdchart(:,:,1:3),2),3))./ctemp;
avBD(1,l+1:vl,4) = 0;
aveMaxED(1,1:l,4) = nonzeros(sum(sum(maxEdchart(:,:,1:3),2),3))./ctemp;
aveMaxED(1,l+1:vl,4) = 0;

vfAve(1,1:l,4)  = nonzeros(vfAve(1,:,4))./ctemp;
vfAve(1,l+1:vl,4) = 0;

counter(1,:,4) = sum(counter(1,:,:)-1,3);

%------------------for loop ended here ------------%
%% EFFECTIVE PERMITTIVITY PLOTS
figure(perm)
c = get(gca,'colororder') % variation of permittivity as scatterplot
scatter(btVFsct(:,1),effSct(:,1),[],c(2,:),'.');
scatter(btVFsct(:,2),effSct(:,2),[],c(2,:),'o');
scatter(btVFsct(:,3),effSct(:,3),[],c(2,:),'x');
errorbar(nonzeros(vfAve(:,:,4)), nonzeros(effAve(1,:,4))', 3*stdev(counter(1,:,4)~=0,4), '-r', 'linewidth',1);
ylabel("Effective Relative Permittivity [F/M]");
xlabel("Volume Fraction of Nanoparticles [%]");
legend('Percolating Geometries','Near Percolating Geometries', 'Non-Percolating Geometries', 'All Geometries', 'Location', 'southeast')

figure(perm2)% mean of permittivity with perry and bruggeman
plot(nonzeros(vfAve(:,:,4)), nonzeros(effAve(1,:,4))', '-r', 'linewidth',1.5);
%plot(perrydat(1,:), perrydat(2,:),'ok','linewidth',1.5)
%plot(vfs, brug,'--k','linewidth',1.5)
ylabel("Effective Relative Permittivity [F/M]");
xlabel("Volume Fraction of Nanoparticles [%]");
legend('Mean of tested geometries','Experimental Data (Perry 2009)', 'Bruggeman Model', 'All Geometries', 'Location', 'southeast')



%% ENERGY DENSITY  PLOTS

figure(eden)
errorbar(nonzeros(vfAve(:,:,4)), nonzeros(aveMaxED(1,:,4))', 3*stdevMED(counter(1,:,4)~=0,4), dottps(4), 'linewidth',2);

ylabel("Maximum Energy Density [J/cm^3]");
xlabel("Volume Fraction of Nanoparticles [%]");
legend('Percolating Geometries','Near Percolating Geometries', 'Non-Percolating Geometries', 'All Geometries', 'Location', 'southeast')

figure(eden2)
%plot(nonzeros(vfAve(:,:,4)), nonzeros(aveMaxED(1,:,4))', linetps(4), 'linewidth',1.5);
c = get(gca,'colororder')
scatter(nonzeros(btVFsct(:,1)),nonzeros(maxEdSct(:,1)),[],c(2,:),'.');
scatter(nonzeros(btVFsct(:,2)),nonzeros(maxEdSct(:,2)),[],c(2,:),'o');
scatter(nonzeros(btVFsct(:,3)),nonzeros(maxEdSct(:,3)),[],c(2,:),'x');
plot(vfAve(:,:,4),aveMaxED(:,:,4),'.-.k','linewidth',1.5) % Plot mean line
%plot(perryed(1,:), perryed(2,:),'.-.','linewidth',1.5)
ylabel("Maximum Energy Density [J/cm^3]");
xlabel("Volume Fraction of Nanoparticles [%]");
legend('Percolating Geometries','Near Percolating Geometries', 'Non-Percolating Geometries', 'Mean', 'Location', 'southeast')
 
%% BREAKDOWN FIELD STRENGTH PLOTS
figure(bdfs)
c = get(gca,'colororder'); % variation of permittivity as scatterplot
scatter(nonzeros(btVFsct(:,3)),nonzeros(bdSct(:,3)).*1E-3,[],c(2,:),'.');
scatter(nonzeros(btVFsct(:,2)),nonzeros(bdSct(:,2)).*1E-3,[],c(2,:),'o');
scatter(nonzeros(btVFsct(:,1)),nonzeros(bdSct(:,1)).*1E-3,[],c(2,:),'x');

errorbar(nonzeros(vfAve(:,:,4)), nonzeros(avBD(1,:,4))'.*1E-3, 3*stdevBD(counter(1,:,4)~=0,4).*1E-3, '.-.k', 'linewidth',1);
%plot(perrybd(1,:), perrybd(2,:),'ok','linewidth',1.5)
ylabel("Breakdown Field Strength [MV/M]");
xlabel("Volume Fraction of Nanoparticles [%]");
legend('Percolating Geometries','Near Percolating Geometries', 'Non-Percolating Geometries', 'All Geometries','Measured Data (Perry)', 'Location', 'northeast')

figure(bdfs2)
plot(nonzeros(vfAve(:,:,4)), nonzeros(avBD(1,:,4))', '.-.k', 'linewidth',1.5);
plot(perrybd(1,:), perrybd(2,:),'ok','linewidth',1.5)
ylabel("Breakdown Field Strength [MV/M]");
xlabel("Volume Fraction of Nanoparticles[%]");
legend('Simulated Geometries','Measured Data (Perry)', 'Non-Percolating Geometries', 'All Geometries', 'Location', 'northeast')


%% MEAN E-FIELD PLOTS
figure(efld)
scatter(btVFsct(:,1),avESct(:,1)/Em2p,[],c(2,:),'.');
scatter(btVFsct(:,2),avESct(:,2)/Em2p,[],c(2,:),'o');
scatter(btVFsct(:,3),avESct(:,3)/Em2p,[],c(2,:),'x');
scatter(btVFsct(:,1),avENPSct(:,1)/Em2p,[],c(2,:),'.');
scatter(btVFsct(:,2),avENPSct(:,2)/Em2p,[],c(2,:),'o');
scatter(btVFsct(:,3),avENPSct(:,3)/Em2p,[],c(2,:),'x');
errorbar(nonzeros(vfAve(:,:,4)), nonzeros(aveE(1,:,4)/Em2p)', 3*stdevE(counter(1,:,4)~=0,4)/Em2p, '.-.k', 'linewidth',1);
errorbar(nonzeros(vfAve(:,:,4)), nonzeros(aveENP(1,:,4)/Em2p)', 3*stdevENP(counter(1,:,4)~=0,4)/Em2p, '.-.k', 'linewidth',1);
ylabel("Average Normalized Electric Field");
xlabel("Volume Fraction of Nanoparticles [%]");
%legend('Inside the Host Material','Inside the Nanoparticles', 'Location', 'southeast')

% figure(efldBT)
% scatter(btVFsct(:,1),avENPSct(:,1),[],c(1,:),'o');
% scatter(btVFsct(:,2),avENPSct(:,2),[],c(2,:),'o');
% scatter(btVFsct(:,3),avENPSct(:,3),[],c(3,:),'o');
% ylabel("Average Electric Field within Nanoparticles");
% xlabel("Volume Fraction of BT Nanoparticles");
% errorbar(nonzeros(vfAve(:,:,4)), nonzeros(aveENP(1,:,4))', stdevENP(counter(1,:,4)~=0,4), linetps(4), 'linewidth',2);

% legend('Percolating Geometries','Near Percolating Geometries', 'Non-Percolating Geometries', 'All Geometries', 'Location', 'southeast')
%% OTHER MISC. PLOTS
% max V vs. permittivity
figure(volts)
%errorbar(nonzeros(vfAve(:,:,4)), nonzeros(avBD(1,:,4))', 3*stdevBD(counter(1,:,4)~=0,4), linetps(4), 'linewidth',1);
% c = get(gca,'colororder') % variation of permittivity as scatterplot
scatter(nonzeros(btVFsct(:,1)),nonzeros(MaxVolt2p(:,1)),[],c(2,:),'o');
scatter(nonzeros(btVFsct(:,2)),nonzeros(MaxVolt2p(:,2)),[],c(2,:),'o');
scatter(nonzeros(btVFsct(:,3)),nonzeros(MaxVolt2p(:,3)),[],c(2,:),'o');

scatter(nonzeros(btVFsct(:,1)),nonzeros(maxVolt2pNP(:,1)),[],c(2,:),'.');
scatter(nonzeros(btVFsct(:,2)),nonzeros(maxVolt2pNP(:,2)),[],c(2,:),'.');
scatter(nonzeros(btVFsct(:,3)),nonzeros(maxVolt2pNP(:,3)),[],c(2,:),'.');
ylabel("Maximum Voltage [V]");
xlabel("Volume Fraction of Nanoparticles [%]");
legend('Percolating in Host','Near Percolating in Host', 'Non-Percolating in Host', 'Percolating in NP','Near Percolating in NP', 'Non-Percolating in NP', 'Location', 'southeast')

figure(volts2)
%errorbar(nonzeros(vfAve(:,:,4)), nonzeros(avBD(1,:,4))', 3*stdevBD(counter(1,:,4)~=0,4), linetps(4), 'linewidth',1);
c = get(gca,'colororder') % variation of permittivity as scatterplot
scatter(nonzeros(effSct(:,1)),nonzeros(MaxVolt2p(:,1)),[],c(2,:),'o');
scatter(nonzeros(effSct(:,2)),nonzeros(MaxVolt2p(:,2)),[],c(2,:),'o');
scatter(nonzeros(effSct(:,3)),nonzeros(MaxVolt2p(:,3)),[],c(2,:),'o');

scatter(nonzeros(effSct(:,1)),nonzeros(maxVolt2pNP(:,1)),[],c(2,:),'.');
scatter(nonzeros(effSct(:,2)),nonzeros(maxVolt2pNP(:,2)),[],c(2,:),'.');
scatter(nonzeros(effSct(:,3)),nonzeros(maxVolt2pNP(:,3)),[],c(2,:),'.');
ylabel("Maximum Voltage [V]");
xlabel("Effective Permittivity [F/M]");
legend('Percolating in Host','Near Percolating in Host', 'Non-Percolating in Host','Percolating in NP','Near Percolating in NP', 'Non-Percolating in NP', 'Location', 'northeast')


%% COUNTER BAR PLOT
figure(countbars)
allcounts = [counter(:,:,1)'-1, counter(:,:,2)'-1, counter(:,:,3)'-1];% counter(:,:,4)'];
%hb = bar(vfAve(:,:,4)', allcounts, 'stacked');
ylabel("Number of Microstructures");
xlabel("Volume Fraction of Nanoparticles [%]");
ylim([0,80]);
legend('Percolaging Geometries', 'Near Percolating Geometries', 'Non-Percolating Geometries')
saveas(countbars,'bars2','fig');
saveas(countbars,'bars2','jpeg');
saveas(perm,'permValid2','fig');
saveas(perm,'permValid2','jpeg');
saveas(eden,'eden2','fig');
saveas(eden,'eden2','jpeg');
saveas(efldBT,'efldBT2','fig');
saveas(efldBT,'efldBT2','jpeg');
saveas(bdfs,'bdfs2','fig');
saveas(bdfs,'bdfs2','jpeg');
saveas(efld,'efld_hst2','fig');
saveas(efld,'efld_hst2','jpeg');
saveas(perm2,'perm3','fig');
saveas(perm2,'perm3','jpeg');
saveas(eden2,'eden3','fig');
saveas(eden2,'eden3','jpeg');
saveas(efldBT2,'efldBT3','fig');
saveas(efldBT2,'efldBT3','jpeg');
saveas(bdfs2,'bdfs3','fig');
saveas(bdfs2,'bdfs3','jpeg');
saveas(efld2,'efld_hst3','fig');
saveas(efld2,'efld_hst3','jpeg');
saveas(volts,'volts1','fig');
saveas(volts,'volts1','jpeg');
saveas(volts2,'volts2','fig');
saveas(volts2,'volts2','jpeg');

figure(MCP)
hold on
%MCPdat = [3.9,3.9,6.3,6.7,7.1,7.1,8.5,8.5,21.0,21.0,21.0,19.0,19.0,23.4,26,26,26.8,29,31,60,95,95,95,95,50,98;1400,1310,630,870,760,920,620,720,470,400,300,250,600,350,320,420,420,380,400,320,300,250,200,100,330]
Mcpdat = [3.9,1400;3.9,1310;6.3,630;6.7,870;7.1,760;7.1,920;8.5,620;8.5,720;21.0,470;21.0,400;21.0,300;19.0,250;19.0,600;23.4,350;26,320;26,420;26.8,420;29,380;31,400;60,320;95,300;95,250;95,200;95,100;50,300;98,330;]
scatter(effSct(:,3),bdSct(:,3).^2,[],c(2,:),'.');
scatter(effSct(:,2),bdSct(:,2).^2,[],c(2,:),'.');
scatter(effSct(:,1),bdSct(:,1).^2,[],c(2,:),'.');
c = get(gca,'colororder') % variation of permittivity as scatterplot

% plot(Mcpdat(:,1), Mcpdat(:,2).^2, 'ok', 'linewidth',1);
ylabel("Breakdown Field Strength^2 [MV/M]");
xlabel("Effective Relative Permittivity [F/M]");
% legend('Percolating Geometries','Near Percolating Geometries', 'Non-Percolating Geometries', 'Data from ', 'Location', 'southeast')
set(gca,'xscale','log')
set(gca,'yscale','log')

% Get a list of all variables
allvars = whos;
% Identify the variables that ARE NOT graphics handles. This uses a regular
% expression on the class of each variable to check if it's a graphics object
tosave = cellfun(@isempty, regexp({allvars.class}, '^matlab\.(ui|graphics)\.'));
% Pass these variable names to save
save('dataset2_nf.mat', allvars(tosave).name)



%% Unpack & Plot Large Datasets
% clc
% clear all
PSeeds = [3,31,40,47];
NPseeds = [12,23,24,37,53,59,62,67,68,71,77,79];
therest = [13,14,18,41,48,55,66,73];
BTbd = 1e6;
dirnm = 'PercRuns/Run_1032/';

perrydat = [0,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6; 11.9, 12.5, 15.2, 16.4, 18.6, 22.3, 26.2, 35.1, 37.5]
perrybd = [0,0.05,0.1,0.15,0.2,0.3,0.4,0.5;380,344,345,300,235,230,218,213];
perryed = [0,0.1,0.2,0.3,0.4,0.5;1.6E6,1.6E6,2E6,2.3E6,2.35E6,3.2E6];
% perm = figure;
% hold on
% efld = figure;
% hold on
% eden = figure;
% hold on
% bdfs = figure;
% hold on
% efldBT = figure;
% hold on
% countbars = figure;
% hold on
% perm2 = figure;
% hold on
% efld2 = figure;
% hold on
% eden2 = figure;
% hold on
% bdfs2 = figure;
% hold on
% efldBT2 = figure;
% hold on
% volts = figure;
% hold on
% volts2 = figure;
% hold on


Em2p = (1)/(65*2e-9)
vfs = [0.05:0.05:0.6]

vl = length(vfs)
files = dir(dirnm);
effchart = zeros(vl,100,4);
avEchart = effchart;
avENPchart = effchart;
bdchart = effchart;
brugchart = effchart;
maxEdchart = effchart;
[numF,~] = size(files);
numF = numF -2;
files = files(3:numF+2);
vfAve = zeros(1,length(vfs),4);
effAve = vfAve;
aveMaxED = vfAve;
counter = vfAve+1;
countvf = vfAve;
aveE = vfAve;
aveENP = vfAve;
avBD = vfAve;


clear VFrac
clear eff
clear btVFr
clear bdv2pBT
clear BTstdE
clear MaxVolt2p
clear maxVolt2pNP
clear MaxEd2pNP
inv = 0;
alltheVF = zeros(1,numF);
for j = 2:numF
    Fn = j - inv ;
    fnm = files(Fn).name;
    fnm = [dirnm,fnm];
    importfile(fnm)
    if length(btVF) < 2
        alltheVF(Fn-1) = btVF;
    else
        files(Fn) = [];
        inv = inv + 1;
    end
end
alltheVF = alltheVF(1:(numF-inv-1))
[~,VFsort]=sort(alltheVF); %Get the order of B
files = files(VFsort);

brg = zeros(1,length(vfs))
btVFr = zeros(numF,3);
btVFsct = btVFr;
effSct = btVFr;
avESct = btVFr;
avENPSct = btVFr;
bdSct = btVFr;
maxEdSct = btVFr;
effSctS = btVFr;
avESctS = btVFr;
avENPSctS = btVFr;
bdSctS = btVFr;
maxEdSctS = btVFr;
              
              
for j = 2:numF-1
    i = j;
    fnm = files(i).name;
    fnm = [dirnm,fnm];
    if length(fnm) > 30
        importfile(fnm)
        
        if(btVF > 0)
            
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
            if p~=4
                % Scale Things...
                eff2p = eff2p*2;
%                 Emean2p = Emean2p;
%                 Emean2pNP = Emean2pNP;
                PVDFbd = PVDFbd*1e-6;
                BTbd = BTbd*1e-6;
                % ^^ fix later in Tester3 so that variation (stdev) is 
                % correctly calculated
                
                
                eff(j,1,p) = eff2p;
                eff(j,2,p) = eff3p;
                eff(j,3,p) = eff3pH;
                eff(j,4,p) = eff3pS;
                % use path to compute average electric field along breakdown path
                stE2pNP = std(std(std(E2p(OmegaBT(:,:,1:gd-1)==2)))); %.*(Delta(1)*Delta(2)*Delta(3));
                bdv2pBT(j,p) = (gd*BTbd)./((Emean2pNP+stE2p));
                maxVolt2pNP(j,p) =   BTbd ./ ( Emean2pNP + stE2pNP );
                MaxEd2pNP(j,p) = Edens2p * maxVolt2pNP(j,p).^2;
                MaxVolt2p(j,p) =  PVDFbd ./ ( Emean2p + stE2p );
                avE(j,1,p) = avE2p;
                
                Ed = (E2p.^2 .*(1/2)*e0 .* Mats(OmegaBT(:,:,1:Zsz-1)+1)); % Delta(1)*Delta(2)*Delta(3)
                Ed = sum(Ed,3);
                Ed = sum(Ed,2);
                Ed2p = sum(Ed,1);
                Edens2p(vfi) = Ed2p/(gd*gd*gd);
                
                % mcE(vfs == volfrac) = mcE
                % maxEd(j) = MaxEd2p;
                if ((maxVolt2pNP(j,p) < MaxVolt2p(j,p)) && p~=3)
                    maxEd(j,p) = (maxVolt2pNP(j,p)*1E6)^2*Edens2p;
                    % FIX THIS --->> maxEdS(j,p) = (maxVolt3pNPS(j,p)*1E6)^2*Edens3pS;
                    bd = maxVolt2pNP(j,p)/(Delta(3)*gd);%BTbd/((Emean2pNP+stE2pNP)/(3.1E6)); % Approximate y intercept
%                     bdS = maxVolt2pNPS(j,p)/(Delta(3)*gd);%BTbd/((Emean2pNP+stE2pNP)/(3.1E6)); % Approximate y intercept
                     % Account for breakdown voltage inside nanoparticle
                     % HERE
                else
                    maxEd(j,p) = (MaxVolt2p(j,p)*1E6)^2*Edens2p;
%                     maxEdS(j,p) = (maxVolt3pS(j,p)*1E6)^2*Edens3pS;
                    bd = MaxVolt2p(j,p)/(Delta(3)*gd);%PVDFbd/((Emean2p+stE2p)/(7.8E6)); % divide by approximate y intercept...
%                     bdS = MaxVolt3pS(j,p)/(Delta(3)*gd);%PVDFbd/((Emean2p+stE2p)/(7.8E6)); % divide by approximate y intercept...

                    % AND HERE
                end
                volfracs = round(vfs*1000) == round(volfrac*1000);
                disp(volfracs)
                eff2p = eff2p/2; % Adjust permittivity
                aveE(1,volfracs,p) = aveE(1,volfracs,p) + Emean2p;
                aveENP(1,volfracs,p) = aveENP(1,volfracs,p) + Emean2pNP;
                counter(1,volfracs,p) = counter(1,volfracs,p) + 1;
                effAve(1,volfracs,p) = effAve(1,volfracs,p) + eff2p;
                avBD(1,volfracs,p) = avBD(1,volfracs,p) + bd;
                aveMaxED(1,volfracs,p) = aveMaxED(1,volfracs,p) + maxEd(j,p);
                effchart(volfracs, counter(1,volfracs,p),p) = eff2p;
                avEchart(volfracs, counter(1,volfracs,p),p) = Emean2p;
                avENPchart(volfracs, counter(1,volfracs,p),p) = Emean2pNP;
                bdchart(volfracs, counter(1,volfracs,p),p) = bd;
                maxEdchart(volfracs, counter(1,volfracs,p),p) = maxEd(j,p);
%                 
%                 eff3pS = eff3pS/2; % Adjust permittivity
%                 aveES(1,volfracs,p) = aveES(1,volfracs,p) + Emean3pS;
%                 aveENPS(1,volfracs,p) = aveENPS(1,volfracs,p) + Emean3pNPS;
%                 counterS(1,volfracs,p) = counterS(1,volfracs,p) + 1;
%                 effAveS(1,volfracs,p) = effAveS(1,volfracs,p) + eff3pS;
%                 avBDS(1,volfracs,p) = avBDS(1,volfracs,p) + bd;
%                 aveMaxEDS(1,volfracs,p) = aveMaxEDS(1,volfracs,p) + maxEdS(j,p);
%                 effchartS(volfracs, counterS(1,volfracs,p),p) = eff3pS;
%                 avEchartS(volfracs, counterS(1,volfracs,p),p) = Emean3pS;
%                 avENPchartS(volfracs, counterS(1,volfracs,p),p) = Emean3pNPS;
%                 bdchartS(volfracs, counterS(1,volfracs,p),p) = bdS;
%                 maxEdchartS(volfracs, counterS(1,volfracs,p),p) = maxEdS(j,p);
                
                vfchart(volfracs, counter(1,volfracs,p),p) = btVF;
                
                cSct = sum(counter-1,2)
                btVFr(   cSct(p),p) = btVF;
                btVFsct( cSct(p),p) = btVF;
                effSct(  cSct(p),p) = eff2p;
                avESct(  cSct(p),p) = Emean2p;
                avENPSct(cSct(p),p) = Emean2pNP;
                bdSct(   cSct(p),p) = bd;
                maxEdSct(j,p) = maxEd(j,p);
                                
%                 cSctS = sum(counterS-1,2)
%                 btVFrS(   cSct(p),p) = btVF;
%                 btVFsctS( cSct(p),p) = btVF;
%                 effSctS(  cSctS(p),p) = eff3pS;
%                 avESctS(  cSctS(p),p) = Emean3pS;
%                 avENPSctS(cSctS(p),p) = Emean3pNPS;
%                 bdSctS(   cSctS(p),p) = bdS;
%                 maxEdSctS(j,p) = maxEdS(j,p);
                
                
                vfAve(1,volfracs,p) = vfAve(1,volfracs,p) + btVF;
                vfAve(1,volfracs,4) = vfAve(1,volfracs,4) + btVF;
%                 vfAveS(1,volfracs,p) = vfAveS(1,volfracs,p) + btVF;
%                 vfAveS(1,volfracs,4) = vfAveS(1,volfracs,4) + btVF;
                countvf(1,volfracs,p) = countvf(1,volfracs,p) + 1;
                
                countvf(1,volfracs,4) = countvf(1,volfracs,4) + 1;
                counter(1,volfracs,4) = counter(1,volfracs,4) + 1;
%                 effchartS(volfracs, counter(1,volfracs,4),4) = eff2p;
%                 avEchartS(volfracs, counter(1,volfracs,4),4) = Emean2p;
%                 avENPchartS(volfracs, counter(1,volfracs,4),4) = Emean2pNP;
%                 bdchartS(volfracs, counter(1,volfracs,4),4) = bd;
%                 maxEdchartS(volfracs, counter(1,volfracs,4),4) = maxEd(j,p);
                
                                syms x % effective permittivity of solid
            if((phases==2) && (brg(mci,volfracs)==0))
                brgEq = (1-volfrac)*((CxPvty(1)-x)/(CxPvty(1)+x)) + (volfrac)*((CxPvty(3)-x)/(CxPvty(3)+x)) == 0;
            elseif((phases == 3)  && (brg(mci,volfracs)==0))
                brgEq = (1-volfrac)*((CxPvty(1)-x)/(CxPvty(1)+x)) + (volfrac)*((CxPvty(2)-x)/(CxPvty(2)+x)) + (volfrac2)*((CxPvty(3)-x)/(CxPvty(3)+x))==0;
            end
            soln = solve(brgEq,x);
            brg(mci,volfracs) = max(soln);


            end
        end
    end
end

for bcnt = 1:length(volfracs)
    brgEq = (1-vfs(bcnt))*((Mats(1)-x)/(Mats(1)+x)) + (vfs(bcnt))*((Mats(3)-x)/(Mats(3)+x)) == 0;
    soln = solve(brgEq,x);
    brugg = max(soln);
    brug(bcnt) = brugg;
end
%     figure(perm);
%     plot(vfs, brug, '--', 'linewidth',2);
    
% line identifiers
linetps = ['-o';'-x';'-*';'-v'];
dottps = ['o';'x';'*';'v'];
for i = 1:3
    ctemp = nonzeros(counter(1,:,i)-1)
    l = length(ctemp)
    
    effAve(1,1:l,i) = nonzeros(effAve(1,:,i))./ctemp;
    effAve(1,l+1:vl,i) = 0;
    aveE(1,1:l,i) = nonzeros(aveE(1,:,i))./ctemp;
    aveE(1,l+1:vl,i) = 0;
    aveENP(1,1:l,i) = nonzeros(aveENP(1,:,i))./ctemp;
    aveENP(1,l+1:vl,i) = 0;
    avBD(1,1:l,i) = nonzeros(avBD(1,:,i))./ctemp;
    avBD(1,l+1:vl,i) = 0;
    vfAve(1,1:l,i)  = nonzeros(vfAve(1,:,i))./ctemp;
    vfAve(1,l+1:vl,i) = 0;
    aveMaxED(1,1:l,i) = nonzeros(aveMaxED(1,:,i))./ctemp;
    aveMaxED(1,l+1:vl,i) = 0;
    for num = 1:length(vfs)
        stdev(num,i) = std(nonzeros(effchart(num,1:nonzeros(counter(1,num,i)-1),i)));
        stdevE(num,i) = std(nonzeros(avEchart(num,1:nonzeros(counter(1,num,i)-1),i)));
        stdevENP(num,i) = std(nonzeros(avENPchart(num,1:nonzeros(counter(1,num,i)-1),i)));
        stdevBD(num,i) = std(nonzeros(bdchart(num,1:nonzeros(counter(1,num,i)-1),i)))
        stdevMED(num,i) = std(nonzeros(maxEdchart(num,1:nonzeros(counter(1,num,i)-1),i)));
        
        stdev(num,4) = std(nonzeros(effchart(num,1:nonzeros(counter(1,num,4)-1),4)));
        stdevE(num,4) = std(nonzeros(avEchart(num,1:nonzeros(counter(1,num,4)-1),4)));
        stdevENP(num,4) = std(nonzeros(avENPchart(num,1:nonzeros(counter(1,num,4)-1),4)));
        stdevBD(num,4) = std(nonzeros(bdchart(num,1:nonzeros(counter(1,num,4)-1),4)))
        stdevMED(num,4) = std(nonzeros(maxEdchart(num,1:nonzeros(counter(1,num,4)-1),4)));
        
%         stdevS(num,i) = std(nonzeros(effchartS(num,1:nonzeros(counterS(1,num,i)-1),i)));
%         stdevES(num,i) = std(nonzeros(avEchartS(num,1:nonzeros(counterS(1,num,i)-1),i)));
%         stdevENPS(num,i) = std(nonzeros(avENPchartS(num,1:nonzeros(counterS(1,num,i)-1),i)));
%         stdevBDS(num,i) = std(nonzeros(bdchartS(num,1:nonzeros(counterS(1,num,i)-1),i)))
%         stdevMEDS(num,i) = std(nonzeros(maxEdchartS(num,1:nonzeros(counterS(1,num,i)-1),i)));
%         
%         stdevS(num,4) = std(nonzeros(effchartS(num,1:nonzeros(counterS(1,num,4)-1),4)));
%         stdevES(num,4) = std(nonzeros(avEchartS(num,1:nonzeros(counterS(1,num,4)-1),4)));
%         stdevENPS(num,4) = std(nonzeros(avENPchartS(num,1:nonzeros(counterS(1,num,4)-1),4)));
%         stdevBDS(num,4) = std(nonzeros(bdchartS(num,1:nonzeros(counterS(1,num,4)-1),4)))
%         stdevMEDS(num,4) = std(nonzeros(maxEdchartS(num,1:nonzeros(counterS(1,num,4)-1),4)));
        

    end
    
%     figure(perm);
%     errorbar(nonzeros(vfAve(:,:,i)), nonzeros(effAve(1,:,i))', 3*stdev(counter(1,:,i)-1~=0,i), dottps(i), 'linewidth',1.5);
%     
%     figure(efld)
%     errorbar(nonzeros(vfAve(:,:,i)), nonzeros(aveE(1,:,i))', 3*stdevE(counter(1,:,i)-1~=0,i), linetps(i), 'linewidth',1.5);
% 
%     figure(efld2)
%     errorbar(nonzeros(vfAve(:,:,i)), nonzeros(aveE(1,:,i))', 3*stdevE(counter(1,:,i)-1~=0,i), dottps(i), 'linewidth',1.5);
%     
%     figure(efldBT)
%     errorbar(nonzeros(vfAve(:,:,i)), nonzeros(aveENP(1,:,i))', 3*stdevENP(counter(1,:,i)-1~=0,i), dottps(i), 'linewidth',1.5);
%     
%     figure(eden);
%     errorbar(nonzeros(vfAve(:,:,i)), nonzeros(aveMaxED(1,:,i))', 3*stdevMED(counter(1,:,i)-1~=0,i), dottps(i),'linewidth',1.5);
%     
%     figure(bdfs);
%     errorbar(nonzeros(vfAve(:,:,i)), nonzeros(avBD(1,:,i))', stdevBD(counter(1,:,i)-1~=0,i), dottps(i),'linewidth',1.5);
%     
end
ctemp = nonzeros(counter(1,:,4)-1);
l = length(ctemp)

ctemp = nonzeros(counter(1,:,4)-1);
effAve(1,1:l,4) = nonzeros(sum(sum(effchart(:,:,1:3),3),2))./ctemp;
effAve(1,l+1:vl,4) = 0;
aveE(1,1:l,4) = nonzeros(sum(sum(avEchart(:,:,1:3),2),3))./ctemp;
aveE(1,l+1:vl,4) = 0;
aveENP(1,1:l,4) = nonzeros(sum(sum(avENPchart(:,:,1:3),2),3))./ctemp;
aveENP(1,l+1:vl,4) = 0;
avBD(1,1:l,4) = nonzeros(sum(sum(bdchart(:,:,1:3),2),3))./ctemp;
avBD(1,l+1:vl,4) = 0;
aveMaxED(1,1:l,4) = nonzeros(sum(sum(maxEdchart(:,:,1:3),2),3))./ctemp;
aveMaxED(1,l+1:vl,4) = 0;

vfAve(1,1:l,4)  = nonzeros(vfAve(1,:,4))./ctemp;
vfAve(1,l+1:vl,4) = 0;

counter(1,:,4) = sum(counter(1,:,:)-1,3);

%------------------for loop ended here ------------%
%% EFFECTIVE PERMITTIVITY PLOTS
figure(perm)
c = get(gca,'colororder') % variation of permittivity as scatterplot
scatter(btVFsct(:,1),effSct(:,1),[],c(3,:),'.');
scatter(btVFsct(:,2),effSct(:,2),[],c(3,:),'o');
scatter(btVFsct(:,3),effSct(:,3),[],c(3,:),'x');
errorbar(nonzeros(vfAve(:,:,4)), nonzeros(effAve(1,:,4))', 3*stdev(counter(1,:,4)~=0,4), '--k', 'linewidth',1);
ylabel("Effective Relative Permittivity [F/M]");
xlabel("Volume Fraction of Nanoparticles [%]");
legend('Percolating Geometries','Near Percolating Geometries', 'Non-Percolating Geometries', 'All Geometries', 'Location', 'southeast')

figure(perm2)% mean of permittivity with perry and bruggeman
plot(nonzeros(vfAve(:,:,4)), nonzeros(effAve(1,:,4))', '--k', 'linewidth',1.5);
plot(perrydat(1,:), perrydat(2,:),'ok','linewidth',1.5)
%plot(vfs, brug,'--k','linewidth',1.5)
ylabel("Effective Relative Permittivity [F/M]");
xlabel("Volume Fraction of Nanoparticles [%]");
legend('Mean of tested geometries','Experimental Data (Perry 2009)', 'Bruggeman Model', 'All Geometries', 'Location', 'southeast')



%% ENERGY DENSITY  PLOTS

figure(eden)
errorbar(nonzeros(vfAve(:,:,4)), nonzeros(aveMaxED(1,:,4))', 3*stdevMED(counter(1,:,4)~=0,4), dottps(4), 'linewidth',2);

ylabel("Maximum Energy Density [J/cm^3]");
xlabel("Volume Fraction of Nanoparticles [%]");
legend('Percolating Geometries','Near Percolating Geometries', 'Non-Percolating Geometries', 'All Geometries', 'Location', 'southeast')

figure(eden2)
%plot(nonzeros(vfAve(:,:,4)), nonzeros(aveMaxED(1,:,4))', linetps(4), 'linewidth',1.5);
c = get(gca,'colororder')
scatter(nonzeros(btVFsct(:,1)),nonzeros(maxEdSct(:,1)),[],c(3,:),'.');
scatter(nonzeros(btVFsct(:,2)),nonzeros(maxEdSct(:,2)),[],c(3,:),'o');
scatter(nonzeros(btVFsct(:,3)),nonzeros(maxEdSct(:,3)),[],c(3,:),'x');
plot(vfAve(:,:,4),aveMaxED(:,:,4),'--k','linewidth',1.5) % Plot mean line
%plot(perryed(1,:), perryed(2,:),'.-.','linewidth',1.5)
ylabel("Maximum Energy Density [J/cm^3]");
xlabel("Volume Fraction of Nanoparticles [%]");
legend('Percolating Geometries','Near Percolating Geometries', 'Non-Percolating Geometries', 'Mean', 'Location', 'southeast')
 
%% BREAKDOWN FIELD STRENGTH PLOTS
figure(bdfs)
c = get(gca,'colororder') % variation of permittivity as scatterplot
scatter(nonzeros(btVFsct(:,3)),nonzeros(bdSct(:,3)).*1E-3,[],c(3,:),'x');
scatter(nonzeros(btVFsct(:,2)),nonzeros(bdSct(:,2)).*1E-3,[],c(3,:),'o');
scatter(nonzeros(btVFsct(:,1)),nonzeros(bdSct(:,1)).*1E-3,[],c(3,:),'.');

errorbar(nonzeros(vfAve(:,:,4)), nonzeros(avBD(1,:,4))'.*1E-3, 3*stdevBD(counter(1,:,4)~=0,4).*1E-3, '-k', 'linewidth',1);
%plot(perrybd(1,:), perrybd(2,:),'ok','linewidth',1.5)
ylabel("Breakdown Field Strength [MV/M]");
xlabel("Volume Fraction of Nanoparticles [%]");
legend('Percolating Geometries','Near Percolating Geometries', 'Non-Percolating Geometries', 'All Geometries','Measured Data (Perry)', 'Location', 'northeast')

figure(bdfs2)
plot(nonzeros(vfAve(:,:,4)), nonzeros(avBD(1,:,4))', '--k', 'linewidth',1.5);
plot(perrybd(1,:), perrybd(2,:),'ok','linewidth',1.5)
ylabel("Breakdown Field Strength [MV/M]");
xlabel("Volume Fraction of Nanoparticles[%]");
legend('Simulated Geometries','Measured Data (Perry)', 'Non-Percolating Geometries', 'All Geometries', 'Location', 'northeast')


%% MEAN E-FIELD PLOTS
figure(efld)
scatter(btVFsct(:,1),avESct(:,1)/Em2p,[],c(3,:),'.');
scatter(btVFsct(:,2),avESct(:,2)/Em2p,[],c(3,:),'o');
scatter(btVFsct(:,3),avESct(:,3)/Em2p,[],c(3,:),'x');
scatter(btVFsct(:,1),avENPSct(:,1)/Em2p,[],c(3,:),'.');
scatter(btVFsct(:,2),avENPSct(:,2)/Em2p,[],c(3,:),'o');
scatter(btVFsct(:,3),avENPSct(:,3)/Em2p,[],c(3,:),'x');
errorbar(nonzeros(vfAve(:,:,4)), nonzeros(aveE(1,:,4)/Em2p)', 3*stdevE(counter(1,:,4)~=0,4)/Em2p, '--k', 'linewidth',1);
errorbar(nonzeros(vfAve(:,:,4)), nonzeros(aveENP(1,:,4)/Em2p)', 3*stdevENP(counter(1,:,4)~=0,4)/Em2p, '--k', 'linewidth',1);
ylabel("Average Normalized Electric Field");
xlabel("Volume Fraction of Nanoparticles [%]");
%legend('Inside the Host Material','Inside the Nanoparticles', 'Location', 'southeast')

% figure(efldBT)
% scatter(btVFsct(:,1),avENPSct(:,1),[],c(1,:),'o');
% scatter(btVFsct(:,2),avENPSct(:,2),[],c(2,:),'o');
% scatter(btVFsct(:,3),avENPSct(:,3),[],c(3,:),'o');
% ylabel("Average Electric Field within Nanoparticles");
% xlabel("Volume Fraction of BT Nanoparticles");
% errorbar(nonzeros(vfAve(:,:,4)), nonzeros(aveENP(1,:,4))', stdevENP(counter(1,:,4)~=0,4), linetps(4), 'linewidth',2);

% legend('Percolating Geometries','Near Percolating Geometries', 'Non-Percolating Geometries', 'All Geometries', 'Location', 'southeast')
%% OTHER MISC. PLOTS
% max V vs. permittivity
figure(volts)
%errorbar(nonzeros(vfAve(:,:,4)), nonzeros(avBD(1,:,4))', 3*stdevBD(counter(1,:,4)~=0,4), linetps(4), 'linewidth',1);
c = get(gca,'colororder') % variation of permittivity as scatterplot
scatter(nonzeros(btVFsct(:,1)),nonzeros(MaxVolt2p(:,1)),[],c(3,:),'o');
scatter(nonzeros(btVFsct(:,2)),nonzeros(MaxVolt2p(:,2)),[],c(3,:),'o');
scatter(nonzeros(btVFsct(:,3)),nonzeros(MaxVolt2p(:,3)),[],c(3,:),'o');

scatter(nonzeros(btVFsct(:,1)),nonzeros(maxVolt2pNP(:,1)),[],c(3,:),'.');
scatter(nonzeros(btVFsct(:,2)),nonzeros(maxVolt2pNP(:,2)),[],c(3,:),'.');
scatter(nonzeros(btVFsct(:,3)),nonzeros(maxVolt2pNP(:,3)),[],c(3,:),'.');
ylabel("Maximum Voltage [V]");
xlabel("Volume Fraction of Nanoparticles [%]");
legend('Percolating in Host','Near Percolating in Host', 'Non-Percolating in Host', 'Percolating in NP','Near Percolating in NP', 'Non-Percolating in NP', 'Location', 'southeast')

figure(volts2)
%errorbar(nonzeros(vfAve(:,:,4)), nonzeros(avBD(1,:,4))', 3*stdevBD(counter(1,:,4)~=0,4), linetps(4), 'linewidth',1);
c = get(gca,'colororder') % variation of permittivity as scatterplot
scatter(nonzeros(effSct(:,1)),nonzeros(MaxVolt2p(:,1)),[],c(3,:),'o');
scatter(nonzeros(effSct(:,2)),nonzeros(MaxVolt2p(:,2)),[],c(3,:),'o');
scatter(nonzeros(effSct(:,3)),nonzeros(MaxVolt2p(:,3)),[],c(3,:),'o');

scatter(nonzeros(effSct(:,1)),nonzeros(maxVolt2pNP(:,1)),[],c(3,:),'.');
scatter(nonzeros(effSct(:,2)),nonzeros(maxVolt2pNP(:,2)),[],c(3,:),'.');
scatter(nonzeros(effSct(:,3)),nonzeros(maxVolt2pNP(:,3)),[],c(3,:),'.');
ylabel("Maximum Voltage [V]");
xlabel("Effective Permittivity [F/M]");
legend('Percolating in Host','Near Percolating in Host', 'Non-Percolating in Host','Percolating in NP','Near Percolating in NP', 'Non-Percolating in NP', 'Location', 'northeast')


%% COUNTER BAR PLOT
figure(countbars)
allcounts = [counter(:,:,1)'-1, counter(:,:,2)'-1, counter(:,:,3)'-1];% counter(:,:,4)'];
%hb = bar(vfAve(:,:,4)', allcounts, 'stacked');
ylabel("Number of Microstructures");
xlabel("Volume Fraction of Nanoparticles [%]");
ylim([0,80]);
legend('Percolaging Geometries', 'Near Percolating Geometries', 'Non-Percolating Geometries')
saveas(countbars,'bars2','fig');
saveas(countbars,'bars2','jpeg');
saveas(perm,'permValid2','fig');
saveas(perm,'permValid2','jpeg');
saveas(eden,'eden2','fig');
saveas(eden,'eden2','jpeg');
saveas(efldBT,'efldBT2','fig');
saveas(efldBT,'efldBT2','jpeg');
saveas(bdfs,'bdfs2','fig');
saveas(bdfs,'bdfs2','jpeg');
saveas(efld,'efld_hst2','fig');
saveas(efld,'efld_hst2','jpeg');
saveas(perm2,'perm3','fig');
saveas(perm2,'perm3','jpeg');
saveas(eden2,'eden3','fig');
saveas(eden2,'eden3','jpeg');
saveas(efldBT2,'efldBT3','fig');
saveas(efldBT2,'efldBT3','jpeg');
saveas(bdfs2,'bdfs3','fig');
saveas(bdfs2,'bdfs3','jpeg');
saveas(efld2,'efld_hst3','fig');
saveas(efld2,'efld_hst3','jpeg');
saveas(volts,'volts1','fig');
saveas(volts,'volts1','jpeg');
saveas(volts2,'volts2','fig');
saveas(volts2,'volts2','jpeg');

figure(MCP)
hold on
%MCPdat = [3.9,3.9,6.3,6.7,7.1,7.1,8.5,8.5,21.0,21.0,21.0,19.0,19.0,23.4,26,26,26.8,29,31,60,95,95,95,95,50,98;1400,1310,630,870,760,920,620,720,470,400,300,250,600,350,320,420,420,380,400,320,300,250,200,100,330]
Mcpdat = [3.9,1400;3.9,1310;6.3,630;6.7,870;7.1,760;7.1,920;8.5,620;8.5,720;21.0,470;21.0,400;21.0,300;19.0,250;19.0,600;23.4,350;26,320;26,420;26.8,420;29,380;31,400;60,320;95,300;95,250;95,200;95,100;50,300;98,330;]
scatter(effSct(:,3),bdSct(:,3).^2,[],c(3,:),'x');
scatter(effSct(:,2),bdSct(:,2).^2,[],c(3,:),'o');
scatter(effSct(:,1),bdSct(:,1).^2,[],c(3,:),'.');
c = get(gca,'colororder') % variation of permittivity as scatterplot

% plot(Mcpdat(:,1), Mcpdat(:,2).^2, 'ok', 'linewidth',1);
ylabel("Breakdown Field Strength^2 [MV/M]");
xlabel("Effective Relative Permittivity [F/M]");
legend('Percolating Geometries','Near Percolating Geometries', 'Non-Percolating Geometries', 'Data from ', 'Location', 'southeast')
set(gca,'xscale','log')
set(gca,'yscale','log')

% Get a list of all variables
allvars = whos;
% Identify the variables that ARE NOT graphics handles. This uses a regular
% expression on the class of each variable to check if it's a graphics object
tosave = cellfun(@isempty, regexp({allvars.class}, '^matlab\.(ui|graphics)\.'));
% Pass these variable names to save
save('dataset2_nf.mat', allvars(tosave).name)


%% Unpack & Plot Large Datasets
% clc
% clear all
PSeeds = [3,31,40,47];
NPseeds = [12,23,24,37,53,59,62,67,68,71,77,79];
therest = [13,14,18,41,48,55,66,73];
BTbd = 1e6;
dirnm = 'PercRuns/Run_1033/';

perrydat = [0,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6; 11.9, 12.5, 15.2, 16.4, 18.6, 22.3, 26.2, 35.1, 37.5]
perrybd = [0,0.05,0.1,0.15,0.2,0.3,0.4,0.5;380,344,345,300,235,230,218,213];
perryed = [0,0.1,0.2,0.3,0.4,0.5;1.6E6,1.6E6,2E6,2.3E6,2.35E6,3.2E6];
% perm = figure;
% hold on
% efld = figure;
% hold on
% eden = figure;
% hold on
% bdfs = figure;
% hold on
% efldBT = figure;
% hold on
% countbars = figure;
% hold on
% perm2 = figure;
% hold on
% efld2 = figure;
% hold on
% eden2 = figure;
% hold on
% bdfs2 = figure;
% hold on
% efldBT2 = figure;
% hold on
% volts = figure;
% hold on
% volts2 = figure;
% hold on


Em2p = (1)/(65*2e-9)
vfs = [0.05:0.05:0.6]

vl = length(vfs)
files = dir(dirnm);
effchart = zeros(vl,100,4);
avEchart = effchart;
avENPchart = effchart;
bdchart = effchart;
brugchart = effchart;
maxEdchart = effchart;
[numF,~] = size(files);
numF = numF -2;
files = files(3:numF+2);
vfAve = zeros(1,length(vfs),4);
effAve = vfAve;
aveMaxED = vfAve;
counter = vfAve+1;
countvf = vfAve;
aveE = vfAve;
aveENP = vfAve;
avBD = vfAve;


clear VFrac
clear eff
clear btVFr
clear bdv2pBT
clear BTstdE
clear MaxVolt2p
clear maxVolt2pNP
clear MaxEd2pNP
inv = 0;
alltheVF = zeros(1,numF);
for j = 2:numF
    Fn = j - inv ;
    fnm = files(Fn).name;
    fnm = [dirnm,fnm];
    importfile(fnm)
    if length(btVF) < 2
        alltheVF(Fn-1) = btVF;
    else
        files(Fn) = [];
        inv = inv + 1;
    end
end
alltheVF = alltheVF(1:(numF-inv-1))
[~,VFsort]=sort(alltheVF); %Get the order of B
files = files(VFsort);

brg = zeros(1,length(vfs))
btVFr = zeros(numF,3);
btVFsct = btVFr;
effSct = btVFr;
avESct = btVFr;
avENPSct = btVFr;
bdSct = btVFr;
maxEdSct = btVFr;
effSctS = btVFr;
avESctS = btVFr;
avENPSctS = btVFr;
bdSctS = btVFr;
maxEdSctS = btVFr;
              
              
for j = 2:numF-1
    i = j;
    fnm = files(i).name;
    fnm = [dirnm,fnm];
    if length(fnm) > 30
        importfile(fnm)
        
        if(btVF > 0)
            
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
            if p~=4
                % Scale Things...
                eff2p = eff2p*2;
%                 Emean2p = Emean2p;
%                 Emean2pNP = Emean2pNP;
                PVDFbd = PVDFbd*1e-6;
                BTbd = BTbd*1e-6;
                % ^^ fix later in Tester3 so that variation (stdev) is 
                % correctly calculated
                
                
                eff(j,1,p) = eff2p;
                eff(j,2,p) = eff3p;
                eff(j,3,p) = eff3pH;
                eff(j,4,p) = eff3pS;
                % use path to compute average electric field along breakdown path
                stE2pNP = std(std(std(E2p(OmegaBT(:,:,1:gd-1)==2)))); %.*(Delta(1)*Delta(2)*Delta(3));
                bdv2pBT(j,p) = (gd*BTbd)./((Emean2pNP+stE2p));
                maxVolt2pNP(j,p) =   BTbd ./ ( Emean2pNP + stE2pNP );
                MaxEd2pNP(j,p) = Edens2p * maxVolt2pNP(j,p).^2;
                MaxVolt2p(j,p) =  PVDFbd ./ ( Emean2p + stE2p );
                avE(j,1,p) = avE2p;
                
                Ed = (E2p.^2 .*(1/2)*e0 .* Mats(OmegaBT(:,:,1:Zsz-1)+1)); % Delta(1)*Delta(2)*Delta(3)
                Ed = sum(Ed,3);
                Ed = sum(Ed,2);
                Ed2p = sum(Ed,1);
                Edens2p(vfi) = Ed2p/(gd*gd*gd);
                
                % mcE(vfs == volfrac) = mcE
                % maxEd(j) = MaxEd2p;
                if ((maxVolt2pNP(j,p) < MaxVolt2p(j,p)) && p~=3)
                    maxEd(j,p) = (maxVolt2pNP(j,p)*1E6)^2*Edens2p;
                    % FIX THIS --->> maxEdS(j,p) = (maxVolt3pNPS(j,p)*1E6)^2*Edens3pS;
                    bd = maxVolt2pNP(j,p)/(Delta(3)*gd);%BTbd/((Emean2pNP+stE2pNP)/(3.1E6)); % Approximate y intercept
%                     bdS = maxVolt2pNPS(j,p)/(Delta(3)*gd);%BTbd/((Emean2pNP+stE2pNP)/(3.1E6)); % Approximate y intercept
                     % Account for breakdown voltage inside nanoparticle
                     % HERE
                else
                    maxEd(j,p) = (MaxVolt2p(j,p)*1E6)^2*Edens2p;
%                     maxEdS(j,p) = (maxVolt3pS(j,p)*1E6)^2*Edens3pS;
                    bd = MaxVolt2p(j,p)/(Delta(3)*gd);%PVDFbd/((Emean2p+stE2p)/(7.8E6)); % divide by approximate y intercept...
%                     bdS = MaxVolt3pS(j,p)/(Delta(3)*gd);%PVDFbd/((Emean2p+stE2p)/(7.8E6)); % divide by approximate y intercept...

                    % AND HERE
                end
                volfracs = round(vfs*1000) == round(volfrac*1000);
                disp(volfracs)
                eff2p = eff2p/2; % Adjust permittivity
                aveE(1,volfracs,p) = aveE(1,volfracs,p) + Emean2p;
                aveENP(1,volfracs,p) = aveENP(1,volfracs,p) + Emean2pNP;
                counter(1,volfracs,p) = counter(1,volfracs,p) + 1;
                effAve(1,volfracs,p) = effAve(1,volfracs,p) + eff2p;
                avBD(1,volfracs,p) = avBD(1,volfracs,p) + bd;
                aveMaxED(1,volfracs,p) = aveMaxED(1,volfracs,p) + maxEd(j,p);
                effchart(volfracs, counter(1,volfracs,p),p) = eff2p;
                avEchart(volfracs, counter(1,volfracs,p),p) = Emean2p;
                avENPchart(volfracs, counter(1,volfracs,p),p) = Emean2pNP;
                bdchart(volfracs, counter(1,volfracs,p),p) = bd;
                maxEdchart(volfracs, counter(1,volfracs,p),p) = maxEd(j,p);
%                 
%                 eff3pS = eff3pS/2; % Adjust permittivity
%                 aveES(1,volfracs,p) = aveES(1,volfracs,p) + Emean3pS;
%                 aveENPS(1,volfracs,p) = aveENPS(1,volfracs,p) + Emean3pNPS;
%                 counterS(1,volfracs,p) = counterS(1,volfracs,p) + 1;
%                 effAveS(1,volfracs,p) = effAveS(1,volfracs,p) + eff3pS;
%                 avBDS(1,volfracs,p) = avBDS(1,volfracs,p) + bd;
%                 aveMaxEDS(1,volfracs,p) = aveMaxEDS(1,volfracs,p) + maxEdS(j,p);
%                 effchartS(volfracs, counterS(1,volfracs,p),p) = eff3pS;
%                 avEchartS(volfracs, counterS(1,volfracs,p),p) = Emean3pS;
%                 avENPchartS(volfracs, counterS(1,volfracs,p),p) = Emean3pNPS;
%                 bdchartS(volfracs, counterS(1,volfracs,p),p) = bdS;
%                 maxEdchartS(volfracs, counterS(1,volfracs,p),p) = maxEdS(j,p);
                
                vfchart(volfracs, counter(1,volfracs,p),p) = btVF;
                
                cSct = sum(counter-1,2)
                btVFr(   cSct(p),p) = btVF;
                btVFsct( cSct(p),p) = btVF;
                effSct(  cSct(p),p) = eff2p;
                avESct(  cSct(p),p) = Emean2p;
                avENPSct(cSct(p),p) = Emean2pNP;
                bdSct(   cSct(p),p) = bd;
                maxEdSct(j,p) = maxEd(j,p);
                                
%                 cSctS = sum(counterS-1,2)
%                 btVFrS(   cSct(p),p) = btVF;
%                 btVFsctS( cSct(p),p) = btVF;
%                 effSctS(  cSctS(p),p) = eff3pS;
%                 avESctS(  cSctS(p),p) = Emean3pS;
%                 avENPSctS(cSctS(p),p) = Emean3pNPS;
%                 bdSctS(   cSctS(p),p) = bdS;
%                 maxEdSctS(j,p) = maxEdS(j,p);
                
                
                vfAve(1,volfracs,p) = vfAve(1,volfracs,p) + btVF;
                vfAve(1,volfracs,4) = vfAve(1,volfracs,4) + btVF;
%                 vfAveS(1,volfracs,p) = vfAveS(1,volfracs,p) + btVF;
%                 vfAveS(1,volfracs,4) = vfAveS(1,volfracs,4) + btVF;
                countvf(1,volfracs,p) = countvf(1,volfracs,p) + 1;
                
                countvf(1,volfracs,4) = countvf(1,volfracs,4) + 1;
                counter(1,volfracs,4) = counter(1,volfracs,4) + 1;
%                 effchartS(volfracs, counter(1,volfracs,4),4) = eff2p;
%                 avEchartS(volfracs, counter(1,volfracs,4),4) = Emean2p;
%                 avENPchartS(volfracs, counter(1,volfracs,4),4) = Emean2pNP;
%                 bdchartS(volfracs, counter(1,volfracs,4),4) = bd;
%                 maxEdchartS(volfracs, counter(1,volfracs,4),4) = maxEd(j,p);
                
                                syms x % effective permittivity of solid
            if((phases==2) && (brg(mci,volfracs)==0))
                brgEq = (1-volfrac)*((CxPvty(1)-x)/(CxPvty(1)+x)) + (volfrac)*((CxPvty(3)-x)/(CxPvty(3)+x)) == 0;
            elseif((phases == 3)  && (brg(mci,volfracs)==0))
                brgEq = (1-volfrac)*((CxPvty(1)-x)/(CxPvty(1)+x)) + (volfrac)*((CxPvty(2)-x)/(CxPvty(2)+x)) + (volfrac2)*((CxPvty(3)-x)/(CxPvty(3)+x))==0;
            end
            soln = solve(brgEq,x);
            brg(mci,volfracs) = max(soln);


            end
        end
    end
end

for bcnt = 1:length(volfracs)
    brgEq = (1-vfs(bcnt))*((Mats(1)-x)/(Mats(1)+x)) + (vfs(bcnt))*((Mats(3)-x)/(Mats(3)+x)) == 0;
    soln = solve(brgEq,x);
    brugg = max(soln);
    brug(bcnt) = brugg;
end
%     figure(perm);
%     plot(vfs, brug, '--', 'linewidth',2);
    
% line identifiers
linetps = ['-o';'-x';'-*';'-v'];
dottps = ['o';'x';'*';'v'];
for i = 1:3
    ctemp = nonzeros(counter(1,:,i)-1)
    l = length(ctemp)
    
    effAve(1,1:l,i) = nonzeros(effAve(1,:,i))./ctemp;
    effAve(1,l+1:vl,i) = 0;
    aveE(1,1:l,i) = nonzeros(aveE(1,:,i))./ctemp;
    aveE(1,l+1:vl,i) = 0;
    aveENP(1,1:l,i) = nonzeros(aveENP(1,:,i))./ctemp;
    aveENP(1,l+1:vl,i) = 0;
    avBD(1,1:l,i) = nonzeros(avBD(1,:,i))./ctemp;
    avBD(1,l+1:vl,i) = 0;
    vfAve(1,1:l,i)  = nonzeros(vfAve(1,:,i))./ctemp;
    vfAve(1,l+1:vl,i) = 0;
    aveMaxED(1,1:l,i) = nonzeros(aveMaxED(1,:,i))./ctemp;
    aveMaxED(1,l+1:vl,i) = 0;
    for num = 1:length(vfs)
        stdev(num,i) = std(nonzeros(effchart(num,1:nonzeros(counter(1,num,i)-1),i)));
        stdevE(num,i) = std(nonzeros(avEchart(num,1:nonzeros(counter(1,num,i)-1),i)));
        stdevENP(num,i) = std(nonzeros(avENPchart(num,1:nonzeros(counter(1,num,i)-1),i)));
        stdevBD(num,i) = std(nonzeros(bdchart(num,1:nonzeros(counter(1,num,i)-1),i)))
        stdevMED(num,i) = std(nonzeros(maxEdchart(num,1:nonzeros(counter(1,num,i)-1),i)));
        
        stdev(num,4) = std(nonzeros(effchart(num,1:nonzeros(counter(1,num,4)-1),4)));
        stdevE(num,4) = std(nonzeros(avEchart(num,1:nonzeros(counter(1,num,4)-1),4)));
        stdevENP(num,4) = std(nonzeros(avENPchart(num,1:nonzeros(counter(1,num,4)-1),4)));
        stdevBD(num,4) = std(nonzeros(bdchart(num,1:nonzeros(counter(1,num,4)-1),4)))
        stdevMED(num,4) = std(nonzeros(maxEdchart(num,1:nonzeros(counter(1,num,4)-1),4)));
        
%         stdevS(num,i) = std(nonzeros(effchartS(num,1:nonzeros(counterS(1,num,i)-1),i)));
%         stdevES(num,i) = std(nonzeros(avEchartS(num,1:nonzeros(counterS(1,num,i)-1),i)));
%         stdevENPS(num,i) = std(nonzeros(avENPchartS(num,1:nonzeros(counterS(1,num,i)-1),i)));
%         stdevBDS(num,i) = std(nonzeros(bdchartS(num,1:nonzeros(counterS(1,num,i)-1),i)))
%         stdevMEDS(num,i) = std(nonzeros(maxEdchartS(num,1:nonzeros(counterS(1,num,i)-1),i)));
%         
%         stdevS(num,4) = std(nonzeros(effchartS(num,1:nonzeros(counterS(1,num,4)-1),4)));
%         stdevES(num,4) = std(nonzeros(avEchartS(num,1:nonzeros(counterS(1,num,4)-1),4)));
%         stdevENPS(num,4) = std(nonzeros(avENPchartS(num,1:nonzeros(counterS(1,num,4)-1),4)));
%         stdevBDS(num,4) = std(nonzeros(bdchartS(num,1:nonzeros(counterS(1,num,4)-1),4)))
%         stdevMEDS(num,4) = std(nonzeros(maxEdchartS(num,1:nonzeros(counterS(1,num,4)-1),4)));
        

    end
    
%     figure(perm);
%     errorbar(nonzeros(vfAve(:,:,i)), nonzeros(effAve(1,:,i))', 3*stdev(counter(1,:,i)-1~=0,i), dottps(i), 'linewidth',1.5);
%     
%     figure(efld)
%     errorbar(nonzeros(vfAve(:,:,i)), nonzeros(aveE(1,:,i))', 3*stdevE(counter(1,:,i)-1~=0,i), linetps(i), 'linewidth',1.5);
% 
%     figure(efld2)
%     errorbar(nonzeros(vfAve(:,:,i)), nonzeros(aveE(1,:,i))', 3*stdevE(counter(1,:,i)-1~=0,i), dottps(i), 'linewidth',1.5);
%     
%     figure(efldBT)
%     errorbar(nonzeros(vfAve(:,:,i)), nonzeros(aveENP(1,:,i))', 3*stdevENP(counter(1,:,i)-1~=0,i), dottps(i), 'linewidth',1.5);
%     
%     figure(eden);
%     errorbar(nonzeros(vfAve(:,:,i)), nonzeros(aveMaxED(1,:,i))', 3*stdevMED(counter(1,:,i)-1~=0,i), dottps(i),'linewidth',1.5);
%     
%     figure(bdfs);
%     errorbar(nonzeros(vfAve(:,:,i)), nonzeros(avBD(1,:,i))', stdevBD(counter(1,:,i)-1~=0,i), dottps(i),'linewidth',1.5);
%     
end
ctemp = nonzeros(counter(1,:,4)-1);
l = length(ctemp)

ctemp = nonzeros(counter(1,:,4)-1);
effAve(1,1:l,4) = nonzeros(sum(sum(effchart(:,:,1:3),3),2))./ctemp;
effAve(1,l+1:vl,4) = 0;
aveE(1,1:l,4) = nonzeros(sum(sum(avEchart(:,:,1:3),2),3))./ctemp;
aveE(1,l+1:vl,4) = 0;
aveENP(1,1:l,4) = nonzeros(sum(sum(avENPchart(:,:,1:3),2),3))./ctemp;
aveENP(1,l+1:vl,4) = 0;
avBD(1,1:l,4) = nonzeros(sum(sum(bdchart(:,:,1:3),2),3))./ctemp;
avBD(1,l+1:vl,4) = 0;
aveMaxED(1,1:l,4) = nonzeros(sum(sum(maxEdchart(:,:,1:3),2),3))./ctemp;
aveMaxED(1,l+1:vl,4) = 0;

vfAve(1,1:l,4)  = nonzeros(vfAve(1,:,4))./ctemp;
vfAve(1,l+1:vl,4) = 0;

counter(1,:,4) = sum(counter(1,:,:)-1,3);

%------------------for loop ended here ------------%
%% EFFECTIVE PERMITTIVITY PLOTS
figure(perm)
c = get(gca,'colororder') % variation of permittivity as scatterplot
scatter(btVFsct(:,1),effSct(:,1),[],c(4,:),'.');
scatter(btVFsct(:,2),effSct(:,2),[],c(4,:),'o');
scatter(btVFsct(:,3),effSct(:,3),[],c(4,:),'x');
errorbar(nonzeros(vfAve(:,:,4)), nonzeros(effAve(1,:,4))', 3*stdev(counter(1,:,4)~=0,4), '--k', 'linewidth',1);
ylabel("Effective Relative Permittivity [F/M]");
xlabel("Volume Fraction of Nanoparticles [%]");
legend('Percolating Geometries','Near Percolating Geometries', 'Non-Percolating Geometries', 'All Geometries', 'Location', 'southeast')

figure(perm2)% mean of permittivity with perry and bruggeman
plot(nonzeros(vfAve(:,:,4)), nonzeros(effAve(1,:,4))', '--k', 'linewidth',1.5);
plot(perrydat(1,:), perrydat(2,:),'ok','linewidth',1.5)
%plot(vfs, brug,'--k','linewidth',1.5)
ylabel("Effective Relative Permittivity [F/M]");
xlabel("Volume Fraction of Nanoparticles [%]");
legend('Mean of tested geometries','Experimental Data (Perry 2009)', 'Bruggeman Model', 'All Geometries', 'Location', 'southeast')



%% ENERGY DENSITY  PLOTS

figure(eden)
errorbar(nonzeros(vfAve(:,:,4)), nonzeros(aveMaxED(1,:,4))', 3*stdevMED(counter(1,:,4)~=0,4), dottps(4), 'linewidth',2);

ylabel("Maximum Energy Density [J/cm^3]");
xlabel("Volume Fraction of Nanoparticles [%]");
legend('Percolating Geometries','Near Percolating Geometries', 'Non-Percolating Geometries', 'All Geometries', 'Location', 'southeast')

figure(eden2)
%plot(nonzeros(vfAve(:,:,4)), nonzeros(aveMaxED(1,:,4))', linetps(4), 'linewidth',1.5);
c = get(gca,'colororder')
scatter(nonzeros(btVFsct(:,1)),nonzeros(maxEdSct(:,1)),[],c(4,:),'.');
scatter(nonzeros(btVFsct(:,2)),nonzeros(maxEdSct(:,2)),[],c(4,:),'o');
scatter(nonzeros(btVFsct(:,3)),nonzeros(maxEdSct(:,3)),[],c(4,:),'x');
plot(vfAve(:,:,4),aveMaxED(:,:,4),'--k','linewidth',1.5) % Plot mean line
%plot(perryed(1,:), perryed(2,:),'.-.','linewidth',1.5)
ylabel("Maximum Energy Density [J/cm^3]");
xlabel("Volume Fraction of Nanoparticles [%]");
legend('Percolating Geometries','Near Percolating Geometries', 'Non-Percolating Geometries', 'Mean', 'Location', 'southeast')
 
%% BREAKDOWN FIELD STRENGTH PLOTS
figure(bdfs)
c = get(gca,'colororder') % variation of permittivity as scatterplot
scatter(nonzeros(btVFsct(:,3)),nonzeros(bdSct(:,3)).*1E-3,[],c(4,:),'x');
scatter(nonzeros(btVFsct(:,2)),nonzeros(bdSct(:,2)).*1E-3,[],c(4,:),'o');
scatter(nonzeros(btVFsct(:,1)),nonzeros(bdSct(:,1)).*1E-3,[],c(4,:),'.');

errorbar(nonzeros(vfAve(:,:,4)), nonzeros(avBD(1,:,4))'.*1E-3, 3*stdevBD(counter(1,:,4)~=0,4).*1E-3, '-k', 'linewidth',1);
%plot(perrybd(1,:), perrybd(2,:),'ok','linewidth',1.5)
ylabel("Breakdown Field Strength [MV/M]");
xlabel("Volume Fraction of Nanoparticles [%]");
legend('Percolating Geometries','Near Percolating Geometries', 'Non-Percolating Geometries', 'All Geometries','Measured Data (Perry)', 'Location', 'northeast')

figure(bdfs2)
plot(nonzeros(vfAve(:,:,4)), nonzeros(avBD(1,:,4))', '--k', 'linewidth',1.5);
plot(perrybd(1,:), perrybd(2,:),'ok','linewidth',1.5)
ylabel("Breakdown Field Strength [MV/M]");
xlabel("Volume Fraction of Nanoparticles[%]");
legend('Simulated Geometries','Measured Data (Perry)', 'Non-Percolating Geometries', 'All Geometries', 'Location', 'northeast')


%% MEAN E-FIELD PLOTS
figure(efld)
scatter(btVFsct(:,1),avESct(:,1)/Em2p,[],c(4,:),'.');
scatter(btVFsct(:,2),avESct(:,2)/Em2p,[],c(4,:),'o');
scatter(btVFsct(:,3),avESct(:,3)/Em2p,[],c(4,:),'x');
scatter(btVFsct(:,1),avENPSct(:,1)/Em2p,[],c(4,:),'.');
scatter(btVFsct(:,2),avENPSct(:,2)/Em2p,[],c(4,:),'o');
scatter(btVFsct(:,3),avENPSct(:,3)/Em2p,[],c(4,:),'x');
errorbar(nonzeros(vfAve(:,:,4)), nonzeros(aveE(1,:,4)/Em2p)', 3*stdevE(counter(1,:,4)~=0,4)/Em2p, '--k', 'linewidth',1);
errorbar(nonzeros(vfAve(:,:,4)), nonzeros(aveENP(1,:,4)/Em2p)', 3*stdevENP(counter(1,:,4)~=0,4)/Em2p, '--k', 'linewidth',1);
ylabel("Average Normalized Electric Field");
xlabel("Volume Fraction of Nanoparticles [%]");
%legend('Inside the Host Material','Inside the Nanoparticles', 'Location', 'southeast')

% figure(efldBT)
% scatter(btVFsct(:,1),avENPSct(:,1),[],c(1,:),'o');
% scatter(btVFsct(:,2),avENPSct(:,2),[],c(2,:),'o');
% scatter(btVFsct(:,3),avENPSct(:,3),[],c(4,:),'o');
% ylabel("Average Electric Field within Nanoparticles");
% xlabel("Volume Fraction of BT Nanoparticles");
% errorbar(nonzeros(vfAve(:,:,4)), nonzeros(aveENP(1,:,4))', stdevENP(counter(1,:,4)~=0,4), linetps(4), 'linewidth',2);

% legend('Percolating Geometries','Near Percolating Geometries', 'Non-Percolating Geometries', 'All Geometries', 'Location', 'southeast')
%% OTHER MISC. PLOTS
% max V vs. permittivity
figure(volts)
%errorbar(nonzeros(vfAve(:,:,4)), nonzeros(avBD(1,:,4))', 3*stdevBD(counter(1,:,4)~=0,4), linetps(4), 'linewidth',1);
c = get(gca,'colororder') % variation of permittivity as scatterplot
scatter(nonzeros(btVFsct(:,1)),nonzeros(MaxVolt2p(:,1)),[],c(4,:),'o');
scatter(nonzeros(btVFsct(:,2)),nonzeros(MaxVolt2p(:,2)),[],c(4,:),'o');
scatter(nonzeros(btVFsct(:,3)),nonzeros(MaxVolt2p(:,3)),[],c(4,:),'o');

scatter(nonzeros(btVFsct(:,1)),nonzeros(maxVolt2pNP(:,1)),[],c(4,:),'.');
scatter(nonzeros(btVFsct(:,2)),nonzeros(maxVolt2pNP(:,2)),[],c(4,:),'.');
scatter(nonzeros(btVFsct(:,3)),nonzeros(maxVolt2pNP(:,3)),[],c(4,:),'.');
ylabel("Maximum Voltage [V]");
xlabel("Volume Fraction of Nanoparticles [%]");
legend('Percolating in Host','Near Percolating in Host', 'Non-Percolating in Host', 'Percolating in NP','Near Percolating in NP', 'Non-Percolating in NP', 'Location', 'southeast')

figure(volts2)
%errorbar(nonzeros(vfAve(:,:,4)), nonzeros(avBD(1,:,4))', 3*stdevBD(counter(1,:,4)~=0,4), linetps(4), 'linewidth',1);
c = get(gca,'colororder') % variation of permittivity as scatterplot
scatter(nonzeros(effSct(:,1)),nonzeros(MaxVolt2p(:,1)),[],c(4,:),'o');
scatter(nonzeros(effSct(:,2)),nonzeros(MaxVolt2p(:,2)),[],c(4,:),'o');
scatter(nonzeros(effSct(:,3)),nonzeros(MaxVolt2p(:,3)),[],c(4,:),'o');

scatter(nonzeros(effSct(:,1)),nonzeros(maxVolt2pNP(:,1)),[],c(4,:),'.');
scatter(nonzeros(effSct(:,2)),nonzeros(maxVolt2pNP(:,2)),[],c(4,:),'.');
scatter(nonzeros(effSct(:,3)),nonzeros(maxVolt2pNP(:,3)),[],c(4,:),'.');
ylabel("Maximum Voltage [V]");
xlabel("Effective Permittivity [F/M]");
legend('Percolating in Host','Near Percolating in Host', 'Non-Percolating in Host','Percolating in NP','Near Percolating in NP', 'Non-Percolating in NP', 'Location', 'northeast')


%% COUNTER BAR PLOT
figure(countbars)
allcounts = [counter(:,:,1)'-1, counter(:,:,2)'-1, counter(:,:,3)'-1];% counter(:,:,4)'];
%hb = bar(vfAve(:,:,4)', allcounts, 'stacked');
ylabel("Number of Microstructures");
xlabel("Volume Fraction of Nanoparticles [%]");
ylim([0,80]);
legend('Percolaging Geometries', 'Near Percolating Geometries', 'Non-Percolating Geometries')
saveas(countbars,'bars2','fig');
saveas(countbars,'bars2','jpeg');
saveas(perm,'permValid2','fig');
saveas(perm,'permValid2','jpeg');
saveas(eden,'eden2','fig');
saveas(eden,'eden2','jpeg');
saveas(efldBT,'efldBT2','fig');
saveas(efldBT,'efldBT2','jpeg');
saveas(bdfs,'bdfs2','fig');
saveas(bdfs,'bdfs2','jpeg');
saveas(efld,'efld_hst2','fig');
saveas(efld,'efld_hst2','jpeg');
saveas(perm2,'perm3','fig');
saveas(perm2,'perm3','jpeg');
saveas(eden2,'eden3','fig');
saveas(eden2,'eden3','jpeg');
saveas(efldBT2,'efldBT3','fig');
saveas(efldBT2,'efldBT3','jpeg');
saveas(bdfs2,'bdfs3','fig');
saveas(bdfs2,'bdfs3','jpeg');
saveas(efld2,'efld_hst3','fig');
saveas(efld2,'efld_hst3','jpeg');
saveas(volts,'volts1','fig');
saveas(volts,'volts1','jpeg');
saveas(volts2,'volts2','fig');
saveas(volts2,'volts2','jpeg');

figure (MCP)
hold on
%MCPdat = [3.9,3.9,6.3,6.7,7.1,7.1,8.5,8.5,21.0,21.0,21.0,19.0,19.0,23.4,26,26,26.8,29,31,60,95,95,95,95,50,98;1400,1310,630,870,760,920,620,720,470,400,300,250,600,350,320,420,420,380,400,320,300,250,200,100,330]
%Mcpdat = [3.9,1400;3.9,1310;6.3,630;6.7,870;7.1,760;7.1,920;8.5,620;8.5,720;21.0,470;21.0,400;21.0,300;19.0,250;19.0,600;23.4,350;26,320;26,420;26.8,420;29,380;31,400;60,320;95,300;95,250;95,200;95,100;50,300;98,330;]
scatter(effSct(:,3),bdSct(:,3).^2,[],c(4,:),'x');
scatter(effSct(:,2),bdSct(:,2).^2,[],c(4,:),'o');
scatter(effSct(:,1),bdSct(:,1).^2,[],c(4,:),'.');
c = get(gca,'colororder') % variation of permittivity as scatterplot

%plot(Mcpdat(:,1), Mcpdat(:,2).^2, 'ok', 'linewidth',1);
ylabel("Breakdown Field Strength^2 [MV/M]");
xlabel("Effective Relative Permittivity [F/M]");
legend('Percolating Geometries','Near Percolating Geometries', 'Non-Percolating Geometries', 'Data from ', 'Location', 'southeast')
set(gca,'xscale','log')
set(gca,'yscale','log')

% Get a list of all variables
allvars = whos;
% Identify the variables that ARE NOT graphics handles. This uses a regular
% expression on the class of each variable to check if it's a graphics object
tosave = cellfun(@isempty, regexp({allvars.class}, '^matlab\.(ui|graphics)\.'));
% Pass these variable names to save
save('dataset2_nf.mat', allvars(tosave).name)


%% Unpack & Plot Large Datasets
% clc
% clear all
PSeeds = [3,31,40,47];
NPseeds = [12,23,24,37,53,59,62,67,68,71,77,79];
therest = [13,14,18,41,48,55,66,73];
BTbd = 1e6;
dirnm = 'PercRuns/Run_1034/';

perrydat = [0,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6; 11.9, 12.5, 15.2, 16.4, 18.6, 22.3, 26.2, 35.1, 37.5]
perrybd = [0,0.05,0.1,0.15,0.2,0.3,0.4,0.5;380,344,345,300,235,230,218,213];
perryed = [0,0.1,0.2,0.3,0.4,0.5;1.6E6,1.6E6,2E6,2.3E6,2.35E6,3.2E6];
% perm = figure;
% hold on
% efld = figure;
% hold on
% eden = figure;
% hold on
% bdfs = figure;
% hold on
% efldBT = figure;
% hold on
% countbars = figure;
% hold on
% perm2 = figure;
% hold on
% efld2 = figure;
% hold on
% eden2 = figure;
% hold on
% bdfs2 = figure;
% hold on
% efldBT2 = figure;
% hold on
% volts = figure;
% hold on
% volts2 = figure;
% hold on


Em2p = (1)/(65*2e-9)
vfs = [0.05:0.05:0.6]

vl = length(vfs)
files = dir(dirnm);
effchart = zeros(vl,100,4);
avEchart = effchart;
avENPchart = effchart;
bdchart = effchart;
brugchart = effchart;
maxEdchart = effchart;
[numF,~] = size(files);
numF = numF -2;
files = files(3:numF+2);
vfAve = zeros(1,length(vfs),4);
effAve = vfAve;
aveMaxED = vfAve;
counter = vfAve+1;
countvf = vfAve;
aveE = vfAve;
aveENP = vfAve;
avBD = vfAve;


clear VFrac
clear eff
clear btVFr
clear bdv2pBT
clear BTstdE
clear MaxVolt2p
clear maxVolt2pNP
clear MaxEd2pNP
inv = 0;
alltheVF = zeros(1,numF);
for j = 2:numF
    Fn = j - inv ;
    fnm = files(Fn).name;
    fnm = [dirnm,fnm];
    importfile(fnm)
    if length(btVF) < 2
        alltheVF(Fn-1) = btVF;
    else
        files(Fn) = [];
        inv = inv + 1;
    end
end
alltheVF = alltheVF(1:(numF-inv-1))
[~,VFsort]=sort(alltheVF); %Get the order of B
files = files(VFsort);

brg = zeros(1,length(vfs))
btVFr = zeros(numF,3);
btVFsct = btVFr;
effSct = btVFr;
avESct = btVFr;
avENPSct = btVFr;
bdSct = btVFr;
maxEdSct = btVFr;
effSctS = btVFr;
avESctS = btVFr;
avENPSctS = btVFr;
bdSctS = btVFr;
maxEdSctS = btVFr;
              
              
for j = 2:numF-1
    i = j;
    fnm = files(i).name;
    fnm = [dirnm,fnm];
    if length(fnm) > 30
        importfile(fnm)
        
        if(btVF > 0)
            
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
            if p~=4
                % Scale Things...
                eff2p = eff2p*2;
%                 Emean2p = Emean2p;
%                 Emean2pNP = Emean2pNP;
                PVDFbd = PVDFbd*1e-6;
                BTbd = BTbd*1e-6;
                % ^^ fix later in Tester3 so that variation (stdev) is 
                % correctly calculated
                
                
                eff(j,1,p) = eff2p;
                eff(j,2,p) = eff3p;
                eff(j,3,p) = eff3pH;
                eff(j,4,p) = eff3pS;
                % use path to compute average electric field along breakdown path
                stE2pNP = std(std(std(E2p(OmegaBT(:,:,1:gd-1)==2)))); %.*(Delta(1)*Delta(2)*Delta(3));
                bdv2pBT(j,p) = (gd*BTbd)./((Emean2pNP+stE2p));
                maxVolt2pNP(j,p) =   BTbd ./ ( Emean2pNP + stE2pNP );
                MaxEd2pNP(j,p) = Edens2p * maxVolt2pNP(j,p).^2;
                MaxVolt2p(j,p) =  PVDFbd ./ ( Emean2p + stE2p );
                avE(j,1,p) = avE2p;
                
                Ed = (E2p.^2 .*(1/2)*e0 .* Mats(OmegaBT(:,:,1:Zsz-1)+1)); % Delta(1)*Delta(2)*Delta(3)
                Ed = sum(Ed,3);
                Ed = sum(Ed,2);
                Ed2p = sum(Ed,1);
                Edens2p(vfi) = Ed2p/(gd*gd*gd);
                
                % mcE(vfs == volfrac) = mcE
                % maxEd(j) = MaxEd2p;
                if ((maxVolt2pNP(j,p) < MaxVolt2p(j,p)) && p~=3)
                    maxEd(j,p) = (maxVolt2pNP(j,p)*1E6)^2*Edens2p;
                    % FIX THIS --->> maxEdS(j,p) = (maxVolt3pNPS(j,p)*1E6)^2*Edens3pS;
                    bd = maxVolt2pNP(j,p)/(Delta(3)*gd);%BTbd/((Emean2pNP+stE2pNP)/(3.1E6)); % Approximate y intercept
%                     bdS = maxVolt2pNPS(j,p)/(Delta(3)*gd);%BTbd/((Emean2pNP+stE2pNP)/(3.1E6)); % Approximate y intercept
                     % Account for breakdown voltage inside nanoparticle
                     % HERE
                else
                    maxEd(j,p) = (MaxVolt2p(j,p)*1E6)^2*Edens2p;
%                     maxEdS(j,p) = (maxVolt3pS(j,p)*1E6)^2*Edens3pS;
                    bd = MaxVolt2p(j,p)/(Delta(3)*gd);%PVDFbd/((Emean2p+stE2p)/(7.8E6)); % divide by approximate y intercept...
%                     bdS = MaxVolt3pS(j,p)/(Delta(3)*gd);%PVDFbd/((Emean2p+stE2p)/(7.8E6)); % divide by approximate y intercept...

                    % AND HERE
                end
                volfracs = round(vfs*1000) == round(volfrac*1000);
                disp(volfracs)
                eff2p = eff2p/2; % Adjust permittivity
                aveE(1,volfracs,p) = aveE(1,volfracs,p) + Emean2p;
                aveENP(1,volfracs,p) = aveENP(1,volfracs,p) + Emean2pNP;
                counter(1,volfracs,p) = counter(1,volfracs,p) + 1;
                effAve(1,volfracs,p) = effAve(1,volfracs,p) + eff2p;
                avBD(1,volfracs,p) = avBD(1,volfracs,p) + bd;
                aveMaxED(1,volfracs,p) = aveMaxED(1,volfracs,p) + maxEd(j,p);
                effchart(volfracs, counter(1,volfracs,p),p) = eff2p;
                avEchart(volfracs, counter(1,volfracs,p),p) = Emean2p;
                avENPchart(volfracs, counter(1,volfracs,p),p) = Emean2pNP;
                bdchart(volfracs, counter(1,volfracs,p),p) = bd;
                maxEdchart(volfracs, counter(1,volfracs,p),p) = maxEd(j,p);
%                 
%                 eff3pS = eff3pS/2; % Adjust permittivity
%                 aveES(1,volfracs,p) = aveES(1,volfracs,p) + Emean3pS;
%                 aveENPS(1,volfracs,p) = aveENPS(1,volfracs,p) + Emean3pNPS;
%                 counterS(1,volfracs,p) = counterS(1,volfracs,p) + 1;
%                 effAveS(1,volfracs,p) = effAveS(1,volfracs,p) + eff3pS;
%                 avBDS(1,volfracs,p) = avBDS(1,volfracs,p) + bd;
%                 aveMaxEDS(1,volfracs,p) = aveMaxEDS(1,volfracs,p) + maxEdS(j,p);
%                 effchartS(volfracs, counterS(1,volfracs,p),p) = eff3pS;
%                 avEchartS(volfracs, counterS(1,volfracs,p),p) = Emean3pS;
%                 avENPchartS(volfracs, counterS(1,volfracs,p),p) = Emean3pNPS;
%                 bdchartS(volfracs, counterS(1,volfracs,p),p) = bdS;
%                 maxEdchartS(volfracs, counterS(1,volfracs,p),p) = maxEdS(j,p);
                
                vfchart(volfracs, counter(1,volfracs,p),p) = btVF;
                
                cSct = sum(counter-1,2)
                btVFr(   cSct(p),p) = btVF;
                btVFsct( cSct(p),p) = btVF;
                effSct(  cSct(p),p) = eff2p;
                avESct(  cSct(p),p) = Emean2p;
                avENPSct(cSct(p),p) = Emean2pNP;
                bdSct(   cSct(p),p) = bd;
                maxEdSct(j,p) = maxEd(j,p);
                                
%                 cSctS = sum(counterS-1,2)
%                 btVFrS(   cSct(p),p) = btVF;
%                 btVFsctS( cSct(p),p) = btVF;
%                 effSctS(  cSctS(p),p) = eff3pS;
%                 avESctS(  cSctS(p),p) = Emean3pS;
%                 avENPSctS(cSctS(p),p) = Emean3pNPS;
%                 bdSctS(   cSctS(p),p) = bdS;
%                 maxEdSctS(j,p) = maxEdS(j,p);
                
                
                vfAve(1,volfracs,p) = vfAve(1,volfracs,p) + btVF;
                vfAve(1,volfracs,4) = vfAve(1,volfracs,4) + btVF;
%                 vfAveS(1,volfracs,p) = vfAveS(1,volfracs,p) + btVF;
%                 vfAveS(1,volfracs,4) = vfAveS(1,volfracs,4) + btVF;
                countvf(1,volfracs,p) = countvf(1,volfracs,p) + 1;
                
                countvf(1,volfracs,4) = countvf(1,volfracs,4) + 1;
                counter(1,volfracs,4) = counter(1,volfracs,4) + 1;
%                 effchartS(volfracs, counter(1,volfracs,4),4) = eff2p;
%                 avEchartS(volfracs, counter(1,volfracs,4),4) = Emean2p;
%                 avENPchartS(volfracs, counter(1,volfracs,4),4) = Emean2pNP;
%                 bdchartS(volfracs, counter(1,volfracs,4),4) = bd;
%                 maxEdchartS(volfracs, counter(1,volfracs,4),4) = maxEd(j,p);
                
                                syms x % effective permittivity of solid
            if((phases==2) && (brg(mci,volfracs)==0))
                brgEq = (1-volfrac)*((CxPvty(1)-x)/(CxPvty(1)+x)) + (volfrac)*((CxPvty(3)-x)/(CxPvty(3)+x)) == 0;
            elseif((phases == 3)  && (brg(mci,volfracs)==0))
                brgEq = (1-volfrac)*((CxPvty(1)-x)/(CxPvty(1)+x)) + (volfrac)*((CxPvty(2)-x)/(CxPvty(2)+x)) + (volfrac2)*((CxPvty(3)-x)/(CxPvty(3)+x))==0;
            end
            soln = solve(brgEq,x);
            brg(mci,volfracs) = max(soln);


            end
        end
    end
end

for bcnt = 1:length(volfracs)
    brgEq = (1-vfs(bcnt))*((Mats(1)-x)/(Mats(1)+x)) + (vfs(bcnt))*((Mats(3)-x)/(Mats(3)+x)) == 0;
    soln = solve(brgEq,x);
    brugg = max(soln);
    brug(bcnt) = brugg;
end
%     figure(perm);
%     plot(vfs, brug, '--', 'linewidth',2);
    
% line identifiers
linetps = ['-o';'-x';'-*';'-v'];
dottps = ['o';'x';'*';'v'];
for i = 1:3
    ctemp = nonzeros(counter(1,:,i)-1)
    l = length(ctemp)
    
    effAve(1,1:l,i) = nonzeros(effAve(1,:,i))./ctemp;
    effAve(1,l+1:vl,i) = 0;
    aveE(1,1:l,i) = nonzeros(aveE(1,:,i))./ctemp;
    aveE(1,l+1:vl,i) = 0;
    aveENP(1,1:l,i) = nonzeros(aveENP(1,:,i))./ctemp;
    aveENP(1,l+1:vl,i) = 0;
    avBD(1,1:l,i) = nonzeros(avBD(1,:,i))./ctemp;
    avBD(1,l+1:vl,i) = 0;
    vfAve(1,1:l,i)  = nonzeros(vfAve(1,:,i))./ctemp;
    vfAve(1,l+1:vl,i) = 0;
    aveMaxED(1,1:l,i) = nonzeros(aveMaxED(1,:,i))./ctemp;
    aveMaxED(1,l+1:vl,i) = 0;
    for num = 1:length(vfs)
        stdev(num,i) = std(nonzeros(effchart(num,1:nonzeros(counter(1,num,i)-1),i)));
        stdevE(num,i) = std(nonzeros(avEchart(num,1:nonzeros(counter(1,num,i)-1),i)));
        stdevENP(num,i) = std(nonzeros(avENPchart(num,1:nonzeros(counter(1,num,i)-1),i)));
        stdevBD(num,i) = std(nonzeros(bdchart(num,1:nonzeros(counter(1,num,i)-1),i)))
        stdevMED(num,i) = std(nonzeros(maxEdchart(num,1:nonzeros(counter(1,num,i)-1),i)));
        
        stdev(num,4) = std(nonzeros(effchart(num,1:nonzeros(counter(1,num,4)-1),4)));
        stdevE(num,4) = std(nonzeros(avEchart(num,1:nonzeros(counter(1,num,4)-1),4)));
        stdevENP(num,4) = std(nonzeros(avENPchart(num,1:nonzeros(counter(1,num,4)-1),4)));
        stdevBD(num,4) = std(nonzeros(bdchart(num,1:nonzeros(counter(1,num,4)-1),4)))
        stdevMED(num,4) = std(nonzeros(maxEdchart(num,1:nonzeros(counter(1,num,4)-1),4)));
        
%         stdevS(num,i) = std(nonzeros(effchartS(num,1:nonzeros(counterS(1,num,i)-1),i)));
%         stdevES(num,i) = std(nonzeros(avEchartS(num,1:nonzeros(counterS(1,num,i)-1),i)));
%         stdevENPS(num,i) = std(nonzeros(avENPchartS(num,1:nonzeros(counterS(1,num,i)-1),i)));
%         stdevBDS(num,i) = std(nonzeros(bdchartS(num,1:nonzeros(counterS(1,num,i)-1),i)))
%         stdevMEDS(num,i) = std(nonzeros(maxEdchartS(num,1:nonzeros(counterS(1,num,i)-1),i)));
%         
%         stdevS(num,4) = std(nonzeros(effchartS(num,1:nonzeros(counterS(1,num,4)-1),4)));
%         stdevES(num,4) = std(nonzeros(avEchartS(num,1:nonzeros(counterS(1,num,4)-1),4)));
%         stdevENPS(num,4) = std(nonzeros(avENPchartS(num,1:nonzeros(counterS(1,num,4)-1),4)));
%         stdevBDS(num,4) = std(nonzeros(bdchartS(num,1:nonzeros(counterS(1,num,4)-1),4)))
%         stdevMEDS(num,4) = std(nonzeros(maxEdchartS(num,1:nonzeros(counterS(1,num,4)-1),4)));
        

    end
    
%     figure(perm);
%     errorbar(nonzeros(vfAve(:,:,i)), nonzeros(effAve(1,:,i))', 3*stdev(counter(1,:,i)-1~=0,i), dottps(i), 'linewidth',1.5);
%     
%     figure(efld)
%     errorbar(nonzeros(vfAve(:,:,i)), nonzeros(aveE(1,:,i))', 3*stdevE(counter(1,:,i)-1~=0,i), linetps(i), 'linewidth',1.5);
% 
%     figure(efld2)
%     errorbar(nonzeros(vfAve(:,:,i)), nonzeros(aveE(1,:,i))', 3*stdevE(counter(1,:,i)-1~=0,i), dottps(i), 'linewidth',1.5);
%     
%     figure(efldBT)
%     errorbar(nonzeros(vfAve(:,:,i)), nonzeros(aveENP(1,:,i))', 3*stdevENP(counter(1,:,i)-1~=0,i), dottps(i), 'linewidth',1.5);
%     
%     figure(eden);
%     errorbar(nonzeros(vfAve(:,:,i)), nonzeros(aveMaxED(1,:,i))', 3*stdevMED(counter(1,:,i)-1~=0,i), dottps(i),'linewidth',1.5);
%     
%     figure(bdfs);
%     errorbar(nonzeros(vfAve(:,:,i)), nonzeros(avBD(1,:,i))', stdevBD(counter(1,:,i)-1~=0,i), dottps(i),'linewidth',1.5);
%     
end
ctemp = nonzeros(counter(1,:,4)-1);
l = length(ctemp)

ctemp = nonzeros(counter(1,:,4)-1);
effAve(1,1:l,4) = nonzeros(sum(sum(effchart(:,:,1:3),3),2))./ctemp;
effAve(1,l+1:vl,4) = 0;
aveE(1,1:l,4) = nonzeros(sum(sum(avEchart(:,:,1:3),2),3))./ctemp;
aveE(1,l+1:vl,4) = 0;
aveENP(1,1:l,4) = nonzeros(sum(sum(avENPchart(:,:,1:3),2),3))./ctemp;
aveENP(1,l+1:vl,4) = 0;
avBD(1,1:l,4) = nonzeros(sum(sum(bdchart(:,:,1:3),2),3))./ctemp;
avBD(1,l+1:vl,4) = 0;
aveMaxED(1,1:l,4) = nonzeros(sum(sum(maxEdchart(:,:,1:3),2),3))./ctemp;
aveMaxED(1,l+1:vl,4) = 0;

vfAve(1,1:l,4)  = nonzeros(vfAve(1,:,4))./ctemp;
vfAve(1,l+1:vl,4) = 0;

counter(1,:,4) = sum(counter(1,:,:)-1,3);

%------------------for loop ended here ------------%
%% EFFECTIVE PERMITTIVITY PLOTS
figure(perm)
c = get(gca,'colororder') % variation of permittivity as scatterplot
scatter(btVFsct(:,1),effSct(:,1),[],c(5,:),'.');
scatter(btVFsct(:,2),effSct(:,2),[],c(5,:),'o');
scatter(btVFsct(:,3),effSct(:,3),[],c(5,:),'x');
errorbar(nonzeros(vfAve(:,:,4)), nonzeros(effAve(1,:,4))', 3*stdev(counter(1,:,4)~=0,4), '--k', 'linewidth',1);
ylabel("Effective Relative Permittivity [F/M]");
xlabel("Volume Fraction of Nanoparticles [%]");
legend('Percolating Geometries','Near Percolating Geometries', 'Non-Percolating Geometries', 'All Geometries', 'Location', 'southeast')

figure(perm2)% mean of permittivity with perry and bruggeman
plot(nonzeros(vfAve(:,:,4)), nonzeros(effAve(1,:,4))', '--k', 'linewidth',1.5);
plot(perrydat(1,:), perrydat(2,:),'ok','linewidth',1.5)
%plot(vfs, brug,'--k','linewidth',1.5)
ylabel("Effective Relative Permittivity [F/M]");
xlabel("Volume Fraction of Nanoparticles [%]");
legend('Mean of tested geometries','Experimental Data (Perry 2009)', 'Bruggeman Model', 'All Geometries', 'Location', 'southeast')



%% ENERGY DENSITY  PLOTS

figure(eden)
errorbar(nonzeros(vfAve(:,:,4)), nonzeros(aveMaxED(1,:,4))', 3*stdevMED(counter(1,:,4)~=0,4), dottps(4), 'linewidth',2);

ylabel("Maximum Energy Density [J/cm^3]");
xlabel("Volume Fraction of Nanoparticles [%]");
legend('Percolating Geometries','Near Percolating Geometries', 'Non-Percolating Geometries', 'All Geometries', 'Location', 'southeast')

figure(eden2)
%plot(nonzeros(vfAve(:,:,4)), nonzeros(aveMaxED(1,:,4))', linetps(4), 'linewidth',1.5);
c = get(gca,'colororder')
scatter(nonzeros(btVFsct(:,1)),nonzeros(maxEdSct(:,1)),[],c(5,:),'.');
scatter(nonzeros(btVFsct(:,2)),nonzeros(maxEdSct(:,2)),[],c(5,:),'o');
scatter(nonzeros(btVFsct(:,3)),nonzeros(maxEdSct(:,3)),[],c(5,:),'x');
plot(vfAve(:,:,4),aveMaxED(:,:,4),'--k','linewidth',1.5) % Plot mean line
%plot(perryed(1,:), perryed(2,:),'.-.','linewidth',1.5)
ylabel("Maximum Energy Density [J/cm^3]");
xlabel("Volume Fraction of Nanoparticles [%]");
legend('Percolating Geometries','Near Percolating Geometries', 'Non-Percolating Geometries', 'Mean', 'Location', 'southeast')
 
%% BREAKDOWN FIELD STRENGTH PLOTS
figure(bdfs)
c = get(gca,'colororder') % variation of permittivity as scatterplot
scatter(nonzeros(btVFsct(:,3)),nonzeros(bdSct(:,3)).*1E-3,[],c(5,:),'x');
scatter(nonzeros(btVFsct(:,2)),nonzeros(bdSct(:,2)).*1E-3,[],c(5,:),'o');
scatter(nonzeros(btVFsct(:,1)),nonzeros(bdSct(:,1)).*1E-3,[],c(5,:),'.');

errorbar(nonzeros(vfAve(:,:,4)), nonzeros(avBD(1,:,4))'.*1E-3, 3*stdevBD(counter(1,:,4)~=0,4).*1E-3, '-k', 'linewidth',1);
%plot(perrybd(1,:), perrybd(2,:),'ok','linewidth',1.5)
ylabel("Breakdown Field Strength [MV/M]");
xlabel("Volume Fraction of Nanoparticles [%]");
legend('Percolating Geometries','Near Percolating Geometries', 'Non-Percolating Geometries', 'All Geometries','Measured Data (Perry)', 'Location', 'northeast')

figure(bdfs2)
plot(nonzeros(vfAve(:,:,4)), nonzeros(avBD(1,:,4))', '--k', 'linewidth',1.5);
plot(perrybd(1,:), perrybd(2,:),'ok','linewidth',1.5)
ylabel("Breakdown Field Strength [MV/M]");
xlabel("Volume Fraction of Nanoparticles[%]");
legend('Simulated Geometries','Measured Data (Perry)', 'Non-Percolating Geometries', 'All Geometries', 'Location', 'northeast')


%% MEAN E-FIELD PLOTS
figure(efld)
scatter(btVFsct(:,1),avESct(:,1)/Em2p,[],c(5,:),'.');
scatter(btVFsct(:,2),avESct(:,2)/Em2p,[],c(5,:),'o');
scatter(btVFsct(:,3),avESct(:,3)/Em2p,[],c(5,:),'x');
scatter(btVFsct(:,1),avENPSct(:,1)/Em2p,[],c(5,:),'.');
scatter(btVFsct(:,2),avENPSct(:,2)/Em2p,[],c(5,:),'o');
scatter(btVFsct(:,3),avENPSct(:,3)/Em2p,[],c(5,:),'x');
errorbar(nonzeros(vfAve(:,:,4)), nonzeros(aveE(1,:,4)/Em2p)', 3*stdevE(counter(1,:,4)~=0,4)/Em2p, '--k', 'linewidth',1);
errorbar(nonzeros(vfAve(:,:,4)), nonzeros(aveENP(1,:,4)/Em2p)', 3*stdevENP(counter(1,:,4)~=0,4)/Em2p, '--k', 'linewidth',1);
ylabel("Average Normalized Electric Field");
xlabel("Volume Fraction of Nanoparticles [%]");
%legend('Inside the Host Material','Inside the Nanoparticles', 'Location', 'southeast')

% figure(efldBT)
% scatter(btVFsct(:,1),avENPSct(:,1),[],c(1,:),'o');
% scatter(btVFsct(:,2),avENPSct(:,2),[],c(2,:),'o');
% scatter(btVFsct(:,3),avENPSct(:,3),[],c(5,:),'o');
% ylabel("Average Electric Field within Nanoparticles");
% xlabel("Volume Fraction of BT Nanoparticles");
% errorbar(nonzeros(vfAve(:,:,4)), nonzeros(aveENP(1,:,4))', stdevENP(counter(1,:,4)~=0,4), linetps(4), 'linewidth',2);

% legend('Percolating Geometries','Near Percolating Geometries', 'Non-Percolating Geometries', 'All Geometries', 'Location', 'southeast')
%% OTHER MISC. PLOTS
% max V vs. permittivity
figure(volts)
%errorbar(nonzeros(vfAve(:,:,4)), nonzeros(avBD(1,:,4))', 3*stdevBD(counter(1,:,4)~=0,4), linetps(4), 'linewidth',1);
c = get(gca,'colororder') % variation of permittivity as scatterplot
scatter(nonzeros(btVFsct(:,1)),nonzeros(MaxVolt2p(:,1)),[],c(5,:),'o');
scatter(nonzeros(btVFsct(:,2)),nonzeros(MaxVolt2p(:,2)),[],c(5,:),'o');
scatter(nonzeros(btVFsct(:,3)),nonzeros(MaxVolt2p(:,3)),[],c(5,:),'o');

scatter(nonzeros(btVFsct(:,1)),nonzeros(maxVolt2pNP(:,1)),[],c(5,:),'.');
scatter(nonzeros(btVFsct(:,2)),nonzeros(maxVolt2pNP(:,2)),[],c(5,:),'.');
scatter(nonzeros(btVFsct(:,3)),nonzeros(maxVolt2pNP(:,3)),[],c(5,:),'.');
ylabel("Maximum Voltage [V]");
xlabel("Volume Fraction of Nanoparticles [%]");
legend('Percolating in Host','Near Percolating in Host', 'Non-Percolating in Host', 'Percolating in NP','Near Percolating in NP', 'Non-Percolating in NP', 'Location', 'southeast')

figure(volts2)
%errorbar(nonzeros(vfAve(:,:,4)), nonzeros(avBD(1,:,4))', 3*stdevBD(counter(1,:,4)~=0,4), linetps(4), 'linewidth',1);
c = get(gca,'colororder') % variation of permittivity as scatterplot
scatter(nonzeros(effSct(:,1)),nonzeros(MaxVolt2p(:,1)),[],c(5,:),'o');
scatter(nonzeros(effSct(:,2)),nonzeros(MaxVolt2p(:,2)),[],c(5,:),'o');
scatter(nonzeros(effSct(:,3)),nonzeros(MaxVolt2p(:,3)),[],c(5,:),'o');

scatter(nonzeros(effSct(:,1)),nonzeros(maxVolt2pNP(:,1)),[],c(5,:),'.');
scatter(nonzeros(effSct(:,2)),nonzeros(maxVolt2pNP(:,2)),[],c(5,:),'.');
scatter(nonzeros(effSct(:,3)),nonzeros(maxVolt2pNP(:,3)),[],c(5,:),'.');
ylabel("Maximum Voltage [V]");
xlabel("Effective Permittivity [F/M]");
legend('Percolating in Host','Near Percolating in Host', 'Non-Percolating in Host','Percolating in NP','Near Percolating in NP', 'Non-Percolating in NP', 'Location', 'northeast')


%% COUNTER BAR PLOT
figure(countbars)
allcounts = [counter(:,:,1)'-1, counter(:,:,2)'-1, counter(:,:,3)'-1];% counter(:,:,4)'];
%hb = bar(vfAve(:,:,4)', allcounts, 'stacked');
ylabel("Number of Microstructures");
xlabel("Volume Fraction of Nanoparticles [%]");
ylim([0,80]);
legend('Percolaging Geometries', 'Near Percolating Geometries', 'Non-Percolating Geometries')
saveas(countbars,'bars2','fig');
saveas(countbars,'bars2','jpeg');
saveas(perm,'permValid2','fig');
saveas(perm,'permValid2','jpeg');
saveas(eden,'eden2','fig');
saveas(eden,'eden2','jpeg');
saveas(efldBT,'efldBT2','fig');
saveas(efldBT,'efldBT2','jpeg');
saveas(bdfs,'bdfs2','fig');
saveas(bdfs,'bdfs2','jpeg');
saveas(efld,'efld_hst2','fig');
saveas(efld,'efld_hst2','jpeg');
saveas(perm2,'perm3','fig');
saveas(perm2,'perm3','jpeg');
saveas(eden2,'eden3','fig');
saveas(eden2,'eden3','jpeg');
saveas(efldBT2,'efldBT3','fig');
saveas(efldBT2,'efldBT3','jpeg');
saveas(bdfs2,'bdfs3','fig');
saveas(bdfs2,'bdfs3','jpeg');
saveas(efld2,'efld_hst3','fig');
saveas(efld2,'efld_hst3','jpeg');
saveas(volts,'volts1','fig');
saveas(volts,'volts1','jpeg');
saveas(volts2,'volts2','fig');
saveas(volts2,'volts2','jpeg');

figure(MCP)
hold on
%MCPdat = [3.9,3.9,6.3,6.7,7.1,7.1,8.5,8.5,21.0,21.0,21.0,19.0,19.0,23.4,26,26,26.8,29,31,60,95,95,95,95,50,98;1400,1310,630,870,760,920,620,720,470,400,300,250,600,350,320,420,420,380,400,320,300,250,200,100,330]
Mcpdat = [3.9,1400;3.9,1310;6.3,630;6.7,870;7.1,760;7.1,920;8.5,620;8.5,720;21.0,470;21.0,400;21.0,300;19.0,250;19.0,600;23.4,350;26,320;26,420;26.8,420;29,380;31,400;60,320;95,300;95,250;95,200;95,100;50,300;98,330;]
scatter(effSct(:,3),bdSct(:,3).^2,[],c(5,:),'x');
scatter(effSct(:,2),bdSct(:,2).^2,[],c(5,:),'o');
scatter(effSct(:,1),bdSct(:,1).^2,[],c(5,:),'.');
c = get(gca,'colororder') % variation of permittivity as scatterplot

%plot(Mcpdat(:,1), Mcpdat(:,2).^2, 'ok', 'linewidth',1);
ylabel("Breakdown Field Strength^2 [MV/M]");
xlabel("Effective Relative Permittivity [F/M]");
legend('Percolating Geometries','Near Percolating Geometries', 'Non-Percolating Geometries', 'Data from ', 'Location', 'southeast')
set(gca,'xscale','log')
set(gca,'yscale','log')

% Get a list of all variables
allvars = whos;
% Identify the variables that ARE NOT graphics handles. This uses a regular
% expression on the class of each variable to check if it's a graphics object
tosave = cellfun(@isempty, regexp({allvars.class}, '^matlab\.(ui|graphics)\.'));
% Pass these variable names to save
save('dataset2_nf.mat', allvars(tosave).name)


%% Unpack & Plot Large Datasets
% clc
% clear all
PSeeds = [3,31,40,47];
NPseeds = [12,23,24,37,53,59,62,67,68,71,77,79];
therest = [13,14,18,41,48,55,66,73];
BTbd = 1e6;
dirnm = 'PercRuns/Run_1035/';

perrydat = [0,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6; 11.9, 12.5, 15.2, 16.4, 18.6, 22.3, 26.2, 35.1, 37.5]
perrybd = [0,0.05,0.1,0.15,0.2,0.3,0.4,0.5;380,344,345,300,235,230,218,213];
perryed = [0,0.1,0.2,0.3,0.4,0.5;1.6E6,1.6E6,2E6,2.3E6,2.35E6,3.2E6];
% perm = figure;
% hold on
% efld = figure;
% hold on
% eden = figure;
% hold on
% bdfs = figure;
% hold on
% efldBT = figure;
% hold on
% countbars = figure;
% hold on
% perm2 = figure;
% hold on
% efld2 = figure;
% hold on
% eden2 = figure;
% hold on
% bdfs2 = figure;
% hold on
% efldBT2 = figure;
% hold on
% volts = figure;
% hold on
% volts2 = figure;
% hold on


Em2p = (1)/(65*2e-9)
vfs = [0.05:0.05:0.6]

vl = length(vfs)
files = dir(dirnm);
effchart = zeros(vl,100,4);
avEchart = effchart;
avENPchart = effchart;
bdchart = effchart;
brugchart = effchart;
maxEdchart = effchart;
[numF,~] = size(files);
numF = numF -2;
files = files(3:numF+2);
vfAve = zeros(1,length(vfs),4);
effAve = vfAve;
aveMaxED = vfAve;
counter = vfAve+1;
countvf = vfAve;
aveE = vfAve;
aveENP = vfAve;
avBD = vfAve;


clear VFrac
clear eff
clear btVFr
clear bdv2pBT
clear BTstdE
clear MaxVolt2p
clear maxVolt2pNP
clear MaxEd2pNP
inv = 0;
alltheVF = zeros(1,numF);
for j = 2:numF
    Fn = j - inv ;
    fnm = files(Fn).name;
    fnm = [dirnm,fnm];
    importfile(fnm)
    if length(btVF) < 2
        alltheVF(Fn-1) = btVF;
    else
        files(Fn) = [];
        inv = inv + 1;
    end
end
alltheVF = alltheVF(1:(numF-inv-1))
[~,VFsort]=sort(alltheVF); %Get the order of B
files = files(VFsort);

brg = zeros(1,length(vfs))
btVFr = zeros(numF,3);
btVFsct = btVFr;
effSct = btVFr;
avESct = btVFr;
avENPSct = btVFr;
bdSct = btVFr;
maxEdSct = btVFr;
effSctS = btVFr;
avESctS = btVFr;
avENPSctS = btVFr;
bdSctS = btVFr;
maxEdSctS = btVFr;
              
              
for j = 2:numF-1
    i = j;
    fnm = files(i).name;
    fnm = [dirnm,fnm];
    if length(fnm) > 30
        importfile(fnm)
        
        if(btVF > 0)
            
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
            if p~=4
                % Scale Things...
                eff2p = eff2p*2;
%                 Emean2p = Emean2p;
%                 Emean2pNP = Emean2pNP;
                PVDFbd = PVDFbd*1e-6;
                BTbd = BTbd*1e-6;
                % ^^ fix later in Tester3 so that variation (stdev) is 
                % correctly calculated
                
                
                eff(j,1,p) = eff2p;
                eff(j,2,p) = eff3p;
                eff(j,3,p) = eff3pH;
                eff(j,4,p) = eff3pS;
                % use path to compute average electric field along breakdown path
                stE2pNP = std(std(std(E2p(OmegaBT(:,:,1:gd-1)==2)))); %.*(Delta(1)*Delta(2)*Delta(3));
                bdv2pBT(j,p) = (gd*BTbd)./((Emean2pNP+stE2p));
                maxVolt2pNP(j,p) =   BTbd ./ ( Emean2pNP + stE2pNP );
                MaxEd2pNP(j,p) = Edens2p * maxVolt2pNP(j,p).^2;
                MaxVolt2p(j,p) =  PVDFbd ./ ( Emean2p + stE2p );
                avE(j,1,p) = avE2p;
                
                Ed = (E2p.^2 .*(1/2)*e0 .* Mats(OmegaBT(:,:,1:Zsz-1)+1)); % Delta(1)*Delta(2)*Delta(3)
                Ed = sum(Ed,3);
                Ed = sum(Ed,2);
                Ed2p = sum(Ed,1);
                Edens2p(vfi) = Ed2p/(gd*gd*gd);
                
                % mcE(vfs == volfrac) = mcE
                % maxEd(j) = MaxEd2p;
                if ((maxVolt2pNP(j,p) < MaxVolt2p(j,p)) && p~=3)
                    maxEd(j,p) = (maxVolt2pNP(j,p)*1E6)^2*Edens2p;
                    % FIX THIS --->> maxEdS(j,p) = (maxVolt3pNPS(j,p)*1E6)^2*Edens3pS;
                    bd = maxVolt2pNP(j,p)/(Delta(3)*gd);%BTbd/((Emean2pNP+stE2pNP)/(3.1E6)); % Approximate y intercept
%                     bdS = maxVolt2pNPS(j,p)/(Delta(3)*gd);%BTbd/((Emean2pNP+stE2pNP)/(3.1E6)); % Approximate y intercept
                     % Account for breakdown voltage inside nanoparticle
                     % HERE
                else
                    maxEd(j,p) = (MaxVolt2p(j,p)*1E6)^2*Edens2p;
%                     maxEdS(j,p) = (maxVolt3pS(j,p)*1E6)^2*Edens3pS;
                    bd = MaxVolt2p(j,p)/(Delta(3)*gd);%PVDFbd/((Emean2p+stE2p)/(7.8E6)); % divide by approximate y intercept...
%                     bdS = MaxVolt3pS(j,p)/(Delta(3)*gd);%PVDFbd/((Emean2p+stE2p)/(7.8E6)); % divide by approximate y intercept...

                    % AND HERE
                end
                volfracs = round(vfs*1000) == round(volfrac*1000);
                disp(volfracs)
                eff2p = eff2p/2; % Adjust permittivity
                aveE(1,volfracs,p) = aveE(1,volfracs,p) + Emean2p;
                aveENP(1,volfracs,p) = aveENP(1,volfracs,p) + Emean2pNP;
                counter(1,volfracs,p) = counter(1,volfracs,p) + 1;
                effAve(1,volfracs,p) = effAve(1,volfracs,p) + eff2p;
                avBD(1,volfracs,p) = avBD(1,volfracs,p) + bd;
                aveMaxED(1,volfracs,p) = aveMaxED(1,volfracs,p) + maxEd(j,p);
                effchart(volfracs, counter(1,volfracs,p),p) = eff2p;
                avEchart(volfracs, counter(1,volfracs,p),p) = Emean2p;
                avENPchart(volfracs, counter(1,volfracs,p),p) = Emean2pNP;
                bdchart(volfracs, counter(1,volfracs,p),p) = bd;
                maxEdchart(volfracs, counter(1,volfracs,p),p) = maxEd(j,p);
%                 
%                 eff3pS = eff3pS/2; % Adjust permittivity
%                 aveES(1,volfracs,p) = aveES(1,volfracs,p) + Emean3pS;
%                 aveENPS(1,volfracs,p) = aveENPS(1,volfracs,p) + Emean3pNPS;
%                 counterS(1,volfracs,p) = counterS(1,volfracs,p) + 1;
%                 effAveS(1,volfracs,p) = effAveS(1,volfracs,p) + eff3pS;
%                 avBDS(1,volfracs,p) = avBDS(1,volfracs,p) + bd;
%                 aveMaxEDS(1,volfracs,p) = aveMaxEDS(1,volfracs,p) + maxEdS(j,p);
%                 effchartS(volfracs, counterS(1,volfracs,p),p) = eff3pS;
%                 avEchartS(volfracs, counterS(1,volfracs,p),p) = Emean3pS;
%                 avENPchartS(volfracs, counterS(1,volfracs,p),p) = Emean3pNPS;
%                 bdchartS(volfracs, counterS(1,volfracs,p),p) = bdS;
%                 maxEdchartS(volfracs, counterS(1,volfracs,p),p) = maxEdS(j,p);
                
                vfchart(volfracs, counter(1,volfracs,p),p) = btVF;
                
                cSct = sum(counter-1,2)
                btVFr(   cSct(p),p) = btVF;
                btVFsct( cSct(p),p) = btVF;
                effSct(  cSct(p),p) = eff2p;
                avESct(  cSct(p),p) = Emean2p;
                avENPSct(cSct(p),p) = Emean2pNP;
                bdSct(   cSct(p),p) = bd;
                maxEdSct(j,p) = maxEd(j,p);
                                
%                 cSctS = sum(counterS-1,2)
%                 btVFrS(   cSct(p),p) = btVF;
%                 btVFsctS( cSct(p),p) = btVF;
%                 effSctS(  cSctS(p),p) = eff3pS;
%                 avESctS(  cSctS(p),p) = Emean3pS;
%                 avENPSctS(cSctS(p),p) = Emean3pNPS;
%                 bdSctS(   cSctS(p),p) = bdS;
%                 maxEdSctS(j,p) = maxEdS(j,p);
                
                
                vfAve(1,volfracs,p) = vfAve(1,volfracs,p) + btVF;
                vfAve(1,volfracs,4) = vfAve(1,volfracs,4) + btVF;
%                 vfAveS(1,volfracs,p) = vfAveS(1,volfracs,p) + btVF;
%                 vfAveS(1,volfracs,4) = vfAveS(1,volfracs,4) + btVF;
                countvf(1,volfracs,p) = countvf(1,volfracs,p) + 1;
                
                countvf(1,volfracs,4) = countvf(1,volfracs,4) + 1;
                counter(1,volfracs,4) = counter(1,volfracs,4) + 1;
%                 effchartS(volfracs, counter(1,volfracs,4),4) = eff2p;
%                 avEchartS(volfracs, counter(1,volfracs,4),4) = Emean2p;
%                 avENPchartS(volfracs, counter(1,volfracs,4),4) = Emean2pNP;
%                 bdchartS(volfracs, counter(1,volfracs,4),4) = bd;
%                 maxEdchartS(volfracs, counter(1,volfracs,4),4) = maxEd(j,p);
                
                                syms x % effective permittivity of solid
            if((phases==2) && (brg(mci,volfracs)==0))
                brgEq = (1-volfrac)*((CxPvty(1)-x)/(CxPvty(1)+x)) + (volfrac)*((CxPvty(3)-x)/(CxPvty(3)+x)) == 0;
            elseif((phases == 3)  && (brg(mci,volfracs)==0))
                brgEq = (1-volfrac)*((CxPvty(1)-x)/(CxPvty(1)+x)) + (volfrac)*((CxPvty(2)-x)/(CxPvty(2)+x)) + (volfrac2)*((CxPvty(3)-x)/(CxPvty(3)+x))==0;
            end
            soln = solve(brgEq,x);
            brg(mci,volfracs) = max(soln);


            end
        end
    end
end

for bcnt = 1:length(volfracs)
    brgEq = (1-vfs(bcnt))*((Mats(1)-x)/(Mats(1)+x)) + (vfs(bcnt))*((Mats(3)-x)/(Mats(3)+x)) == 0;
    soln = solve(brgEq,x);
    brugg = max(soln);
    brug(bcnt) = brugg;
end
%     figure(perm);
%     plot(vfs, brug, '--', 'linewidth',2);
    
% line identifiers
linetps = ['-o';'-x';'-*';'-v'];
dottps = ['o';'x';'*';'v'];
for i = 1:3
    ctemp = nonzeros(counter(1,:,i)-1)
    l = length(ctemp)
    
    effAve(1,1:l,i) = nonzeros(effAve(1,:,i))./ctemp;
    effAve(1,l+1:vl,i) = 0;
    aveE(1,1:l,i) = nonzeros(aveE(1,:,i))./ctemp;
    aveE(1,l+1:vl,i) = 0;
    aveENP(1,1:l,i) = nonzeros(aveENP(1,:,i))./ctemp;
    aveENP(1,l+1:vl,i) = 0;
    avBD(1,1:l,i) = nonzeros(avBD(1,:,i))./ctemp;
    avBD(1,l+1:vl,i) = 0;
    vfAve(1,1:l,i)  = nonzeros(vfAve(1,:,i))./ctemp;
    vfAve(1,l+1:vl,i) = 0;
    aveMaxED(1,1:l,i) = nonzeros(aveMaxED(1,:,i))./ctemp;
    aveMaxED(1,l+1:vl,i) = 0;
    for num = 1:length(vfs)
        stdev(num,i) = std(nonzeros(effchart(num,1:nonzeros(counter(1,num,i)-1),i)));
        stdevE(num,i) = std(nonzeros(avEchart(num,1:nonzeros(counter(1,num,i)-1),i)));
        stdevENP(num,i) = std(nonzeros(avENPchart(num,1:nonzeros(counter(1,num,i)-1),i)));
        stdevBD(num,i) = std(nonzeros(bdchart(num,1:nonzeros(counter(1,num,i)-1),i)));
        stdevMED(num,i) = std(nonzeros(maxEdchart(num,1:nonzeros(counter(1,num,i)-1),i)));
        
        stdev(num,4) = std(nonzeros(effchart(num,1:nonzeros(counter(1,num,4)-1),4)));
        stdevE(num,4) = std(nonzeros(avEchart(num,1:nonzeros(counter(1,num,4)-1),4)));
        stdevENP(num,4) = std(nonzeros(avENPchart(num,1:nonzeros(counter(1,num,4)-1),4)));
        stdevBD(num,4) = std(nonzeros(bdchart(num,1:nonzeros(counter(1,num,4)-1),4)))
        stdevMED(num,4) = std(nonzeros(maxEdchart(num,1:nonzeros(counter(1,num,4)-1),4)));
        
%         stdevS(num,i) = std(nonzeros(effchartS(num,1:nonzeros(counterS(1,num,i)-1),i)));
%         stdevES(num,i) = std(nonzeros(avEchartS(num,1:nonzeros(counterS(1,num,i)-1),i)));
%         stdevENPS(num,i) = std(nonzeros(avENPchartS(num,1:nonzeros(counterS(1,num,i)-1),i)));
%         stdevBDS(num,i) = std(nonzeros(bdchartS(num,1:nonzeros(counterS(1,num,i)-1),i)))
%         stdevMEDS(num,i) = std(nonzeros(maxEdchartS(num,1:nonzeros(counterS(1,num,i)-1),i)));
%         
%         stdevS(num,4) = std(nonzeros(effchartS(num,1:nonzeros(counterS(1,num,4)-1),4)));
%         stdevES(num,4) = std(nonzeros(avEchartS(num,1:nonzeros(counterS(1,num,4)-1),4)));
%         stdevENPS(num,4) = std(nonzeros(avENPchartS(num,1:nonzeros(counterS(1,num,4)-1),4)));
%         stdevBDS(num,4) = std(nonzeros(bdchartS(num,1:nonzeros(counterS(1,num,4)-1),4)))
%         stdevMEDS(num,4) = std(nonzeros(maxEdchartS(num,1:nonzeros(counterS(1,num,4)-1),4)));
        

    end
    
%     figure(perm);
%     errorbar(nonzeros(vfAve(:,:,i)), nonzeros(effAve(1,:,i))', 3*stdev(counter(1,:,i)-1~=0,i), dottps(i), 'linewidth',1.5);
%     
%     figure(efld)
%     errorbar(nonzeros(vfAve(:,:,i)), nonzeros(aveE(1,:,i))', 3*stdevE(counter(1,:,i)-1~=0,i), linetps(i), 'linewidth',1.5);
% 
%     figure(efld2)
%     errorbar(nonzeros(vfAve(:,:,i)), nonzeros(aveE(1,:,i))', 3*stdevE(counter(1,:,i)-1~=0,i), dottps(i), 'linewidth',1.5);
%     
%     figure(efldBT)
%     errorbar(nonzeros(vfAve(:,:,i)), nonzeros(aveENP(1,:,i))', 3*stdevENP(counter(1,:,i)-1~=0,i), dottps(i), 'linewidth',1.5);
%     
%     figure(eden);
%     errorbar(nonzeros(vfAve(:,:,i)), nonzeros(aveMaxED(1,:,i))', 3*stdevMED(counter(1,:,i)-1~=0,i), dottps(i),'linewidth',1.5);
%     
%     figure(bdfs);
%     errorbar(nonzeros(vfAve(:,:,i)), nonzeros(avBD(1,:,i))', stdevBD(counter(1,:,i)-1~=0,i), dottps(i),'linewidth',1.5);
%     
end
ctemp = nonzeros(counter(1,:,4)-1);
l = length(ctemp)

ctemp = nonzeros(counter(1,:,4)-1);
effAve(1,1:l,4) = nonzeros(sum(sum(effchart(:,:,1:3),3),2))./ctemp;
effAve(1,l+1:vl,4) = 0;
aveE(1,1:l,4) = nonzeros(sum(sum(avEchart(:,:,1:3),2),3))./ctemp;
aveE(1,l+1:vl,4) = 0;
aveENP(1,1:l,4) = nonzeros(sum(sum(avENPchart(:,:,1:3),2),3))./ctemp;
aveENP(1,l+1:vl,4) = 0;
avBD(1,1:l,4) = nonzeros(sum(sum(bdchart(:,:,1:3),2),3))./ctemp;
avBD(1,l+1:vl,4) = 0;
aveMaxED(1,1:l,4) = nonzeros(sum(sum(maxEdchart(:,:,1:3),2),3))./ctemp;
aveMaxED(1,l+1:vl,4) = 0;

vfAve(1,1:l,4)  = nonzeros(vfAve(1,:,4))./ctemp;
vfAve(1,l+1:vl,4) = 0;

counter(1,:,4) = sum(counter(1,:,:)-1,3);

%------------------for loop ended here ------------%
%% EFFECTIVE PERMITTIVITY PLOTS
figure(perm)
c = get(gca,'colororder') % variation of permittivity as scatterplot
scatter(btVFsct(:,1),effSct(:,1),[],c(6,:),'.');
scatter(btVFsct(:,2),effSct(:,2),[],c(6,:),'o');
scatter(btVFsct(:,3),effSct(:,3),[],c(6,:),'x');
errorbar(nonzeros(vfAve(:,:,4)), nonzeros(effAve(1,:,4))', 3*stdev(counter(1,:,4)~=0,4), '--k', 'linewidth',1);
ylabel("Effective Relative Permittivity [F/M]");
xlabel("Volume Fraction of Nanoparticles [%]");
legend('Percolating Geometries','Near Percolating Geometries', 'Non-Percolating Geometries', 'All Geometries', 'Location', 'southeast')

figure(perm2)% mean of permittivity with perry and bruggeman
plot(nonzeros(vfAve(:,:,4)), nonzeros(effAve(1,:,4))', '--k', 'linewidth',1.5);
plot(perrydat(1,:), perrydat(2,:),'ok','linewidth',1.5)
%plot(vfs, brug,'--k','linewidth',1.5)
ylabel("Effective Relative Permittivity [F/M]");
xlabel("Volume Fraction of Nanoparticles [%]");
legend('Mean of tested geometries','Experimental Data (Perry 2009)', 'Bruggeman Model', 'All Geometries', 'Location', 'southeast')


%% ENERGY DENSITY  PLOTS

figure(eden)
errorbar(nonzeros(vfAve(:,:,4)), nonzeros(aveMaxED(1,:,4))', 3*stdevMED(counter(1,:,4)~=0,4), dottps(4), 'linewidth',2);

ylabel("Maximum Energy Density [J/cm^3]");
xlabel("Volume Fraction of Nanoparticles [%]");
legend('Percolating Geometries','Near Percolating Geometries', 'Non-Percolating Geometries', 'All Geometries', 'Location', 'southeast')

figure(eden2)
%plot(nonzeros(vfAve(:,:,4)), nonzeros(aveMaxED(1,:,4))', linetps(4), 'linewidth',1.5);
c = get(gca,'colororder')
scatter(nonzeros(btVFsct(:,1)),nonzeros(maxEdSct(:,1)),[],c(6,:),'.');
scatter(nonzeros(btVFsct(:,2)),nonzeros(maxEdSct(:,2)),[],c(6,:),'o');
scatter(nonzeros(btVFsct(:,3)),nonzeros(maxEdSct(:,3)),[],c(6,:),'x');
plot(vfAve(:,:,4),aveMaxED(:,:,4),'--k','linewidth',1.5) % Plot mean line
%plot(perryed(1,:), perryed(2,:),'.-.','linewidth',1.5)
ylabel("Maximum Energy Density [J/cm^3]");
xlabel("Volume Fraction of Nanoparticles [%]");
legend('Percolating Geometries','Near Percolating Geometries', 'Non-Percolating Geometries', 'Mean', 'Location', 'southeast')
 
%% BREAKDOWN FIELD STRENGTH PLOTS
figure(bdfs)
c = get(gca,'colororder') % variation of permittivity as scatterplot
scatter(nonzeros(btVFsct(:,3)),nonzeros(bdSct(:,3)).*1E-3,[],c(6,:),'x');
scatter(nonzeros(btVFsct(:,2)),nonzeros(bdSct(:,2)).*1E-3,[],c(6,:),'o');
scatter(nonzeros(btVFsct(:,1)),nonzeros(bdSct(:,1)).*1E-3,[],c(6,:),'.');

errorbar(nonzeros(vfAve(:,:,4)), nonzeros(avBD(1,:,4))'.*1E-3, 3*stdevBD(counter(1,:,4)~=0,4).*1E-3, '-k', 'linewidth',1);
%plot(perrybd(1,:), perrybd(2,:),'ok','linewidth',1.5)
ylabel("Breakdown Field Strength [MV/M]");
xlabel("Volume Fraction of Nanoparticles [%]");
legend('Percolating Geometries','Near Percolating Geometries', 'Non-Percolating Geometries', 'All Geometries','Measured Data (Perry)', 'Location', 'northeast')

figure(bdfs2)
plot(nonzeros(vfAve(:,:,4)), nonzeros(avBD(1,:,4))', '--k', 'linewidth',1.5);
plot(perrybd(1,:), perrybd(2,:),'ok','linewidth',1.5)
ylabel("Breakdown Field Strength [MV/M]");
xlabel("Volume Fraction of Nanoparticles[%]");
legend('Simulated Geometries','Measured Data (Perry)', 'Non-Percolating Geometries', 'All Geometries', 'Location', 'northeast')


%% MEAN E-FIELD PLOTS
figure(efld)
scatter(btVFsct(:,1),avESct(:,1)/Em2p,[],c(6,:),'.');
scatter(btVFsct(:,2),avESct(:,2)/Em2p,[],c(6,:),'o');
scatter(btVFsct(:,3),avESct(:,3)/Em2p,[],c(6,:),'x');
scatter(btVFsct(:,1),avENPSct(:,1)/Em2p,[],c(6,:),'.');
scatter(btVFsct(:,2),avENPSct(:,2)/Em2p,[],c(6,:),'o');
scatter(btVFsct(:,3),avENPSct(:,3)/Em2p,[],c(6,:),'x');
errorbar(nonzeros(vfAve(:,:,4)), nonzeros(aveE(1,:,4)/Em2p)', 3*stdevE(counter(1,:,4)~=0,4)/Em2p, '--k', 'linewidth',1);
errorbar(nonzeros(vfAve(:,:,4)), nonzeros(aveENP(1,:,4)/Em2p)', 3*stdevENP(counter(1,:,4)~=0,4)/Em2p, '--k', 'linewidth',1);
ylabel("Average Normalized Electric Field");
xlabel("Volume Fraction of Nanoparticles [%]");
%legend('Inside the Host Material','Inside the Nanoparticles', 'Location', 'southeast')

% figure(efldBT)
% scatter(btVFsct(:,1),avENPSct(:,1),[],c(1,:),'o');
% scatter(btVFsct(:,2),avENPSct(:,2),[],c(2,:),'o');
% scatter(btVFsct(:,3),avENPSct(:,3),[],c(5,:),'o');
% ylabel("Average Electric Field within Nanoparticles");
% xlabel("Volume Fraction of BT Nanoparticles");
% errorbar(nonzeros(vfAve(:,:,4)), nonzeros(aveENP(1,:,4))', stdevENP(counter(1,:,4)~=0,4), linetps(4), 'linewidth',2);

% legend('Percolating Geometries','Near Percolating Geometries', 'Non-Percolating Geometries', 'All Geometries', 'Location', 'southeast')
%% OTHER MISC. PLOTS
% max V vs. permittivity
figure(volts)
%errorbar(nonzeros(vfAve(:,:,4)), nonzeros(avBD(1,:,4))', 3*stdevBD(counter(1,:,4)~=0,4), linetps(4), 'linewidth',1);
c = get(gca,'colororder') % variation of permittivity as scatterplot
scatter(nonzeros(btVFsct(:,1)),nonzeros(MaxVolt2p(:,1)),[],c(6,:),'o');
scatter(nonzeros(btVFsct(:,2)),nonzeros(MaxVolt2p(:,2)),[],c(6,:),'o');
scatter(nonzeros(btVFsct(:,3)),nonzeros(MaxVolt2p(:,3)),[],c(6,:),'o');

scatter(nonzeros(btVFsct(:,1)),nonzeros(maxVolt2pNP(:,1)),[],c(6,:),'.');
scatter(nonzeros(btVFsct(:,2)),nonzeros(maxVolt2pNP(:,2)),[],c(6,:),'.');
scatter(nonzeros(btVFsct(:,3)),nonzeros(maxVolt2pNP(:,3)),[],c(6,:),'.');
ylabel("Maximum Voltage [V]");
xlabel("Volume Fraction of Nanoparticles [%]");
legend('Percolating in Host','Near Percolating in Host', 'Non-Percolating in Host', 'Percolating in NP','Near Percolating in NP', 'Non-Percolating in NP', 'Location', 'southeast')

figure(volts2)
%errorbar(nonzeros(vfAve(:,:,4)), nonzeros(avBD(1,:,4))', 3*stdevBD(counter(1,:,4)~=0,4), linetps(4), 'linewidth',1);
c = get(gca,'colororder') % variation of permittivity as scatterplot
scatter(nonzeros(effSct(:,1)),nonzeros(MaxVolt2p(:,1)),[],c(6,:),'o');
scatter(nonzeros(effSct(:,2)),nonzeros(MaxVolt2p(:,2)),[],c(6,:),'o');
scatter(nonzeros(effSct(:,3)),nonzeros(MaxVolt2p(:,3)),[],c(6,:),'o');

scatter(nonzeros(effSct(:,1)),nonzeros(maxVolt2pNP(:,1)),[],c(6,:),'.');
scatter(nonzeros(effSct(:,2)),nonzeros(maxVolt2pNP(:,2)),[],c(6,:),'.');
scatter(nonzeros(effSct(:,3)),nonzeros(maxVolt2pNP(:,3)),[],c(6,:),'.');
ylabel("Maximum Voltage [V]");
xlabel("Effective Permittivity [F/M]");
legend('Percolating in Host','Near Percolating in Host', 'Non-Percolating in Host','Percolating in NP','Near Percolating in NP', 'Non-Percolating in NP', 'Location', 'northeast')


%% COUNTER BAR PLOT
figure(countbars)
allcounts = [counter(:,:,1)'-1, counter(:,:,2)'-1, counter(:,:,3)'-1];% counter(:,:,4)'];
%hb = bar(vfAve(:,:,4)', allcounts, 'stacked');
ylabel("Number of Microstructures");
xlabel("Volume Fraction of Nanoparticles [%]");
ylim([0,80]);
legend('Percolaging Geometries', 'Near Percolating Geometries', 'Non-Percolating Geometries')
saveas(countbars,'bars2','fig');
saveas(countbars,'bars2','jpeg');
saveas(perm,'permValid2','fig');
saveas(perm,'permValid2','jpeg');
saveas(eden,'eden2','fig');
saveas(eden,'eden2','jpeg');
saveas(efldBT,'efldBT2','fig');
saveas(efldBT,'efldBT2','jpeg');
saveas(bdfs,'bdfs2','fig');
saveas(bdfs,'bdfs2','jpeg');
saveas(efld,'efld_hst2','fig');
saveas(efld,'efld_hst2','jpeg');
saveas(perm2,'perm3','fig');
saveas(perm2,'perm3','jpeg');
saveas(eden2,'eden3','fig');
saveas(eden2,'eden3','jpeg');
saveas(efldBT2,'efldBT3','fig');
saveas(efldBT2,'efldBT3','jpeg');
saveas(bdfs2,'bdfs3','fig');
saveas(bdfs2,'bdfs3','jpeg');
saveas(efld2,'efld_hst3','fig');
saveas(efld2,'efld_hst3','jpeg');
saveas(volts,'volts1','fig');
saveas(volts,'volts1','jpeg');
saveas(volts2,'volts2','fig');
saveas(volts2,'volts2','jpeg');

figure(MCP)
hold on
%MCPdat = [3.9,3.9,6.3,6.7,7.1,7.1,8.5,8.5,21.0,21.0,21.0,19.0,19.0,23.4,26,26,26.8,29,31,60,95,95,95,95,50,98;1400,1310,630,870,760,920,620,720,470,400,300,250,600,350,320,420,420,380,400,320,300,250,200,100,330]
Mcpdat = [3.9,1400;3.9,1310;6.3,630;6.7,870;7.1,760;7.1,920;8.5,620;8.5,720;21.0,470;21.0,400;21.0,300;19.0,250;19.0,600;23.4,350;26,320;26,420;26.8,420;29,380;31,400;60,320;95,300;95,250;95,200;95,100;50,300;98,330;]
scatter(effSct(:,3),bdSct(:,3).^2,[],c(6,:),'x');
scatter(effSct(:,2),bdSct(:,2).^2,[],c(6,:),'o');
scatter(effSct(:,1),bdSct(:,1).^2,[],c(6,:),'.');
c = get(gca,'colororder') % variation of permittivity as scatterplot
saveas(MCP,'MCP','fig');


%plot(Mcpdat(:,1), Mcpdat(:,2).^2, 'ok', 'linewidth',1);
ylabel("Breakdown Field Strength^2 [MV/M]");
xlabel("Effective Relative Permittivity [F/M]");
legend('Percolating Geometries','Near Percolating Geometries', 'Non-Percolating Geometries', 'Data from ', 'Location', 'southeast')
set(gca,'xscale','log')
set(gca,'yscale','log')

% Get a list of all variables
allvars = whos;
% Identify the variables that ARE NOT graphics handles. This uses a regular
% expression on the class of each variable to check if it's a graphics object
tosave = cellfun(@isempty, regexp({allvars.class}, '^matlab\.(ui|graphics)\.'));
% Pass these variable names to save
save('dataset2_nf.mat', allvars(tosave).name)
