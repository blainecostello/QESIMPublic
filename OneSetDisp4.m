%% Display Single Set
% TO DO:  
%   1. Eliminate duplicates
%       (Use Errorbar instead of points) - consolidate VF clusters
%       OR
%       (induce randomness with nonzero standard deviation)
%       [TEMP SOL: get average of points at equal VF]
%   2. 
% Set framerate...
fr = 10;

maxE2p = max(max(max(E2p)))
minE2p = min(min(min(E2p)))


gd = GD(1)

Xsz = gd
Ysz = gd 
Zsz = gd-1

vidnm1 = char("Efld_"+Seed+"_"+round(VF*1000)+"_"+runID+".avi")
vidnm2 = char("Geom_"+Seed+"_"+round(VF*1000)+"_"+runID+".avi")


X1 = VideoWriter(vidnm1);
X1.FrameRate = fr;
open(X1);

X2 = VideoWriter(vidnm2);
X2.FrameRate = fr;
open(X2);

for n = 1:length(bdp)
E2p2(bdp(n,1),bdp(n,2),bdp(n,3)) = (1.75*E2p(bdp(n,1),bdp(n,2),bdp(n,3)));
end

    x = [0, gd+1, gd+1, 0]
    z = [0, 0, gd+1, gd+1]
parfor k = 1:gd
% 2 phase
    figure
    imagesc(squeeze(rot90(abs(E2p(k,:,:)))), [real(minE2p),real(maxE2p)])
    e2p(k) = getframe(gcf);
    writeVideo(X1, e2p(k).cdata);
    

% point in volume
    y = [k, k, k, k]
    hold on   
    plotSparse(Shells, Voxels_bt, GD(1), GD(2), GD(3), 'fast');
    view([40,150])

    patch(y,x,z,'blue','FaceAlpha',0.4);
    
    scr(k) = getframe(gcf);
    writeVideo(X2, scr(k).cdata);

    
close
close

end
% 
% e2p = [e2p;e2p];
% e3p = [e3p;e3p];
% e3pH = [e3pH;e3pH];
% e3pS = [e3pS;e3pS];
% writeVideo(X1, e2p.cdata);
% writeVideo(X2, e3p.cdata);
% writeVideo(X3, e3pH.cdata);
% writeVideo(X4, e3pS.cdata);


close(X1)
close(X2)


plotSparse(Shells, Voxels_bt, Xsz, Ysz, Zsz, 'fast');
view([40,150])

%Voxels_one = npMatPop([22,22,22], [11,11,11,10], Radii, Shells, StDev, runID,Seed);

%plotSparse(Shells, Voxels_one, 22, 22, 22, 'fast');

% for i = 1:180
%     npPos(:,4) = 2;
%     plotSparse(Shells, npPos, gd,gd,gd,'fast');
%     view([2*i,150])
%     scr(i) = getframe(gcf);
%     writeVideo(X3, scr(i).cdata);
%     disp(i)
%     close
% 
%     plotSparse(Shells, Voxels_bt, gd,gd,gd,'fast');
%     view([2*i,150])
%     scr(i) = getframe(gcf);
%     writeVideo(X4, scr(i).cdata);
%     disp(i)
%     close
%     
% end

