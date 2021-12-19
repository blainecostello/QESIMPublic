% Show geometry 

% prompt user to input

k = round(GD(3)/2);
%k = input("Specify the Y-position for the cross-section you wish to examine...\n")


% Set up some stuff...
maxE2pR = max(max(max(real(E2p))))
minE2pR = min(min(min(real(E2p))))

maxE2pI = max(max(max(imag(E2p))))
minE2pI = min(min(min(imag(E2p))))

gd = GD(1);
x = [0, gd+1, gd+1, 0];
z = [0, 0, gd+1, gd+1];

% 2 phase (real)

figure

imagesc(squeeze(rot90(imag(E2p(k,:,:)))), [(minE2pI),(maxE2pI)])
set(gcf, "Position", [600,400,500,400])
colorbar

figure

imagesc(squeeze(rot90(real(E2p(k,:,:)))), [(minE2pR),(maxE2pR)])
set(gcf, "Position", [600,100,500,400])
colorbar


% point in volume
y = [k, k, k, k];
plotSparse(Shells, Voxels_bt, GD(1), GD(2), GD(3), 'fast');
view([40,150])

patch(y,x,z,'blue','FaceAlpha',0.4);
set(gcf, "Position", [100,100,450,400])
