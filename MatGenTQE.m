function MatGenTQE(GD, Radii, Shells, Ovlp, VF, runID, Seed)
%PaceParalleltoolbox_r2016b
%% MatGenT1 - Material Generator [Test1]
%  This code generates a 3D finite difference grid of dielectric constants
%  to set up the inputs for the electromagnetic analysis of these simulated
%  geometries.
%
% ARGS:
%   GD should be a the grid dimensions in the form of [X Y Z] limit values
%   Radii is the radius of each shell, the first entry is the core and
%       should always be zero.
%   Shells should be [0 1 2 ..] to the number of shells.  The order can be
%       changed to switch the order of shelling
%   Ovlp is an integer representing the allowed overlap between nanoparticles
%   VF is a double between 0 and 0.6 representing the volume fraction of
%       nanoparticle inclusions.
%   RunID is the experiment ID and defines the subfolder of data storage
%   Seed is the random seed used to generate the geometry. unique
%       combinations of seed, vf, and GD will produce unique geometries that
%       can be easily recreated with varying parameters.
%------===-----==----------------===--------------===------------==-------=

% Construct easily identifiable, unique filename
dir_nm = char("Run_" + runID);
mkdir(dir_nm)
fnm = "" + dir_nm + "/" + GD(1) + "_" + round(VF*1000) + "E-3_" + Seed + ".mat";



if exist(fnm, 'file')
    disp(fnm)
    disp("File exists, skipping to next");
else
    
    % Declare necessary variables
    maxRad = Radii(length(Radii)-1);
    coreRad = Radii(2);
    stdev = 0;
    Omega = zeros(GD);
    btVox = 0;
    %%=========================================%%
    % USE OLD TESTER FILES TO FILL IN THIS PART %
    
    % 1. Generate nanoparticle position list
    
%     GD
%     
%     VF
%     
%     Ovlp
%     
%     maxRad
%     
%     stdev
%     
%     Seed
    
    npPos = npMatGen(GD, VF, Ovlp, coreRad, stdev, 'iterateNP',Seed);
    npPos(:,4) = ones([length(npPos(:,4)),1]).*Radii(length(Radii));
    % 2. Populate nanoparticle voxel list
    Voxels_bt = npMatPop(GD, npPos, Radii, Shells, [stdev;stdev], runID, Seed);
    Voxels_bt = unique(Voxels_bt, 'rows');
    % 3. Write voxel list to 3D matrix
    btVox = 0;
    shVox = 0;
    for i = 1:length(Voxels_bt)
                
       % Omega formation is updated to prioritize core over shells, provides
       % physical representation of interphase formation around clusters of
       % NP material.
        if(Omega(round(Voxels_bt(i,1)),round(Voxels_bt(i,2)),round(Voxels_bt(i,3))) == 0 || Omega(round(Voxels_bt(i,1)),round(Voxels_bt(i,2)),round(Voxels_bt(i,3))) > Voxels_bt(i,4))
            Omega(round(Voxels_bt(i,1)),round(Voxels_bt(i,2)),round(Voxels_bt(i,3))) = Voxels_bt(i,4);
            if(Voxels_bt(i,4) == 1)
                btVox = btVox + 1;
            elseif(Voxels_bt(i,4) == 2)
                shVox = shVox + 1;
            end
        end
    end
    btVF = btVox/(GD(1)*GD(2)*GD(3));
    shVF = shVox/(GD(1)*GD(2)*GD(3));
    btVFsh = (btVox + (shVox/2))/(GD(1)*GD(2)*GD(3)) ;
    if(length(Radii) == 3)
        ro = Radii(3);
        ri = Radii(2);
        btVFshHR = (btVox + shVox*(((1/2*(ro+ri))^3 - ri^3)/(ro^3-ri^3)))/(GD(1)*GD(2)*GD(3))
    end
    % 4. Determine percolating configuration (prox; 0=perc, 1=nearperc, 2=nonperc)
    [prox,percPath] = percAnalyzer(npPos, GD(3));
    
    %          plotSparse([0,1], Voxels_bt, GD(1), GD(2), GD(3), 'fast');
    %          title("Nanoparticle Geometry")
    %          view([40,150]);
    
    %%=========================================%%
    % (5)
    % Save generated geometry in specified filename 'fnm'
    save(fnm) % Add list of specific variables to save (memory efficiency)
    pause(0.1)
    disp("File has been saved: ./" + fnm)
    
  plotSparse(Shells, Voxels_bt, GD(1), GD(2), GD(3), 'fast');
  view([20,10])  
  set(gcf, "Position",[100,100,800,500])
  pause(0.5)
  print(gcf, 'geom1', '-djpeg', '-r600');
    
end

end



