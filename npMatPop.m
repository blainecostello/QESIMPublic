function [Voxels] = npMatPop(gridDims,npMat,radii,materials,sigma,geoID,seed)
% NP Matrix Populator
%   generates 3D matrix of specified grid dimensions and [Nx4] vector of 
%   3-space coordinates representing voxels to make up nanoparticle spheres 
%   Outputs a file including voxel coordinates and centerpoint 
%
%==================={ ARGUMENTS }=====================%
%   gridDims  = dimensions of simulated 3-space grid      [X, Y, Z]
%   npMat     = matrix of np centerpoints and max radius  [4 x N]
%   shells    = mean maximum outer-shell radius           [units: voxels]
%   sigma     = standard deviation of outer-shell radius  [units: voxels]
%   geoID .   = geometry ID (used in file generation)     [units: integer]
%   
%==================={ OUTPUTS }=====================%
%   Omega     = 3D matrix representing material value at each index 
%   Voxels    = Matrix of nanoparticle coordinates and mat # at each [Nx4]
%
%
%
i = 0;
sorted = 0;
rng(seed+3);
Voxels = zeros([0,4]);

nnp = length(npMat(:,1));
nrad = length(radii(:,1));
% For each nanoparticle defined in NP Mat... if any nanoparticles at all
if(nnp >= 1)
for(i=1:nnp)
    %disp(i);
    XYZ = [npMat(i,1),npMat(i,2),npMat(i,3)];
    sorted = 0;
    while(sorted == 0)
        radiiAdj = normrnd(radii(1:nrad-1), sigma(1:(length(sigma)-1))');
        if(~issorted(radiiAdj))
            sort(radiiAdj);
        end
        if(radiiAdj(nrad-1) < npMat(i,4))
            radiiAdj = [radiiAdj;npMat(i,4)];
        end
        sorted = 1;
    end
    % adjust shell radii
    
    [voxNP] = npGen(XYZ, [gridDims,8], radiiAdj, materials,seed);
    Voxels = [Voxels; voxNP];
    %disp("NP Generation Time:")
    
%     x = npMat(i,1);
%     y = npMat(i,2);
%     z = npMat(i,3);
%     r = ceil(radiiAdj());

    %Omega(x-r:x+r, y-r:y+r, z-r:z+r)  = Omega(x-r:x+r, y-r:y+r, z-r:z+r) + omNP;
end
end


end

