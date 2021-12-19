function [SmoothEfld] = CorrectField(Efld, Omega, filterSize, filterType)
%UNTITLED2 Summary of this function goes here
% Example Usage:  SmoothedField = CorrectField(RawField, Omega, [5,5,5])
% This code takes a 'jagged' or raw electric field from high discrepancy in 
% permittivity or conductivity and smooths the field to more accurately 
% capture the continuity. 
% In essence, this adapts the shape of the electric field to better align 
% with the voxellated nanoparticles it contours

% This is achieved with a conditional / weighted 27-point 3D averaging filter 



SmoothEfld = Efld;

aveBox = zeros(filterSize);
% get 
boxRad = (filterSize(1)-1)/2;

[X,Y,Z] = size(Efld);

% Notes:  1. Make circular about X and Y
%         2. Make sure to do a 5x5x3 and 5x5x4 at top and bottom two layers
%         3. divide by number of valid voxels in averaging cluster 

% iterate through each voxel
for i = 1:X
    for j = 1:Y
        for k = 1:Z
            % at each voxel, check to see how many of the 27 adjacent
            % voxels are valid (And which ones that is...)
            %% 0. reset averaging box
            aveBox = zeros(filterSize);
            nnz = 0;
            %% 0.5. check if voxel is within np or not...
            if Omega(i,j,k) == 0
                %% 1. Populate AveBox (leave zero where invalid, i.e. inside NP)
                for ii = -boxRad:boxRad
                    for jj = -boxRad:boxRad
                        for kk = -boxRad:boxRad
                            % set indices for Average box
                            iib = ii + boxRad + 1;
                            jjb = jj + boxRad + 1;
                            kkb = kk + boxRad + 1;
                            % Set indices for E-field
                            ie = circ(i + ii, X);
                            je = circ(j + jj, Y);
                            ke = circ(k + kk, Z);
                            % Madke sure z-index is within bounds...
                            % AND location is inside host (not NP...)
                            if Omega(ie,je,ke) == 0  && (ke > 0 || ke <= Z)
                                aveBox(iib,jjb,kkb) = Efld(ie,je,ke);
                                nnz = nnz + 1;
                            else
                                aveBox(iib,jjb,kkb) = 0;
                            end
                        end
                    end
                end
                %% 2. Compute average (based only on valid voxels)
                if filterType == 'G'  % Geometric Average
                    aveBox = aveBox *1E-7;
                    SmoothEfld(i,j,k) = prod(real(nonzeros(aveBox())))^(1/nnz)*1E7;
                elseif filterType == 'A'  % Arithmetic Average
                    SmoothEfld(i,j,k) = sum(sum(sum(aveBox)))/nnz;
                end
                
                %disp("[" + i +", "+ j +", "+ k + "]")
            end
           
        end
    end
end

