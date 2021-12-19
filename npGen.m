function [coords] = npGen(position, dimension, radii, shells, seed)
%   Generates single spherical nanoparticle of N shells in a given 3d
%             ndSparse matrix.
%
% Position:  vector of length 3 -> Center of nanoparticle = [x,y,z]
% Dimension: vector of length 3 -> Matrix boundaries = [Xmax,Ymax,Zmax]
% Radii:     vector of length N -> Diam. of shells = [0;coreS; ... ; NthS]
% Shells:    vector of length N -> dielectric const = [0;coreDC; ... ; NthDC]
%
%
% new:      Omega
% coords:   List of voxel coordinates and material value at each.
rng(seed+1);

if ((length(position) == 3) && (length(dimension) == 4))
    numSh = length(shells);
    
    % Populate Core Centers
    coreMat = zeros([1,length(position)]) + shells(2);
    
    if (sum(position <= 0)  + sum(position >= dimension(3))) == 0
        coords = [position,coreMat(1)];
    else
        coords = zeros([0,4]);
    end
    
    shi = 1; % Indicates current shell being constructed [SHell Index]
    
    % Define counters to corner of nanoparticle construction box
    % 2R x 2R x 2R;  X,Y,Z Indexing vectors for iteration
    xi = 1;
    yi = 1;
    zi = 1;
    ci = 1;
    
    err = 0;
    % initialize coordinates list
    
    % loop through each material type at x,y,z, plotting concentric spheres
    % to form single nanoparticle as defined
    for(shi = length(radii):-1:2)
        ci = shi;
        xblow = position(1)-round(radii(shi));
        xbhi  = position(1)+round(radii(shi));
        if(xblow < 1)
            xblow = 1;
        end
        if(xbhi > dimension(1))
            xbhi = dimension(1);
        end
        for(xi = xblow:xbhi)
            %disp((xi-xblow)/(xbhi-xblow));
            yblow = position(2)-round(radii(shi));
            ybhi  = position(2)+round(radii(shi));
            if(yblow < 1)
                yblow = 1;
            end
            if(ybhi > dimension(2))
                ybhi = dimension(2);
            end
            for(yi = yblow:ybhi)
                zblow = position(3)-round(radii(shi));
                zbhi  = position(3)+round(radii(shi));
                if(zblow < 1)
                    zblow = 1; 
                end
                if(zbhi > dimension(3))
                    zbhi = dimension(3);
                end
                for(zi = zblow:zbhi)
                    
                    XYZ = [int16(xi),int16(yi),int16(zi)];
                    
                    %voxel = accumarray(XYZ, shells(ci), dimension(1:3));
                    
                        
                        dist = ((xi-position(1))^2 + (yi-position(2))^2 + (zi-position(3))^2)^(0.5);
                        if(((dist <= radii(ci)) && (radii(ci-1) < dist)))
                            %np = voxel;  % logical OR with nanoparticle sparse matrix to add to space
                            coords = [coords;xi,yi,zi,shells(shi)]; % Save coordinates for later
                        end
                end
            end
        end
    end
else
    disp('error: No Nanoparticle generated, check nanoparticle ')
end

end

