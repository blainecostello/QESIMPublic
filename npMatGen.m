function [npCoords] = npMatGen(gridDims, vf, OL, maxRad, sigma, algType, seed)
% npMatGen (Nanoparticle Matrix Generator)
%   generates [4 by N] vector of 3-space coordinates representing
%   centerpoint of each nanoparticle in grid along with each nanoparticle's
%   maximum radius.  A series of algorithms are supported, and a PDF for
%   each will be determined using a small scale Monte Carlo simulation.
%
%==================={ ARGUMENTS }=====================%
%   gridDims  = dimensions of simulated 3-space grid      [X, Y, Z]
%   vf        = volume fraction                           [0 < vf < 0.45]
%   maxRad    = mean maximum outer-shell radius           [units: voxels]
%   sigma     = standard deviation of outer-shell radius  [units: voxels]
%   algType   = 'iterateNP' || 'brutefrce' || 'selection' || 'annealing' || ...
%
%==================={ OUTPUTS }=====================%
%   npCoords  = matrix of nanoparticle coordinates and radius at each [Nx4]
%   actVF     = Actual volume fraction computed after generating all NP
%
%=============={ ABOUT THE ALGORITHMS }===============%
%  'iterateNP'   = Iterative Algorithm
%          This algorithm is the initial concept for the randomization of
%          nanoparticle positions and radii, but is hypothesized to produce
%          a non-uniform random distribution of coordinates along the
%          simulated space.
%
%  'brutefrce'   = Brute Force Algorithm
%          This algorithm is the simplest to code, but can result in a null
%          return case or take an unnecessarily long time if volume
%          fraction exceeds a certain value.
%
%  'selection'  = Selective Re-randomization Algorithm
%          This algorithm is similar to the brute force algorithm, but
%          in the event of an invalid geometry, only the intersecting
%          particles are re-computed instead of recomuting the entire
%          nanoparticle field.  A counter will be implemented to avoid
%          infinite looping through impossible geometries.
%
%  'annealing'  = Annealing Algorithm
%          Similar to the thermal process of metallurgical Annealing,
%          this algorithm begins with a randomized set of points and loops
%          through a 'cooling' process.  During the cooling process, each
%          particle will gradually 'settle' into a position such that a
%          valid, non-intersecting geometry results.
%
%
%% Initialize
% Set Flags
valid = 0;
rng(seed+2);
% Initialize Counters
i = 0;   % nanoparticle generator loop counter
j = 0;   % loop - safety counter (avoid infinite looping)
c = 0;   % cooling loop counter
o = 0;
% Overlap (OL) of nanoparticles 

% Estimate number of nanoparticles (round up for sliced particles)

eVolNP = (4/3 * pi * (maxRad).^3);
simVol = gridDims(1) * gridDims(2) * gridDims(3);
numNP = floor((vf * simVol)./eVolNP)

% Initialize Vectors
npCoords = zeros(numNP,4);
npCirc = zeros(0,5);

%% Check if valid geometry is possible & Run Selected Alg
if(vf > 0.9)
    
    disp('Volume Fraction too high:');
    disp('    (this will result in an invalid geometry or very long computation time)')
elseif(vf == 0)
elseif(algType == 'twoNPssss')
    disp("Generating Two Nanoparticles")
    pause(5)
    npPos = [round(gridDims(1)/2), round(gridDims(2)/2), round(gridDims(3)/2)+round(seed/2), maxRad];
    npPos = [npPos; round(gridDims(1)/2), round(gridDims(2)/2), round(gridDims(3)/2)-round(seed/2), maxRad];
    npCoords = npPos;
elseif(algType == 'singleNPE')
    disp("Generating One Nanoparticle")
    pause(5)
    npPos = [round(gridDims(1)/2), round(gridDims(2)/2), round(gridDims(2)/2), maxRad];
    npCoords = npPos;
elseif(algType == 'iterateNP')
    
    zbuffer = 3
    %% Iterative:   1. The first Random x,y,z coordinate is generated along
    %                  with gaussian normal max-radius value.
    %               2. Subseqent random x,y,z coordinate and max-radius is
    %                  computed and checked against existing points.  if
    %                  invalid geometry results, recalculate and overwrite
    %                  the point.
    %               3. Repeat until specified volume fraction is reached.
    
    % 1.0 - Generate initial XYZ coordinate and radius
    % X and Y axes will be circular, and thus only need to cover
    % centerpoints inside of simulated space.
        npPosX = round((rand() .* gridDims(1)));
        npPosY = round((rand() .* gridDims(2)));
    % Z axis will be 
        % npPosZ = round((rand() .* (gridDims(3) + 2*maxRad)) - maxRad);
        npPosZ = round((rand()*(gridDims(3) - 2*maxRad - 2*zbuffer)) + maxRad + zbuffer)
    % Compute adjusted radius on normal distribution
        npRad = normrnd(maxRad, sigma);
        npMat = [npPosX, npPosY, npPosZ, npRad];
        
        % Check for edge overlap (var 'o' is used as integer case flag)
        if(npPosX > (gridDims(1)-npRad))
            npCirc = [npCirc; (npPosX-gridDims(1)),npPosY,npPosZ,npRad,1];
            o = 1;
        elseif(npPosX <= npRad)
            npCirc = [npCirc; (npPosX+gridDims(1)),npPosY,npPosZ,npRad,1];
            o = 2;
        else
            o = 0;
        end
        % Y axis handling (detect X & Y overlap and populate accordingly)
        if(npPosY > (gridDims(2)-npRad))
            npCirc = [npCirc; npPosX,(npPosY-gridDims(2)),npPosZ,npRad,1];
            % handle multi-plane overlap
            if(o == 1)
                npCirc = [npCirc; (npPosX-gridDims(1)),(npPosY-gridDims(2)),npPosZ,npRad,1];
            elseif(o == 2)
                npCirc = [npCirc; (npPosX+gridDims(1)),(npPosY-gridDims(2)),npPosZ,npRad,1];
            end
        elseif(npPosY <= npRad)
            npCirc = [npCirc; npPosX,(npPosY+gridDims(2)),npPosZ,npRad,1];
            % handle multi-plane overlap
            if(o == 1)
                npCirc = [npCirc; (npPosX-gridDims(1)),(npPosY+gridDims(2)),npPosZ,npRad,1];
            elseif(o == 2)
                npCirc = [npCirc; (npPosX+gridDims(1)),(npPosY+gridDims(2)),npPosZ,npRad,1];
            end
        end
                
        % 2.0 - Generate remaining points
        for i = 1:numNP-1
            cnt = 0;
            lmt = round((gridDims(1)*gridDims(2)*gridDims(3))/(3)); % one third the size of the total space
            valid = 0;
            while (valid == 0 && cnt < lmt)
                cnt = cnt + 1;
                %   2.1 - Generate one set of coordinates
                % centerpoints inside of simulated space.
                npPosX = round((rand() .* gridDims(1)));
                npPosY = round((rand() .* gridDims(2)));
                % Z axis will be
                npPosZ = round((rand()*(gridDims(3) - 2*maxRad - 2*zbuffer)) + maxRad + zbuffer)
                npRad = normrnd(maxRad, sigma);
                orig = [npPosX,npPosY,npPosZ,npRad];
                %       2.2 - Check against existing points
                dists = sqrt((npPosX - npMat(:,1)).^2 + (npPosY - npMat(:,2)).^2 + (npPosZ - npMat(:,3)).^2);
                minDists = npMat(:,4) + npRad - (OL);
                if((minDists < dists))
                    valid = 1;
                    copies = 0;
                    % Check for edge overlap (var 'o' is used as integer case flag)
                    if(npPosX > (gridDims(1)-npRad))
                        npCirc = [npCirc; (npPosX-gridDims(1)),npPosY,npPosZ,npRad,i+1];
                        o = 1;
                        copies = copies + 1;
                    elseif(npPosX < npRad)
                        npCirc = [npCirc; (npPosX+gridDims(1)),npPosY,npPosZ,npRad,i+1];
                        o = 2;
                        copies = copies + 1;
                    else
                        o = 0;
                    end
                    % Y axis handling (detect X & Y overlap and populate accordingly)
                    if(npPosY > (gridDims(2)-npRad))
                        npCirc = [npCirc; npPosX,(npPosY-gridDims(2)),npPosZ,npRad,i+1];
                        copies = copies + 1;
                        % handle multi-plane overlap
                        if(o == 1)
                            npCirc = [npCirc; (npPosX-gridDims(1)),(npPosY-gridDims(2)),npPosZ,npRad,i+1];
                            copies = copies + 1;
                        elseif(o == 2)
                            npCirc = [npCirc; (npPosX+gridDims(1)),(npPosY-gridDims(2)),npPosZ,npRad,i+1];
                            copies = copies + 1;
                        end
                    elseif(npPosY < npRad)
                        npCirc = [npCirc; npPosX,(npPosY+gridDims(2)),npPosZ,npRad,i+1];
                        copies = copies + 1;
                        % handle multi-plane overlap
                        if(o == 1)
                            npCirc = [npCirc; (npPosX-gridDims(1)),(npPosY+gridDims(2)),npPosZ,npRad,i+1];
                            copies = copies + 1;
                        elseif(o == 2)
                            npCirc = [npCirc; (npPosX+gridDims(1)),(npPosY+gridDims(2)),npPosZ,npRad,i+1];
                            copies = copies + 1;
                        end
                    end
                % Check distances
                    npCheck = [npMat;npCirc(:,1:4)];
                    dists = sqrt((npPosX - npCheck(:,1)).^2 + (npPosY - npCheck(:,2)).^2 + (npPosZ - npCheck(:,3)).^2);
                    minDists = npCheck(:,4) + npRad;                    
                    % append to end of vector
                    if(dists < minDists)
                        valid = 0;
                        len = length(npCirc);
                        npCirc = npCirc(1:len-copies,:);
                        %disp("circular intersection...")
                    end

                end
            end
            npMat = [npMat; orig];
        end
        
        
    % 3.0 - Return generated coordinate matrix OR error messages...
    npCoords = [npMat;npCirc(:,1:4)];
 
elseif(algType == 'brutefrce') 
    %% Brute force: 1. Random x,y,z coordinates are generated with gaussian
    %                  normal max-radii values for each.
    %               2. Distances between centerpoints are checked against
    %                  radii to determine if geometry is valid.
    %               3. Once valid geometry is randomly generated,
    %                  nanoparticle vector is returned.
    while (valid == 0)
        % 1 - Generate XYZ Coordinates and radii
        % 1.1 - Compute XYZ coordinates.
        npPosX = round((rand(numNP, 1) .* (gridDims(1) + 2*maxRad)) - maxRad);
        npPosY = round((rand(numNP, 1) .* (gridDims(2) + 2*maxRad)) - maxRad);
        npPosZ = round((rand(numNP, 1) .* (gridDims(3) + 2*maxRad)) - maxRad);
        
        % 1.2 - Compute gaussian normal max-radius for each set of coords.
        npRad = normrnd(maxRad, sigma, [numNP, 1]);
        npMat = [npPosX, npPosY, npPosZ, npRad];
        
        % 2.0 - Check for any invalid geometry
        % (assume valid until invalid geometry is detected)
        valid = 1;
        i = 1;
        while ((i < numNP) && (valid == 1))
            dists = sqrt((npPosX - npPosX(i)).^2 + (npPosY - npPosY(i)).^2 + (npPosZ - npPosZ(i)).^2);
            minDists = npRad + npRad(i);
            
            if((minDists >= dists))
                % 2.1 - if invalid geometry is detected, recompute XYZ coordinates
                valid = 0;
            end
            i = i + 1;
        end
        % 2.2 - Repeat while invalid geometry is present (WHILE LOOP)
    end
    % 3.0 - Return valid geometry or time-out error message..
    npCoords = npMat;
    
    
elseif(algType == 'selection')
    %% Selective Re-rand: 1. Random x,y,z coordinate are generated with
    %                        gaussian normal max-radii values for each.
    %                     2. Distances between centerpoints are checked
    %                        against radii to determine if geometry is
    %                        valid.
    %                     3. Once valid geometry is randomly reached,
    %                        nanoparticle vector is returned.
    
    % 1.0 - Generate XYZ Coordinates and radii
        npPosX = round((rand(numNP, 1) .* (gridDims(1) + 2*maxRad)) - maxRad);
        npPosY = round((rand(numNP, 1) .* (gridDims(2) + 2*maxRad)) - maxRad);
        npPosZ = round((rand(numNP, 1) .* (gridDims(3) + 2*maxRad)) - maxRad);
    
        % 1.1 - Compute XYZ coordinates.
        npRad = normrnd(maxRad, sigma, [numNP, 1]);
        % 1.2 - Compute gaussian normal max-radius for each set of coords.
        npMat = [npPosX, npPosY, npPosZ, npRad];
    % 2.0 - Check distances vector for each particle against min-dist vect.
    
    for(i = 1:numNP)
        dists = sqrt((npPosX - npPosX(i)).^2 + (npPosY - npPosY(i)).^2 + (npPosZ - npPosZ(i)).^2);
        minDists = npRad + npRad(i);
        minDists(i) = 0;
        
        % 2.1 - iterate through invalid particles, recompute XYZ, overwrite
        if((minDists >= dists))
            % 2.2 - while invalid geometry is detected, recompute XYZ coordinates
            valid = 0;
            while ((valid == 0))
                npPosX(i) = round((rand() .* (gridDims(1) + 2*maxRad)) - maxRad);
                npPosY(i) = round((rand() .* (gridDims(2) + 2*maxRad)) - maxRad);
                npPosZ(i) = round((rand() .* (gridDims(3) + 2*maxRad)) - maxRad);
                
                npRad(i) = normrnd(maxRad, sigma, [numNP, 1]);
                
                dists = sqrt((npPosX - npPosX(i)).^2 + (npPosY - npPosY(i)).^2 + (npPosZ - npPosZ(i)).^2);
                minDists = npRad + npRad(i);
                minDists(i) = 0;

                
                if((minDists <= dists))
                    % 2.1 - if invalid geometry is detected, recompute XYZ coordinates
                    npMat = [npPosX, npPosY, npPosZ, npRad];
                    valid = 1;
                    
                end
            end
        end
    end
    
    % 3.0 - Check for valid geometry and return OR display error message.
    npCoords = npMat;
    
elseif(algType == 'annealing')
    %% Annealing: 1. Random x,y,z coordinates are generated with gaussian
    %                normal max-radii values for each.
    %             2. Particles are moved away from each other at magnitude
    %                based on distance to nearest particle
    %             3. Repeat, decreasing magnitude of motion as particle
    %                matrix 'cools'.
    %             4. Once valid geometry is reached, nanoparticle vector is
    %                returned.
    
    % 1.0 - Generate XYZ Coordinates and radii
      % 1.1 - Compute XYZ coordinates.
        npPosX = round((rand(numNP, 1) .* (gridDims(1) + 2*maxRad)) - maxRad);
        npPosY = round((rand(numNP, 1) .* (gridDims(2) + 2*maxRad)) - maxRad);
        npPosZ = round((rand(numNP, 1) .* (gridDims(3) + 2*maxRad)) - maxRad);
        
      % 1.2 - Compute gaussian normal max-radius for each set of coords.
        npRad = normrnd(maxRad, sigma, [numNP, 1]);
        npMat = [npPosX, npPosY, npPosZ, npRad];
        mvmt = npMat .* 0;
        i = 100
        while(i < 100 && valid == 0)
            % 2.0 - Cooling loop - 3.0
            % 2.1 - Compute particle vectors and resulting motion vector from x,y,z
            for(i = 1:numNP)
                mX = (1/i) * sum(maxRad./(npPosX - npPosX(i)));
                mY = (1/i) * sum(maxRad./(npPosY - npPosY(i)));
                mZ = (1/i) * sum(maxRad./(npPosZ - npPosZ(i)));
                mvmt(i) = [mX, mY, mZ, 0];
                
                dists = sqrt((npPosX - npPosX(i)).^2 + (npPosY - npPosY(i)).^2 + (npPosZ - npPosZ(i)).^2);
                minDists = npRad + npRad(i);
                % 2.2 - check for valid geometry
                
            end
            i = i-1;
            npMat = npMat + mvmt;
        end
    % 4.0 - Return nanoparticle vector OR error message if no valid
    %       geometry is possible...
    % Refine coordinates, prioritize core over shells when forming omega 
    
    
    npCoords = npMat;
else
    %% Input Parameter Error Handling
    disp('algorithm type not recognized');
end

%% Output Placeholders
end

