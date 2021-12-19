function [pCond,foundPath, bulkPath] = bdPath3(Omega, Mats, BDFS,Efld)
% initialize
%  Binary Omega Clones: Conductivity,
%
% Determine all starting voxels (top electrode)

%[gd,~,~] = size(Omega);

numVox = 1;
[gd, gd2, gd3] = size(Omega);
GD = [gd gd2 gd3];
Cond = zeros(gd,gd,gd-1);        % Initialize conductivity matrix
Visi = Cond;
maxP = 0.95;                      % maximum percentage to explore...
incP = 0.050;
minP = 0.00;
maxNV = round(maxP*gd*gd*gd);    % maximum number of voxels
evox = zeros(maxNV,1);            % set matrix to store energy of all voxels included in path
condVox = zeros(maxNV,4); % form conductivity matrix
temp = Efld(:,:,1:gd-1)./BDFS(Omega(:,:,1:gd-1)+1);
pFound = 0;

bdWeight = temp;


numSV = 0;

firstPathCond = 0;


% Increment percentage of included voxels
pCond = minP;
while pCond < maxP && pFound == 0
    pCond = pCond + incP;
    % Determine Number of Voxels included
    prevNumVox = numVox;
    numVox = round(GD(1)*GD(2)*GD(3)*pCond);
    
            
    % POPULATE CONDUCTIVITY MATRIX
    for(i = prevNumVox:numVox)
        [evox(i), Idx] = max(temp(:));
        [ix, iy, iz] = ind2sub(size(temp),Idx);
        
        % FIND ALL START VOXELS
        if(iz == gd-1)
            numSV = numSV + 1;
            startVox(numSV,:) = [ix,iy,iz];
        end
        
        % set conductivity to true at 'this' voxel
        Cond(ix,iy,iz) = 1;
        condVox(i,:) = [ix,iy,iz,1]';
        temp(ix,iy,iz) = -1;
    end
    
    disp("Plotting conductive voxel space")
    
    plotSparse([0,1], condVox(:,1:4), GD(1), GD(2), GD(3), 'fast');
    view([40,150]);
    title("Voxel Path-Search Space: " + pCond*100 + "% Conducting")
    pause(0.5)
    
    % Determine minP quickly with break detection algorithm
    layer = gd-1; % Z layer index 
    broken = 0;   % flag describing whether break in Cond has been detected
    nvoxVect = zeros(layer, 1); % vector containing the number of voxels for each Z layer
    botNeck = zeros(0,3);  % vector containing indices of all bottleneck points at several layers
    % Look through each Z layer of Cond
    while layer > 0 && broken == 0
        nvoxLay = sum(sum(Cond(:,:,layer)));
        nvoxVect(layer) = nvoxLay;
        disp("Num vox in Layer " + layer + ": " + nvoxLay);
        
        if nvoxLay < 1
            broken = 1;
        elseif nvoxLay == 1
            CondLay = Cond(:,:,layer);
            [evox(i), Idx] = max(CondLay(:));
            [ibx, iby] = ind2sub(size(CondLay),Idx);
            %botNeck = [botNeck; ibx,iby,layer]
        end
        % decrement  layer
        layer = layer - 1;
    end
    
    
    
    
    % construct temp variable and conductivity matrix
    
    % Display surface plot of top electrode.
    %surf(temp(:,:,gd-1))
    
    % Find location of all first voxels included in top electrode
   
    % (Debugging & testing)
    %plotSparse([0;1;2], condVox, gd, gd, gd-1, 'fast');
    %view([40,150])
    
    
    if(broken == 0)
        if firstPathCond == 0
            firstPathCond = pCond;
            fpcp = [nonzeros(condVox(:,1)),nonzeros(condVox(:,2)),nonzeros(condVox(:,3)),nonzeros(condVox(:,4))];% ones(length(condVox),1)];;
        end
        % ============ ======== =============
        % Initialize searching data structures
        % stack of ALL previously visited voxels (allocate large space for this)
%           Moved below for testing robustness
            visiVox = zeros(gd*gd,3);
            visiMat = zeros(gd,gd,gd-1);
%             % stack of voxels in single path (allocate less space for this)
        
        
        % thisvox (CURRENT POSITION)
        tv = [0,0,0];
        % DETERMINE STARTING VOXEL
        % 1. Check each Z layer for conducting voxels
        % 2. if All Z layers have at least one conducting voxel
        % 2.1. Number of layers with one single conducting voxel
        pathVox = zeros(0,3);
        pvOut = zeros(0,3);
        for SV = 1:numSV

            % stack of voxels in single path (allocate less space for this)
            
            
            % last voxel direction.  (initialized to 1 for top electrode)
            lvd = 1; % LEGEND: - [z+ = 1; z- = 2; y+ = 3; y- = 4; x+ = 5; x- = 6;];
            tv = startVox(SV,:);
            pathVox(1,:) = tv;
            numVis = 1;
            numPV = 1;
            cmax = (GD(1)*GD(2)*4);
            count = 0;
            %% Search for Path through "condMat"
            while pFound == 0 && ~isempty(pathVox) && count < cmax
                %             disp(tv(1))
                %             disp(tv(2))
                %             disp(tv(3))
                count = count + 1;
                numVis = numVis + 1;
                visiVox(numVis,:) = tv;
                disp("This voxel at layer: " + tv(3))
                visiMat(tv(1),tv(2),tv(3)) = 1;
                if tv(3) == 1
                    disp("We've found a path!  Now go look at it to calc a breakdown voltage!")
                    pFound = 1;
                    
                    pvOut = [pvOut;length(pathVox),length(pathVox),length(pathVox);pathVox];
                    %plotSparse([0;1;2], [pathVox,ones(length(pathVox),1)], gd, gd, gd-1, 'fast');
                    %view([40,150])
                else
                    % no path found yet, determine next voxel move
                    
                    % For comparison, store all adjacent voxel locations in a col vector
                    adjVox = [tv(1),tv(2),tv(3)-1, 1; tv(1),tv(2)+1,tv(3)-1, 2; ...
                        tv(1)+1,tv(2),tv(3)-1, 3; tv(1),tv(2)-1,tv(3)-1, 4; ...
                        tv(1)-1,tv(2),tv(3)-1, 5; ...
                        tv(1),tv(2)+1,tv(3), 6; tv(1)+1,tv(2)+1,tv(3), 7; ...
                        tv(1)+1,tv(2),tv(3), 8; tv(1)+1,tv(2)-1,tv(3), 9; ...
                        tv(1),tv(2)-1,tv(3), 10; tv(1)-1,tv(2)-1,tv(3), 11; ...
                        tv(1)-1,tv(2),tv(3), 12; tv(1)-1,tv(2)+1,tv(3), 13; ...
                        tv(1),tv(2),tv(3)+1, 14; tv(1),tv(2)+1,tv(3)+1, 15; ...
                        tv(1)+1,tv(2),tv(3)+1, 16; tv(1),tv(2)-1,tv(3)+1, 17; ...
                        tv(1)-1,tv(2),tv(3)+1, 18];
                    
                    for i = 1:length(adjVox)
                        adjVox(i,1) = circ(adjVox(i,1),gd);
                        adjVox(i,2) = circ(adjVox(i,2),gd);
                    end
                    
                    % initialize column vector voxWeight
                    
                    % cut off top nextVox's if thisVoxel is at top of window
                    if tv(3) == gd - 1
                        adjVox = adjVox(1:13,:);
                    end
                    
                    voxWeight = zeros(length(adjVox),1);
                    
                    
                    for vxnm = 1:length(adjVox)
                        loc = adjVox(vxnm,:);
                        
                        % check specified adjacent voxel location
                        % check if valid unexplored path voxel
                        if (visiMat(circ(loc(1),gd),circ(loc(2),gd),loc(3)) == 1) || (Cond(circ(loc(1),gd),circ(loc(2),gd),loc(3)) == 0)
                            voxWeight(vxnm) = -1; % -1 denotes invalid voxel
                        else
                            % determine weight (likelihood of taking this direction)
                            voxWeight(vxnm) = bdWeight(circ(loc(1),gd),circ(loc(2),gd),loc(3));
                        end
                        % compute maximum weighted direction
                        
                        
                    end
                    if max(voxWeight) > 0
                        [evox(i), Idx] = max(voxWeight);
                        [ix] = ind2sub(size(voxWeight),Idx);
                        
                        % weight for traverse toward electrode
                        dnVox = adjVox(1:5,:);
                        voxWeight(1:5) = (1./voxWeight(1:5)) .* 3.0;
                        
                        % weight for lateral traversal
                        ltVox = adjVox(6:13,:);
                        voxWeight(6:13) = (1./voxWeight(6:13)) .* 1.0;
                        
                        % if not at top electrode
                        if tv(3) ~= gd - 1
                            % weight for traverse away from electrode
                            upVox = adjVox(14:18,:);
                            voxWeight(14:18) = (1./voxWeight(14:18)) .* 0.33;
                        end
                        
                        % if path is available
                        pathVox = [pathVox; adjVox(ix,1:3)];
                        tv = adjVox(ix,1:3);
                        numPV = numPV + 1;
                    else
                        % if no other voxels are valid in adjacent voxel list, pop
                        % voxels off path stack and repeat
                        numPV = numPV - 1;
                        if numPV < 1
                            pathVox = [];
                        else
                            pathVox = pathVox(1:numPV,:);
                            tv = pathVox(numPV,:);
                        end
                    end
                end
            end
            
            %         pause(0.1)
            %%========================================
            
        end
        disp("SV Complete: " + SV + " of " + numSV);
    end
end

% find maximum value E-field/BDFS in pathVox - vMax is scaled using the
% inverse of this value
%(initial energy density approximations)
% EinP = zeros(length(pathVox),1)
% for i = 1:length(pathVox)
%     loc = pathVox(i,:);
%     EinP(i) = E2p(loc(1),loc(2),loc(3));
% end
% EDmax = (1./min(EinP))^2*Ed2p/(gd*gd*gd)
% EDmin = (1./max(EinP))^2*Ed2p/(gd*gd*gd)

%% Now that we have a path through the solid,
% for each layer in 'Z', find the maximum probabilistic breakdown voxel
% join non-adjacent breakdown voxels with lateral voxels to fill in the
% gaps
% what results is a single voxel line through the material depicting the
% predicted breakdown path
foundPath = [pathVox, ones(length(pathVox),1)];
bulkPath = [nonzeros(condVox(:,1)),nonzeros(condVox(:,2)),nonzeros(condVox(:,3)),nonzeros(condVox(:,4))];% ones(length(condVox),1)];

% [pvW,pvH] = size(pathVox);
% if(pvW == 0) % check if path vox is empty... 0x0
%     
%     % 1. convert condVox vector into 3D matrix
%     pathMat = zeros(gd,gd,gd-1);
%     for v = 1:length(pathVox)
%         % set all path voxels to the bd likelihood from bdWeight
%         pathMat(pathVox(v,1),pathVox(v,2),pathVox(v,3)) = bdWeight(pathVox(v,1),pathVox(v,2),pathVox(v,3));
%     end
%     % 2. iterate through each z layer and find path framework
%     thinPath = zeros(0,4);
%     for layer = 1:gd-1
%         thisLayer = pathMat(:,:,layer);
%         [~, Idx] = max(thisLayer(:));
%         [ix, iy] = ind2sub(size(thisLayer),Idx);
%         
%         thinPath = [thinPath;ix,iy,layer,pathMat(ix,iy,layer)];
%     end
%     
%     plottablePath = [thinPath(:,1:3),ones(length(thinPath),1)];
%     
%     %       plotSparse([0;1;2], bulkPath, gd, gd, gd-1, 'fast');
%     %       view([40,150])
%     
%     %% 3. populate gaps
%     maxVox = max(thinPath(:,4))
%     % for each layer, (connect the dots algorithm)
%     addVox = zeros(0,4);
%     for layer = 1:gd-2
%         % check if there is a gap at this layer
%         % get distance of gap at this layer
%         % if there is a gap, fill the gap in
%         % make sure to check to make sure all gaps are filled by voxels
%         % contained within the 'condVox' list.
%         % if gap is too large, skip layer (continue in direction of greatest bdWeight...)
%         % else, do nothing & continue
%         this = thinPath(layer,1:3);
%         next = thinPath(layer + 1,1:3);
%         dist = [(next(1)-this(1)), (next(2)-this(2)), (next(3)-this(3))];
%         
%         if abs(dist(1)) > gd/2
%             if (this(1) > next(1))
%                 next(1) = next(1) + gd;
%             else
%                 next(1) = next(1) - gd;
%             end
%         end
%         if abs(dist(2)) > gd/2
%             if (this(2) > next(2))
%                 next(2) = next(2) + gd;
%             else
%                 next(2) = next(2) - gd;
%             end
%         end
%         
%         dist = [(next(1)-this(1)), (next(2)-this(2)), (next(3)-this(3))];
%         
%         % calculate number of voxels between these two points
%         % calculate travel distance in X and Y per loop (diagonal line)
%         x1 = this(1);
%         x2 = next(1);
%         y1 = this(2);
%         y2 = next(2);
%         
%         dx = x2 - x1;
%         dy = y2 - y1;
%         if dx ~= 0
%             for x = x1:x2
%                 y = y1 + dy * (x - x1) / dx;
%                 newVox = [x,y,layer,1];
%                 addVox = [addVox;newVox];
%             end
%         end
%         
%     end
%     
%     for i = 1:length(addVox)
%         addVox(i,1) = circ(addVox(i,1),gd);
%         addVox(i,2) = circ(addVox(i,2),gd);
%     end
    
    if pCond == maxP
        % No path found... Using first pass conductive possibility case (at least one vox in each z layer is conducting)
        foundPath = fpcp;
        pCond = firstPathCond;
    end

    %thinPath = [thinPath;addVox];
    %plottablePath = [thinPath(:,1:3),ones(length(thinPath),1)];
    %bulkPath = plottablePath;
    %%
    %condVox;
    
end

