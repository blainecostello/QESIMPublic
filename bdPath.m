function [bulkPath,foundPath] = bdPath(Omega,Mats,BDFS,Efld)
% initialize 
%  Binary Omega Clones: Conductivity, 
% 
% Determine all starting voxels (top electrode)
%load('/Users/darvidbigardly/Desktop/Stuff & Things/PBS_Test/limit_test/NP_parallel/PercRuns/Run_1100/id_203-18-4-3:25-65-32.mat')
%load('/Users/darvidbigardly/Desktop/Stuff & Things/PBS_Test/limit_test/NP_parallel/PercRuns/Run_1100/id_9-1-4-5:25-65-10.mat')
%load('/Users/darvidbigardly/Desktop/Stuff & Things/PBS_Test/limit_test/NP_parallel/PercRuns/Run_1100/id_1007-22-4-40:25-65-2.mat')





%% remove this section (not included in function)
% Efld = E2p;
% Vfld = volt;
% om1 = OmegaBT;
% Omega = OmegaBT
% Mats = [11.8,50,80]
% BDFS = [380,100,100];
% numVox = round(0.05*(gd*gd))
%a = voxbd(om1,mats,Efld,Volt);
%b = bdSim5(om1,mats,bdfs,Efld)
%% 

%[gd,~,~] = size(Omega);

numVox = 1;
[gd gd2 gd3] = size(Omega);
GD = [gd gd2 gd3];
Cond = zeros(gd,gd,gd-1);        % Initialize conductivity matrix
Visi = Cond;
maxP = 1;                      % maximum percentage to explore...
incP = 0.005;
minP = 0.01;
maxNV = round(maxP*gd*gd*gd);    % maximum number of voxels 
evox = zeros(maxNV,1);            % set matrix to store energy of all voxels included in path
condVox = zeros(maxNV,4); % form conductivity matrix
temp = Efld(:,:,1:gd-1)./BDFS(Omega(:,:,1:gd-1)+1);
pFound = 0;


bdWeight = temp;

% Increment percentage of included voxels 
pCond = minP;
while pCond < maxP && pFound == 0
    pCond = pCond + incP;
  % Determine Number of Voxels included
    prevNumVox = numVox;
    numVox = round(gd*gd*gd*pCond);
    
  % construct temp variable and conductivity matrix
    
  % Display surface plot of top electrode.
    %surf(temp(:,:,gd-1)) 
    
  % Find location of all first voxels included in top electrode 
  numSV = 0;
    for(i = prevNumVox:numVox)
        [evox(i), Idx] = max(temp(:));
        [ix, iy, iz] = ind2sub(size(temp),Idx);
        
      % store locations of start Voxels
        if(iz == gd-1)
            numSV = numSV + 1;
            startVox(numSV,:) = [ix,iy,iz];
        end
        
      % set conductivity to true at 'this' voxel
        Cond(ix,iy,iz) = 1;
        condVox(i,:) = [ix,iy,iz,1]';
        temp(ix,iy,iz) = -1;
    end
    
    
 %%
    % Print Conductivity Voxels... (For debugging & testing)
%       plotSparse([0;1;2], condVox, gd, gd, gd-1, 'fast');
%       view([40,150])

 %%
 
 
    % ------------ -------- ------------%  
      condMat = zeros(gd,gd,gd);
      % encode condVox to condMat
      for i = 1:length(condVox)
         condMat(nonzeros(condVox(i,1)),nonzeros(condVox(i,2)),nonzeros(condVox(i,3))) = 1;
      end
    % ============ ======== =============
    % Initialize searching data structures
        % stack of ALL previously visited voxels (allocate large space for this)
        visiVox = zeros(gd*gd,3);
        visiMat = zeros(gd,gd,gd-1);
        % stack of voxels in single path (allocate less space for this)
        pathVox = zeros(0,3);

        % thisvox (CURRENT POSITION)
        tv = [0,0,0];
     % determine starting voxel
     pvOut = zeros(0,3);
     for SV = 1:numSV  
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
            visiMat(tv(1),tv(2),tv(3)) = 1;
            if tv(3) == 1
                disp("We've found a path!  Now go look at it to calc a breakdown voltage!")
                pFound = 1;
                
                pvOut = [pvOut;length(pathVox),length(pathVox),length(pathVox);pathVox];
%                 plotSparse([0;1;2], [pathVox,ones(length(pathVox),1)], gd, gd, gd-1, 'fast');
%                 view([40,150])
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
                    voxWeight(1:5) = voxWeight(1:5) .* 3.0;

                    % weight for lateral traversal
                    ltVox = adjVox(6:13,:);
                    voxWeight(6:13) = voxWeight(6:13) .* 1.0;

                    % if not at top electrode
                    if tv(3) ~= gd - 1
                        % weight for traverse away from electrode
                        upVox = adjVox(14:18,:);
                        voxWeight(14:18) = voxWeight(14:18) .* 0.33;
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
        %disp("SV complete")
%         pause(0.1)
        %%========================================
      
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

[pvW,pvH] = size(pathVox);
if(pvW == 0)
    
% 1. convert condVox vector into 3D matrix
pathMat = zeros(gd,gd,gd-1);
for v = 1:length(pathVox)
   % set all path voxels to the bd likelihood from bdWeight
    pathMat(pathVox(v,1),pathVox(v,2),pathVox(v,3)) = bdWeight(pathVox(v,1),pathVox(v,2),pathVox(v,3));
end
% 2. iterate through each z layer and find path framework
thinPath = zeros(0,4);
for layer = 1:gd-1
        thisLayer = pathMat(:,:,layer);
        [~, Idx] = max(thisLayer(:));
        [ix, iy] = ind2sub(size(thisLayer),Idx);
        
        thinPath = [thinPath;ix,iy,layer,pathMat(ix,iy,layer)];
end

plottablePath = [thinPath(:,1:3),ones(length(thinPath),1)];

%       plotSparse([0;1;2], bulkPath, gd, gd, gd-1, 'fast');
%       view([40,150])

%% 3. populate gaps
maxVox = max(thinPath(:,4))
% for each layer, (connect the dots algorithm)
addVox = zeros(0,4);
for layer = 1:gd-2
    % check if there is a gap at this layer
    % get distance of gap at this layer
        % if there is a gap, fill the gap in
            % make sure to check to make sure all gaps are filled by voxels
            % contained within the 'condVox' list.
        % if gap is too large, skip layer (continue in direction of greatest bdWeight...)
        % else, do nothing & continue
        this = thinPath(layer,1:3);
        next = thinPath(layer + 1,1:3);
        dist = [(next(1)-this(1)), (next(2)-this(2)), (next(3)-this(3))];
        
        if abs(dist(1)) > gd/2
            if (this(1) > next(1))
                next(1) = next(1) + gd;
            else
                next(1) = next(1) - gd;
            end
        end
        if abs(dist(2)) > gd/2
            if (this(2) > next(2))
                next(2) = next(2) + gd;
            else
                next(2) = next(2) - gd;
            end
        end
        
        dist = [(next(1)-this(1)), (next(2)-this(2)), (next(3)-this(3))];
                
        % calculate number of voxels between these two points
        % calculate travel distance in X and Y per loop (diagonal line)
        x1 = this(1);
        x2 = next(1);
        y1 = this(2);
        y2 = next(2);
        
        dx = x2 - x1;
        dy = y2 - y1;
        if dx ~= 0
            for x = x1:x2 
              y = y1 + dy * (x - x1) / dx;
              newVox = [x,y,layer,1];
              addVox = [addVox;newVox];
            end
        end
        
end
 
        for i = 1:length(addVox)
            addVox(i,1) = circ(addVox(i,1),gd);
            addVox(i,2) = circ(addVox(i,2),gd);
        end
        
thinPath = [thinPath;addVox];
plottablePath = [thinPath(:,1:3),ones(length(thinPath),1)];
bulkPath = plottablePath;
%% 
o1 = condVox;

end
end
