% Load Geometry File, specify first np, remove from list
function [perc, path] = isPerc(thisNP,npList)
% isPerc - Is Percolating?  
% Determines whether or not a geometry (defined by NP list "npPos") 
% contains at least one percolating path through recursion from the
% uppermost nanoparticle contacting the top electrode, down to the 
% bottom electrode via whatever available path.  A path is defined as 
% 
% Input:
%   thisNP - Current Nanoparticle  - 4x1 matrix  -  [x,y,z,r]
%   npList - list of nanoparticles - 4xn matrix  -  [...; x,y,z,r ;...]
% 
% Output:
%   perc - Boolean variable true percolating path exists between electrodes
%   path - list of nanoparticles included in first path.
%   
% 

% Initialize & Compute temp variables  
% (Avoid unnecessary variables to save space)
% make all entries unique (delete redundant nanoparticles)
npList = unique(npList,'rows');
% 
if (isempty(thisNP))
    perc = 0;
    path = zeros(0,4);
elseif (thisNP(3) - thisNP(4) < 1)
    perc = 1;
    path = thisNP;
else
% get list of only nanoparticles overlapping with current nanoparticle
    touchInd = npList(:,4)+thisNP(4) >= ((npList(:,1)-thisNP(1)).^2+(npList(:,2)-thisNP(2)).^2+(npList(:,3)-thisNP(3)).^2).^0.5;
    touchInd = find(touchInd);
    numT = length(touchInd);
    path = zeros(0,4);
% Counted While LOOP: 
%   remove  nanoparticle  from List & RECURSE
    c = numT;
    perc = 0;
    if(~isempty(touchInd))
        while ((c > 0) && (perc == 0))
%           save current as temp variable
            next = npList(touchInd(c),:);
%           remove next NP from np list
            npList(touchInd(c),:) = [];
% -----     RECURSION      -----
            [perc, path] = isPerc(next,npList);
            if perc == 1
                path = [path;next]; 
            else
                npList = [npList;next];
            end
%           update count
            c = c - 1;
        end
    end
end
end



