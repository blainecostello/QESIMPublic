function [proxim,pathNP] = percAnalyzer(npList, gd)
% Percolation Analyzer Tool 
%   utilizes the recursive isPerc function to provide quantitative answers 
%   to several questions that arise with any given geometry 
%       1. Is it percolating?
%       2. What is the most likely breakdwown path?
%
%   This code first determines whether or not the geometry is
%   percolating between the electrodes.  If the code is percolating, it
%   simply returns the path and denotes percolation by setting 
%
%                       proxim = 0;
%
%   If the geometry is not found to be percolating, the analyzer continues
%   to search for potential breakdown paths by performing the recursive
%   percolation check with incrementally increasing max radii.  The extra
%   distance outside of the nanoparticle's geometric boundary is returned
%   as 'proxim' and the path returned will represent the most likely
%   breakdown path.  
%
%   The code uses this breakdown path to compute the boundaries of cubic
%   volumes within the simulated space that are 
%
% INPUTS
%   npList - List of nanoparticle positions & their radii in matrix [Nx4] 
%   gd     - Window size in Voxels, scalar representing z-height
%
% OUTPUTS
%   proxim  - proximity of farthest nanoparticles in breakdown path (0 = percolating)
%   pathNP  - list of nanoparticle positions for those included in the 
%
%
perc = 0;  % perc flag
proxim = -1;
% Determine Nanoparticles in contact with top electrode (if any) 
if(~isempty(npList))
    while (perc == 0)   % Only exit when percolating geom is found.
        proxim = proxim + 1;
        % Get temp nanoparticle list & adjust for proximity
        npListT = npList;
        npListT(:,4) = npListT(:,4) + proxim;
        % find top nanoparticles
        ind = (npListT(:,3) >= (gd - npListT(:,4)));
        ind = find(ind);
        i = length(ind);
        
        if ~isempty(ind)
            while (i > 0 && perc == 0) 
                first = npListT(ind(i),:);
                npListT(ind(i),:) = [];
                [perc, path] = isPerc(first,npListT);
                i = i - 1;
            end
        end
        if (perc == 1)
            pathNodes = size(path);
            if(pathNodes(1) > 2)
                pathNP = [path(2:pathNodes(1),:); first];
            else
                pathNP = [path;first];
            end
        end
    end
    
else % Handle geometry with no nanoparticles 
    pathNP = zeros(0,4);
end
end

