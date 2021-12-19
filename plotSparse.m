function plotSparse(mats, coords, Xmax, Ymax, Zmax, mode)
%  Plots a voxelized visual representation of a 3d ndSparse matrix
%   coords represents a Nx3 matrix where each column contains all data for
%   x,y,z respectively.
%   Xmax Ymax and Zmax are the boundaries of the sparse matrix
%   Mode = 'slow' || 'fast'
pos = [1,1,1]; % position to plot single voxel **Loop this to print all entries of sparse array
vsize = [1,1,1];

% Material Colors
cT = [1, 0, 0;
      %1, 0.75, 0;
      0.65, 0.65, 0.65;
      1, 1, 0;
      0, 1, 0;
      0, 1, 1];
  
if(mode == 'slow')
    [r,c] = size(coords)
    

    
    for(i = 1:r)
        x = coords(i,1);
        y = coords(i,2);
        z = coords(i,3);
        c = coords(i,4);
        pos = [x,y,z];
        disp(i/r)
        % Color Transformations
        if(c(i) > 0)
            R = cT(c(i),1);
            G = cT(c(i),2);
            B = cT(c(i),3);
            color = [R,G,B];
            
            alpha = 0.2;
            plotcube(vsize,pos,alpha,color);
        end
    end
    

elseif(mode == 'fast')
    X = round(coords(:,1));
    Y = round(coords(:,2));
    Z = round(coords(:,3));
    M = coords(:,4);
    
    figure
    hold on
    for(i=1:length(mats))
        scatter3(X(M==mats(i)),Y(M==mats(i)),Z(M==mats(i)),'MarkerEdgeColor','k','MarkerFaceColor',cT(i,:));
    end
    
end
    axis([0, Xmax+1, 0, Ymax+1, 0, Zmax+1]);
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    pbaspect([1 1 1]);
end

