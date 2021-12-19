function [Exyz, E, V] = electrifyQE(Omega, Voltage, cP, dVol)
% Electrify - Induces a potential across the material in the Z direction.
%   Computes the electric field using FEM to generate the laplace equation
%   coefficients.
%   
%==================={ ARGUMENTS }=====================%
%   Omega        = 3-space grid of material voxels           [X by Y by Z]
%   Voltage      = [Vx-,Vx+;Vy-,Vy+;Vz-,Vz+];                [3 by 2]
%               (-1 = circular b.c.) Ex: [-1,-1;-1,-1;0,1] 
%   cP           = complex permittivity of each material     [#mats by 2]
%   sigma        = conductivity of each material
%   freq         = frequency in rad/sec
%   dVol         = delta x,y,z, or cubic voxel dimensions    [3 by 1]
%
%==================={ OUTPUTS }=====================%
%   E      = 3D complex electric field matrix     [X by Y by Z by 2] 
%  
%
%=============={ ABOUT THE ALGORITHMS }===============%
%
%   
% Initialize Problem size, Vectors and Matrices
[X, Y, Z] = size(Omega);

N = X*Y*(Z);                              % Problem size
% K = 8.61733035*10^-5;                   % boltzmanns
A = sparse([],[],[],N,N,(N*7));         % SPARSE Matrix A
b = sparse([],[],[],N,1,(X*Y*6));       % SPARSE Matrix b

% complex permittivity
% cP = (K*eps+1i*(sigma/freq));

% 1. Convert Omega [X, Y, Z] to coefficeint vector of length [X*Y*Z, 1]
%    for Omega(i,j,k) for i = 1:X, j = 1:Y, k = 1:Z.  Values at each index
%    location will be equal to _________
% handle boundary index conditions:
%                Dirichlet:  if(k+1 > Z)->k+1 = Vdd ; if(k-1 < 1)->k-1 = 0
%                Circular:   if(j+1 > Y)->j+1 = 1   ; if(j-1 < 1)->j-1 = Y
%                            if(i+1 > X)->i+1 = 1   ; if(i-1 < 1)->i-1 = X
% 
% Compress indexing: I = i + j(X) + k(X*Y)
% Reverse indexing: k = floor(I/X*Y), j = floor((I%(X*Y))/X), k = (I%(X*Y))

% Handle Omega Boundary Conditions
% Concat X Axis Boundary Condition 
if((Voltage(1,1) == -1) && (Voltage(1,2) == -1))
  % Circular Boundary Conditions
  % [ Omega(X,:,:) ; Omega ; Omega(1,:,:) ]
    Omega = cat(1, Omega(X,:,:), Omega);
    Omega = cat(1, Omega, Omega(1,:,:));
else
  % Electrode Boundary Conditions
  % [ zeros(Y,Z) + Voltage(1,1) ; Omega ; zeros(Y,Z) + Voltage(1,2) ]
    Omega = cat(1, zeros(1,Y,Z)+1, Omega);
    Omega = cat(1, Omega, zeros(1,Y,Z)+1);
end
% Concat Y Axis Boundary Condition
[X, Y, Z] = size(Omega);
if((Voltage(2,1) == -1) && (Voltage(2,2) == -1))
  % Circular Boundary Conditions
  % [Omega(:,Y,:),Omega,Omega(:,1,:)];
    Omega = cat(2, Omega(:,Y,:), Omega);
    Omega = cat(2, Omega, Omega(:,1,:));
else
  % Electrode Boundary Conditions
  % [zeros(Y,Z)+Voltage(2,1);Omega;zeros(Y,Z)+Voltage(2,2)];
    Omega = cat(2, zeros(X,1,Z) + Voltage(2,1), Omega);
    Omega = cat(2, Omega, zeros(X,1,Z) + Voltage(2,2));
end
% Concat Z Axis Boundary Condition
[X, Y, Z] = size(Omega);
if((Voltage(3,1) == -1) && (Voltage(3,2) == -1))
  % Circular Boundary Conditions
  % [Omega(:,:,Z);Omega;Omega(:,:,1)];
    Omega = cat(3, Omega(:,:,Z), Omega);
    Omega = cat(3, Omega, Omega(:,:,1));
else
  % Electrode Boundary Conditions
  % [zeros(Y,Z)+Voltage(3,1);Omega;zeros(Y,Z)+Voltage(3,2)];
    Omega = cat(3, zeros(X,Y,1)+1, Omega);
    Omega = cat(3, Omega, zeros(X,Y,1)+1);
end
I = 1;

[X, Y, Z] = size(Omega);
X = X - 2;
Y = Y - 2;
Z = Z - 2;
disp('FEM... (Generating A matrix...)')
for(k = 2:Z+2)
    %disp("FEM Z Index = " + k)
    for(j = 2:Y+1)
        for(i = 2:X+1)
        % disp(i + ", " + j + ", " + k)
      % Get permittivity of each nonzero in I'th row of A matrix 
         eijk    = cP(Omega(i-1,j-1,k-1) + 1);
         eijpk   = cP(Omega(i-1,j,k-1) + 1);
         eijpkp  = cP(Omega(i-1,j,k) + 1);
         eijkp   = cP(Omega(i-1,j-1,k) + 1);
         eipjk   = cP(Omega(i,j-1,k-1) + 1);
         eipjpk  = cP(Omega(i,j,k-1) + 1);
         eipjpkp = cP(Omega(i,j,k) + 1);
         eipjkp  = cP(Omega(i,j-1,k) + 1);  
         
         
         dimCx = ((dVol(2)*dVol(3))/(dVol(1)*4));
         dimCy = ((dVol(1)*dVol(3))/(dVol(2)*4));
         dimCz = ((dVol(1)*dVol(2))/(dVol(3)*4));

         
       val = [(-(dimCx + dimCy + dimCz)*(eijk + eipjk + eijpk + eijkp +  eipjpk  + eipjpkp + eipjkp + eijpkp)), ...
                (dimCx * (eipjkp  + eipjk  + eipjpkp   + eipjpk )), (dimCx*(eijkp + eijk + eijpkp  + eijpk  )), ...
                (dimCy * (eijpkp + eipjpkp + eijpk  + eipjpk)), (dimCy*(eijkp  + eipjkp  + eijk   + eipjk   )), ...
                (dimCz * (eijpkp   + eipjpkp   + eijkp  + eipjkp)), (dimCz*(eijpk  + eipjpk  + eijk + eipjk ))];    
            %val = [1,1,1,1,1,1,1]
            dm = dimCx + dimCy + dimCz;

     % X Boundary Check
         if((Voltage(1,1) ~= -1) || (Voltage(1,2) ~= -1))
             % Electrode Boundary Condition: 
             if(i == X+1)
                 % Omega(:,:,Z+1) = Voltage(3,2);
                  volts = Voltage(1,2) ...
                        + ( Voltage(1,2) * 2 * (X-2/(X-2)) ) ... 
                        + ( Voltage(1,2) * (X-3/(X-2)) ) ;
                  b = b + sparse(I,1,Voltage(3,2), N,1,(X*Y*6));
             elseif(i == 2)
                  %volts = 0;
                  %b = b + sparse(I,1,-Voltage(3,2), N,1,(X*Y*6));
             end
             
         end
     % Y Boundary Check
         if((Voltage(2,1) ~= -1) || (Voltage(2,2) ~= -1))
             % Electrode Boundary Condition: 
             if(j == Y+1)
                 % Omega(:,:,Z+1) = Voltage(3,2);
                  volts = Voltage(2,2) ...
                        + ( Voltage(2,2) * 2 * (Y-2/(Y-2)) ) ... 
                        + ( Voltage(2,2) * (Y-3/(Y-2)) ) ;
                  b = b + sparse(I,1,Voltage(2,2), N,1,(X*Y*6));
             elseif(j == 2)
                  %volts = 0;
                  %b = b + sparse(I,1,-Voltage(3,2), N,1,(X*Y*6));
             end
         end
     % Z Boundary Check
         if((Voltage(3,1) ~= -1) || (Voltage(3,2) ~= -1))
             % Electrode Boundary Condition: 
             dV = Voltage(3,2)-Voltage(3,1);
             if(k == Z+1)
              % Omega(:,:,Z+1) = Voltage(3,2);
                volts = dV;%(dV*(Z))/(Z+1);   %( dV + ( dV * 5 * ((k-2)/(Z)) )   ... 
                         % + ( dV * ((k-3)/(Z)) ) )/7 ;
                %bc = 1.72*10^-8;
                b = b + sparse(I,1,Voltage(3,2), N,1,(X*Y*6));
                val = [1,0,0,0,0,0,0];

             elseif(k == Z+2)
                 % set A matrix entities to zero... ( ijk+ -> top electrode)
                 %val = [1,0,0,0,0,0,0];
             elseif(k == 2)
                % volts = 0;
                volts = 0;%(dV/(Z));  %(dV + ( dV * 5 * (((k-2)/(Z))) ))/6 ;
                b = b + sparse(I,1,Voltage(3,1), N,1,(X*Y*6));
                % set A matrix entities to zero... ( ijk- -> bottom electrode)
                val = [1,0,0,0,0,0,0];

             end
         end
         
 
         
%     Dimension coefficients
%          dimCx = ((dVol(1)*dVol(2)*dVol(3))/(4));
%          dimCy = ((dVol(1)*dVol(2)*dVol(3))/(4));
%          dimCz = ((dVol(1)*dVol(2)*dVol(3))/(4));
         

      % compute index of each nonzero in I'th row of A matrix
      x = i-1;
      y = j-1;
      z = k-1;
%          ijk    = I;
%          ipjk   = circ(x+1,X) + (y-1)*(X) + (z-1)*((X)*(Y));
%          injk   = circ(x-1,X) + (y-1)*(X) + (z-1)*((X)*(Y));
%          ijpk   = (x) + circ( y ,Y)*(X) + (z-1)*((X)*(Y));
%          ijnk   = (x) + circ(y-2,Y)*(X) + (z-1)*((X)*(Y));
%          ijkp   = (x) + (y-1)*(X) + (circ( z ,Z))*((X)*(Y));
%          ijkn   = (x) + (y-1)*(X) + (circ(z-2,Z))*((X)*(Y));
%     working Indices 
      ijk = I;
      ipjk = circ(x+1, X) + (y-1)*(X) + (z-1)*((X)*(Y));
      injk = circ(x-1, X) + (y-1)*(X) + (z-1)*((X)*(Y));
      ijpk = circ(x + ( y )*X, X*Y) + (z-1) * ((X)*(Y));
      ijnk = circ(x + (y-2)*X, X*Y) + (z-1) * ((X)*(Y));
      ijkp = circ(x + (y-1) * X + ( z )*X*Y, X*Y*Z);
      ijkn = circ(x + (y-1) * X + (z-2)*X*Y, X*Y*Z);
       % X Perturbation
   
%          A(I,circ(ijk,I))  = -(dimCx + dimCy + dimCz)*(eijk + eipjk + eijpk + eijkp +  eipjpk  + eipjpkp + eipjkp + eijpkp);
%          A(I,circ(ipjk,I)) = dimCx*(eijpkp  + eijpk  + eijkp   + eijk);
%          A(I,circ(injk,I)) = dimCx*(eipjpkp + eipjpk + eipjkp  + eipjk);
%          A(I,circ(ijpk,I)) = dimCy*(eipjpkp + eijpkp + eipjpk  + eijpk);
%          A(I,circ(ijnk,I)) = dimCy*(eipjkp  + eijkp  + eipjk   + eijk);
%          % Z electrode boundary...
%          A(I,circ(ijkp,I)) = dimCz*(eipjk   + eijk   + eipjpk  + eijpk);
%          A(I,circ(ijkn,I)) = dimCz*(eipjkp  + eijkp  + eipjpkp + eijpkp);
        %disp(i + ", " + j + ", " + k);
        %disp("|c = "+ijk+"|x- = "+injk+"|x+ = "+ipjk+"|y- = "+ijnk+"|y+ = "+ijpk+"|z- = "+ijkn+"|z+ = "+ijkp);
       row = [I;I;I;I;I;I;I];
       col = [ijk, ipjk, injk, ijpk, ijnk, ijkp, ijkn];
       
%        val = [(-(dimCx + dimCy + dimCz)*(eijk + eipjk + eijpk + eijkp +  eipjpk  + eipjpkp + eipjkp + eijpkp)), ...
%                 (dimCx * (eipjkp  + eipjk  + eipjpkp   + eipjpk )), (dimCx*(eijkp + eijk + eijpkp  + eijpk  )), ...
%                 (dimCy * (eijpkp + eipjpkp + eijpk  + eipjpk)), (dimCy*(eijkp  + eipjkp  + eijk   + eipjk   )), ...
%                 (dimCz * (eijpkp   + eipjpkp   + eijkp  + eipjkp)), (dimCz*(eijpk  + eipjpk  + eijk + eipjk ))];
          if(k < Z+2)
                %disp(N)
                %disp(col)
              A = A + sparse(row, col, val, N,N,(N*7));
             % Now increment I -> I = i + j(X) + k(X*Y)
               I = I + 1;
          end
        end
    end
end

%disp("Solving...")
% 2. Solve Transformation for V [X*Y*Z, 1]
    %volts = ThomasLU(A,b);
    %volts = conjgrad(A,b);
    %volts = pcg(A,b);
    clear Omega,
    disp("Crunching Ax=b ... [ Solve for x using mldivide() ]")
    %dA = decomposition(A,'lu');
    tic
    volts = A\b;
    toc
    clear A b
    %figure 
    %image((full(A)~=0)*255)
    disp("Operation complete, Now let's finish up calculating the electric fields")

    
%disp("Converting...")
% 3. Convert V [X*Y*Z, 1] back to E [X, Y, Z]
x = 1; 
y = 1;
z = 1;
V = zeros(X,Y,Z);
%B = zeros(X,Y,Z);
for(i = 1:I-1)
   %disp("I: " + i + " | xyz:" + x + ", " + y + ", " + z)
   % build 3d matrix
   V(x,y,z) = volts(i,1);
%   B(x,y,z) = b(i,1);
   % increment indices
   if(x == X)  
      if(y == Y)
        z = z + 1;
        y = 0;
      end
      x = 0;
      y = y + 1;
   end
   x = x + 1;
end

%figure
%plotSparseE(B);
%figure
%surf(abs(A))


%% Compute Electric Field:
% concat XYZ - [voltage [0,1] across Z]
vtemp = cat(1, V(X,:,:), V);
vtemp = cat(1, vtemp, V(2,:,:));
vtemp = cat(2, vtemp(:,Y,:), vtemp);
vtemp = cat(2, vtemp, vtemp(:,2,:));
vtemp = cat(3, zeros([X+2,Y+2,1])+Voltage(3,1), vtemp);
vtemp = cat(3, vtemp, zeros([X+2,Y+2,1])+Voltage(3,2));

% figure
% plotSparseE(vtemp,[0,1]);
%disp("Gradient Calculation...");
dV = zeros([X,Y,Z,3]);
Emag = zeros([X,Y,Z-1]);
for(x = 2:X+1)
    for(y = 2:Y+1)
        for(z = 2:Z)
            dVdx = (vtemp(x+1,y,z) - vtemp(x,y,z))/dVol(1);
            dVdy = (vtemp(x,y+1,z) - vtemp(x,y,z))/dVol(2);
            dVdz = (vtemp(x,y,z+1) - vtemp(x,y,z))/dVol(3);
            dV(x,y,z,:) = [dVdx,dVdy,dVdz];
            Emag(x-1,y-1,z-1) = (dVdx.^2 + dVdy.^2 + dVdz.^2)^(0.5);
        end
    end
end
E = Emag;%-dV(2:X,2:Y,2:Z,:);
Exyz = dV(2:X+1, 2:Y+1, 2:Z, :);

disp("Ok, all calculated, {electrifyQE} has successfully completed execution")

% plot vector field
% figure
% plotSparseE(Emag);

end


