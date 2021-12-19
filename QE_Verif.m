% New Electrify Tester (AC / Quasi-Electrostatics)

% 1) Load Geometry & Set input variables for ElectrifyQE
        % Material data obtained from: <url>
        Omega = zeros([20 20 20]);
        Volt = [-1,-1;-1,-1;0,1];
        [Xsz, Ysz, Zsz] = size(Omega);
        freq = 2*pi/10;
        sigma = 2.2E-14;  % 1/Ohm-Meters
        
        e0 = 8.8541878128E-12;
        Mats = [3.9 - 1i*(sigma/(freq*e0)), 3.9 - 1i*(sigma/(freq*e0))];

        % Base Thickness of Sample:
        tSamp = 500E-9;
        Delta = [tSamp/Xsz tSamp/Ysz tSamp/Zsz];

% 2) Run ElectrifyQE 
        [Exyz, E, V] = electrifyQE(Omega, Volt, Mats, Delta);
% 3) Calculate Loss Tangent & Q Factor
        % Total Energy:
           EinEa = ((1/2)*e0 .* Mats(Omega(:,:,1:Zsz-1)+1).* E.^2) .* (Delta(1)*Delta(2)*Delta(3)); % energy in each voxel
           Ed = sum(EinEa, 3);
           Ed = sum(Ed, 2);
           Etot = sum(Ed, 1);
           eff = (2*Etot*Delta(3)*(Zsz-1))/((1^2)*e0*Delta(1)*Xsz*Delta(2)*Ysz);
           
           disp("Effective Permittivity:     " + eff)
           disp("Material Permittivity:      " + Mats(1))

           tanDel = -1 * imag(eff)./real(eff);
            
           disp("Loss Tangent:               " + tanDel)
           
           sigeff = -1 * imag(eff) * freq * e0;
        
           
           
           rhoeff = 1/sigeff;
           Reff = rhoeff * (Zsz*Delta(3))/1; 
           
           Rvox = freq * e0 * sum(-1./imag(Mats(Omega(:,:,:)+1)),3);
           
           IleakVox2 = sum(sum(1./Rvox .* sum(abs(Exyz(:,:,:,3)),3),2),1);
           
           Ileak = 1/Reff; % Amps per square meter [m^2]
           
           disp("Effective Leakage Current:  " + Ileak + " [Amps/m^2]")
           
           Vvox = 1/(Zsz);
           
           IleakVox = sum(sum( ( freq * e0 * imag(Mats(Omega(:,:,1)+1)).*abs(Exyz(:,:,1,3))  *  (-1) ) ,2),1); 
           
           disp("Leakage Current (Voxels):   " + IleakVox + " [Amps/m^2]")
           %disp("Leakage Current (Voxels):   " + IleakVox2 + " [Amps/m^2]")
           figure
           imagesc(squeeze(rot90(real(E2p(10,:,:)))), [real(minE2p*0.75),real(maxE2p)])

        % Effective Permittivity
            
        
        % Loss Tangent
        
        