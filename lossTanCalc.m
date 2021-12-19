        

E2p = Eapp 
powTot = Eapp^2 .* 1.5E-8 * GD(1)*GD(2)*GD(3)*Delta(1)*Delta(2)*(Delta(3))
disp(powTot)
cond_fromPow = powTot / (Eapp^2*Delta(3)*(GD(3)-1)) 
disp(cond_fromPow)
% recalculate effective permittivity
% eff_recalc = real(effSp) - 1i*(real(cond_recalc))/(freq * e0)
eff_recalc = real(effRe) - 1i*(abs(cond_fromPow))/(freq * e0)
        
% reccalculate loss tangent 
tanD_recalc = abs(imag(eff_recalc)/real(eff_recalc)) 

%%  Hand calculation** 
%   Back-calculate conductivity from loss tangent and real permittivity
%   from first measured datapoint... 
%   (Loss tan = 0.013)
%   (Real Permittivity = 2.5)
%   ()