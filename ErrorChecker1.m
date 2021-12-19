% Convergence Checker for AnaValid3
clc
clear all
inc = 20;
max = 120;
for i = inc:inc:max
    error(round(i/inc)) = AnaValid3(i,-2)*100;
end

figure
plot([inc:inc:max], error, '-s',"Linewidth", 2);
title("Error between Analytical and Simulated Permittivity")
xlabel("Simulation Size")
ylabel("% error")
set(gcf, 'Position', [100 100 500 300])


