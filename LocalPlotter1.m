% Local Plotter

figure
hold on
plot(VolFrac, Meff, '-ok', 'linewidth', 2)
plot(SVF, SeffExt, 'o', 'linewidth', 2)
title("Effective Composite Permittivity")
xlabel("Volume Fraction of Inclusions")
ylabel("Permittivity")
legend("Measured Data", "Simulated Data", "Location", "Northwest")
set(gcf, 'Position', [0 300 300 200])
%ylim([0,10])


figure
hold on
plot(VolFrac, MTanD, '-ok', 'linewidth', 2)
plot(SVF, StanDExt, 'o', 'linewidth', 2)
title("Loss Tangent")
xlabel("Volume Fraction of Inclusions")
ylabel("Loss Tangent")
%ylim([0, 0.1])

legend("Measured Data", "Simulated Data", "Location", "Northwest")
set(gcf, 'Position', [300 300 300 200])
ylim([0, 0.05])



figure
hold on
plot(VolFrac, MBDFS, '-ok', 'linewidth', 2)
plot(SVF, SBDFSExt, 'o', 'linewidth', 2)
title("Breakdown Field Strength")
xlabel("Volume Fraction of Inclusions")
ylabel("Breakdown Strength")
%ylim([100, 450])
legend("Measured Data", "Simulated Data", "Location","Northeast")
set(gcf, 'Position', [600 300 300 200])
ylim([0,500])
