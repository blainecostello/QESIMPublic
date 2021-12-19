

days = {'Monday'; 'Tuesday'; 'Wednesday'; 'Thursday'; 'Friday'; 'Saturday'; 'Sunday'}

% All Data (50,000 days)
amts = [ -59.17712137672981, -2.3436303608749594,  -6.625477911006305, -21.307592185019587, 19.817919941329265,  20.143359561928328,  -50.24572286443639]
cnts =  [7663, 7401, 7607, 7618, 7618, 7618, 7623]

% ADA Data (181 Days)
amts2 = [0.063265, 0.003401, 0.128614, -0.055608, -0.206682, 0.343057, -0.374262]
cnts2 = [26,25,26,26,26,26,26];

% XLM Data (300 Days)
amts3 = [-0.008889, -0.012238, -0.099159, -0.123974,  0.545891, -0.090670, 0.068198]
cnts3 = [43,42,43,43,43,43,43]


figure
bar(amts./cnts*100)
set(gca, 'xtick', 1:7, 'xticklabel', days)
title("All Cryptos")

figure
bar(amts2./cnts2*100)
set(gca, 'xtick', 1:7, 'xticklabel', days)
title("ADA-USD")

figure
bar(amts3./cnts3*100)
set(gca, 'xtick', 1:7, 'xticklabel', days)
title("XLM-USD")