% bootstrap to have a confidence interval on seabed fraction.

frac = [0.0698203;
0.0503932;
0.0473392;
0.0630356;
0.0592118;
0.0634536;
0.0582511;
0.0640809;
0.0606762;
0.0681157;
0.0653511;
0.0641338;
0.068262;
0.0703444;
0.0654905;
0.068336;
0.0628307;
0.0657623;
0.0622844;
0.0660646;
0.0641322;
0.0716378;
0.0608331;
0.0628977;
0.0581741;
0.0586287];


% -- bootstrap -- %
disp('bootstrap...')
N = length(frac);
nboot = 1000;
% create random sampling
for b = 1:nboot
    r = rand(N,1);
    r = ceil(r*N/1);
    frac_boot_b(b) = nanmean(frac(r));
end

% mean of random sampling
frac_boot_dot = nanmean(frac_boot_b);

% Compute 95% confidence interval
frac_sort = sort(frac_boot_b);

CI_2p5 = round(2.5/100*nboot);
CI_97p5 = round(97.5/100*nboot);

frac_2p5 = frac_sort(CI_2p5);
frac_97p5 = frac_sort(CI_97p5);

frac_ave = frac_boot_dot;


%%%%%%%%%%%%%%%%%%%%%%
% - within the CIL - %
%%%%%%%%%%%%%%%%%%%%%%

disp({'fraction...'})
[frac_ave frac_2p5 frac_97p5]

