% Figure for Fred's thesis Appendix 1.
% Largely inspired from Wolk2002, Fig.13
% (called in appendix_dewey_opti.m)

epsList = [1e-15 1e-10 1e-9 1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e0];
expo = log10(epsList);
nasfft = 129;
viscosity = 1e-6;

ks2 = (epsList./viscosity^3).^0.25/(2*pi); % kolmog. [cpm] 
ks = (epsList./viscosity^3).^0.25; % [rad/m]

figure(1)
clf
set(gcf,'PaperUnits','centimeters', 'paperposition', [0 0 30 20])
subplot(1,3,[1 2])
for i = 1:length(epsList)
    [uphi,uk]=nasmyth(epsList(i),viscosity,nasfft); % [cpm]
    
% $$$     I = find(uk<=ks2(i));
% $$$     loglog(uk(I), uphi(I), 'k', 'linewidth', 2)
    if i>1 & i< length(epsList)
        loglog(uk, uphi, 'k', 'linewidth', 2)
        if i == 2
            hold on
        end
 
        %    text(uk(1), uphi(1), num2str(log10(epsList(i))), 'horizontalAlignment', ...
        if i==length(epsList)-1
            text(uk(1), uphi(1), '\epsilon=10^{-3}', 'horizontalAlignment', ...
                 'right', 'verticalAlignment', 'bottom')
        else       
            text(uk(1), uphi(1), sprintf('10^{%d}', expo(i)), 'horizontalAlignment', ...
                 'right', 'verticalAlignment', 'bottom')
        end
    end
    
    % 90% variance
    [Y, I] = min(abs(uk-.5*ks2(i))); % 0.5Kn = 90% variance
    kCut(i) = uk(I); 
    psdCut(i) = uphi(I);
        
    [Y, I] = min(abs(uk-ks2(i))); 
    kKol(i) = uk(I); 
    psdKol(i) = uphi(I); 
    
    [Y, I] = min(abs(uk-.1*2*pi*ks2(i))); 
    kCut1(i) = uk(I); 
    psdCut1(i) = uphi(I);
    
    [Y, I] = min(abs(uk-.001*2*pi*ks2(i))); 
    kCut2(i) = uk(I); 
    psdCut2(i) = uphi(I);
    
    [Y, I] = min(abs(uk-.125*ks2(i))); 
    kCrete(i) = uk(I); 
    psdCrete(i) = uphi(I);
    
    kEnd(i) = uk(end); 
    psdEnd(i) = uphi(end);
end
hold off


ylabel('\phi (s^{-2} cpm^{-1})')
xlabel('k (cpm)')
xlim([4e-3 2e3])
ylim([1e-9 1e0])

hold on
[uphi,uk]=nasmyth(2e-3,viscosity,nasfft);
I = find(uk>=1e0 & uk<=1e2);
loglog(uk(I), uphi(I), '--k', 'linewidth', 2)
text(4, 5e-1, '~k^{1/3}')
% $$$ plot(kCut, psdCut, '--r')
% $$$ plot(kKol, psdKol, '--r')


patch([kCut fliplr(kKol) kCut(1)], [psdCut fliplr(psdKol) psdCut(1)], ...
      [1 1 1]*.7)
patch([kEnd fliplr(kKol) kEnd(1)], [psdEnd fliplr(psdKol) psdEnd(1)], ...
      [1 1 1]*.4)

for i = 2:length(epsList)-1
    [uphi,uk]=nasmyth(epsList(i),viscosity,nasfft);
    loglog(uk, uphi, 'k', 'linewidth', 2)
% $$$ 
% $$$     I = find(uk<=ks2(i));
% $$$     loglog(uk(I), uphi(I), 'k', 'linewidth', 2)
end

plot(kCut1, psdCut1, '--m', 'linewidth', 2)
plot(kCut2, psdCut2, '--m', 'linewidth', 2)
%plot(kCrete, psdCrete, '--r', 'linewidth', 2)
%plot(kCutCut, psdCutCut, '--r', 'linewidth', 2)

% Add shear spectrum
zzVec = [6.5, 54.5, 72.5];
zz = zzVec(1);
I = find(zs>zz-zbin/2 & zs<zz+zbin/2); % I = Imicrosctruc
Ifine = find(z>zz-zbin/2  & z<zz+zbin/2); % Ifinesctruc

II = find(isnan(shear(I))==1);
shear(I(II))=0; % pad NaNs with zeros


[ps1, f1, AA, UU, UA] = clean_shear_spec_fc(A(I,:), shear(I), 512); % Clean Power Spectrum 
Wm = mean(W(Ifine)); % mean falling speed for this bin 
k1 = f1/Wm;
ps10 = ps1;
ps1 = respcorr(ps1, k1);
plot(k1, ps1, 'color', [1 0 0]*.7, 'linewidth', 3, 'lineStyle','--')
I = find(k1>=2.0917 & k1<=168.553);
plot(k1(I), ps1(I), 'color', [1 0 0]*.7, 'linewidth', 3)
%plot(k1, ps10, 'color', [1 0 0]*.7, 'linewidth', 2, 'lineStyle','--')
text(2.05, 1.2e-2, '\epsilon=3.0x10^{-5}', 'verticalAlignment', 'bottom', 'fontweight', 'bold','color', [1 0 0]*.7)
text(2.05, 1.2e-2, '(20% lost)', 'verticalAlignment', 'top','fontweight', 'bold','color', [1 0 0]*.7)
% Nasymth fit
[uphi,uk]=nasmyth(3.9e-5,viscosity,nasfft);
plot(uk, uphi, 'color', [1 0 0]*.7, 'linewidth', 1, 'linestyle', '--')



zz = zzVec(2);
I = find(zs>zz-zbin/2 & zs<zz+zbin/2); % I = Imicrosctruc
Ifine = find(z>zz-zbin/2  & z<zz+zbin/2); % Ifinesctruc

II = find(isnan(shear(I))==1);
shear(I(II))=0; % pad NaNs with zeros


[ps1, f1, AA, UU, UA] = clean_shear_spec_fc(A(I,:), shear(I), 512); % Clean Power Spectrum 
ps10 = ps1;
Wm = mean(W(Ifine)); % mean falling speed for this bin 
k1 = f1/Wm;
ps1 = respcorr(ps1, k1);
plot(k1, ps1, 'color', [0 1 0]*.7, 'linewidth', 3, 'lineStyle','--')
I = find(k1>=2.4039 & k1<=6.5828);
plot(k1(I), ps1(I), 'color', [0 1 0]*.7, 'linewidth', 3)
%plot(k1, ps10, 'color', [0 1 0]*.7, 'linewidth', 2, 'lineStyle','--')
text(.9, 5e-7, '\epsilon=1.0x10^{-9}', 'verticalAlignment', 'bottom', 'fontweight', 'bold','color', [0 1 0]*.7)
text(.9, 5e-7, '(75% lost)', 'verticalAlignment', 'top','fontweight', 'bold','color', [0 1 0]*.7)
% Nasymth fit
[uphi,uk]=nasmyth(1e-9,viscosity,nasfft);
plot(uk, uphi, 'color', [0 1 0]*.7, 'linewidth', 1, 'linestyle', '--')

zz = zzVec(3);
I = find(zs>zz-zbin/2 & zs<zz+zbin/2); % I = Imicrosctruc
Ifine = find(z>zz-zbin/2  & z<zz+zbin/2); % Ifinesctruc

II = find(isnan(shear(I))==1);
shear(I(II))=0; % pad NaNs with zeros

[ps1, f1, AA, UU, UA] = clean_shear_spec_fc(A(I,:), shear(I), 512); % Clean Power Spectrum 
ps10 = ps1;
Wm = mean(W(Ifine)); % mean falling speed for this bin 
k1 = f1/Wm;
ps1 = respcorr(ps1, k1);
plot(k1, ps1, 'color', [0 0 1]*.7, 'linewidth', 3, 'lineStyle','--')
I = find(k1>=1.4917 & k1<=38.0136);
plot(k1(I), ps1(I), 'color', [0 0 1]*.7, 'linewidth', 3)
%plot(k1, ps10, 'color', [0 0 1]*.7, 'linewidth', 2, 'lineStyle','--')
text(2, 1.2e-4, '\epsilon=1.7x10^{-7}', 'verticalAlignment', 'bottom', 'fontweight', 'bold','color', [0 0 1]*.7)
text(2, 1.2e-4, '(35% lost)', 'verticalAlignment', 'top','fontweight', 'bold','color', [0 0 1]*.7)
% Nasymth fit
[uphi,uk]=nasmyth(1.7e-7,viscosity,nasfft);
plot(uk, uphi, 'color', [0 0 1]*.7, 'linewidth', 1, 'linestyle', '--')
set(gca, 'tickdir', 'out')


subplot(133)
plot(shear1, p, 'k')
XLIM = 1;
xlim([-XLIM XLIM])
set(gca, 'ydir', 'reverse')
xlabel('du/dz')
ylabel('Depth (m)')
hold on
% $$$ patch([-XLIM XLIM XLIM -XLIM -XLIM], [zzVec(1)+zbin/2 zzVec(1)+zbin/2 zzVec(1)-zbin/2 ...
% $$$                     zzVec(1)-zbin/2 zzVec(1)+zbin/2], [1 0 0]*.7)
% $$$ patch([-XLIM XLIM XLIM -XLIM -XLIM], [zzVec(2)+zbin/2 zzVec(2)+zbin/2 zzVec(2)-zbin/2 ...
% $$$                     zzVec(2)-zbin/2 zzVec(2)+zbin/2], [0 1 0]*.7)
% $$$ patch([-XLIM XLIM XLIM -XLIM -XLIM], [zzVec(3)+zbin/2 zzVec(3)+zbin/2 zzVec(3)-zbin/2 ...
% $$$                     zzVec(3)-zbin/2 zzVec(3)+zbin/2], [0 0 1]*.7)
% $$$ plot(shear1, p, 'k')
I = find(p>=zzVec(1)-zbin/2 & p<=zzVec(1)+zbin/2);
plot(shear1(I), p(I), 'color', [1 0 0]*.7)
I = find(p>=zzVec(2)-zbin/2 & p<=zzVec(2)+zbin/2);
plot(shear1(I), p(I), 'color', [0 1 0]*.7)
I = find(p>=zzVec(3)-zbin/2 & p<=zzVec(3)+zbin/2);
plot(shear1(I), p(I), 'color', [0 0 1]*.7)
set(gca, 'tickdir', 'out')


print('-dpng', 'appendix_spec_int.png')
set(gcf, 'renderer', 'painters')
print('-depsc2', 'appendix_spec_int.eps')
