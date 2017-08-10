clear

tic
% K range
K_i = 2e-5:5e-6:1e-4;
size(K_i)
% Comparation time

%for mm = 7:9
%mm=10;
%ncompare = datenum(999, mm, 15);

ncompare = datenum(999, 5:11, 15);

initT = load('T_bootclim_05.dat');
initS = load('S_bootclim_05.dat');

dlmwrite('initT.dat', [initT(:,1), initT(:,2)+initT(:,3)],'delimiter',' ','precision',6)
dlmwrite('initS.dat', [initS(:,1), initS(:,2)+initS(:,3)],'delimiter',' ','precision',6)


for i = 1:length(K_i)
    
    diffus_model_noflux('initT.dat', 'initS.dat' ...
                        ,datenum(999,5 ,15), K_i(i), 'constant_k', 'daily_forcing.dat',100, 3600*24, 3600*24*30*7)
    
    
    % -- load model results -- %
    T=load('T_diffus_daily.dat');
    S=load('S_diffus_daily.dat');
    N=load('N_diffus_daily.dat');
    P = [1:300]';
    
    
    % cummulated misfit
    misfitT = 0;
    misfitS = 0;
    for j = 1:length(ncompare)
        [Y, I] = min(abs(N-ncompare(j)));
        T_mod = T(:,I); % profiles for misfit
        S_mod = S(:,I);

        
        % -- load observations -- %
        month = str2num(datestr(ncompare(j),5));    
        
        if month <10
            Toutname  = sprintf('T_bootclim_0%d.dat', month);
            Soutname  = sprintf('S_bootclim_0%d.dat', month);
        else
            Toutname  = sprintf('T_bootclim_%d.dat', month);
            Soutname  = sprintf('S_bootclim_%d.dat', month);
        end
        
        TT = load(Toutname);
        SS = load(Soutname);
        T_obs = TT(:,2)- TT(:,3);
        S_obs = SS(:,2)- SS(:,3);

        
        % keyboard
        % -- find misfit -- %
        I = find(P>50 & P<150);
        misfitT = misfitT + sum((T_obs(I)-T_mod(I)).^2);
        misfitS = misfitS + sum((S_obs(I)-S_mod(I)).^2);

% $$$         if i==1
% $$$             figure(1)
% $$$             clf
% $$$             set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 15 ...
% $$$                                 18])
% $$$             figure(2)
% $$$             clf
% $$$             set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 15 18])
% $$$         end
        
                
    end
    bestK_T(i,:) = [K_i(i), misfitT];
    bestK_S(i,:) = [K_i(i), misfitS];
% $$$     
% $$$     figure(1)
% $$$     subplot(4,5,i)
% $$$     plot(T_obs, P, 'r')
% $$$     set(gca, 'ydir', 'reverse')
% $$$     title(sprintf('K = %d', K_i(i)))
% $$$     xlabel('T')
% $$$     hold on
% $$$     plot(T_mod, P)
% $$$     hold off
% $$$     set(gca, 'ygrid', 'on')
% $$$     set(gca, 'xgrid', 'on')
% $$$ 
% $$$     if i == length(K_i)
% $$$         print('-depsc2', 'differentK_forT.eps')
% $$$     end
% $$$     
% $$$ 
% $$$     figure(2)
% $$$     subplot(4,5,i)
% $$$     plot(S_obs, P, 'r')
% $$$     set(gca, 'ydir', 'reverse')
% $$$     title(sprintf('K = %d', K_i(i)))
% $$$     xlabel('S')
% $$$     hold on
% $$$     plot(S_mod, P)
% $$$     hold off
% $$$     set(gca, 'ygrid', 'on')
% $$$     set(gca, 'xgrid', 'on')
% $$$     if i == length(K_i)
% $$$         print('-depsc2', 'differentK_forS.eps')
% $$$     end
    
    %pause
     
end


figure(1)
clf
semilogx(bestK_T(:,1), bestK_T(:,2), 'k', 'linewidth', 2)
xlabel('K values')
ylabel('misfit')
title('Best K for T')
print('-deps2', 'bestK_forT.eps')

figure(2)
clf
semilogx(bestK_S(:,1), bestK_S(:,2), 'k', 'linewidth', 2)
xlabel('K values')
ylabel('misfit')
title('Best K for S')
print('-deps2', 'bestK_forS.eps')


[Y, It] = min(bestK_T(:,2));
[Y, Is] = min(bestK_S(:,2));

%disp(sprintf('month %d', mm));
[bestK_T(It,1) bestK_S(Is,1)]


dlmwrite('bestK_T_plusmin.dat', bestK_T,'delimiter',' ','precision',6)


% $$$ pause
% $$$ 
% $$$ end

toc