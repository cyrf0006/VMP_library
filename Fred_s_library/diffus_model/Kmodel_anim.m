function Kmodel_anim(Tfname1, Nfname, varargin)
% usage ex:  Kmodel_anim('T_diffus_daily_kobs.dat',
% 'N_diffus_daily.dat', 'T_diffus_daily_5e-5.dat',
% 'T_diffus_daily_tobs.dat')
% 
% T_diffus_daily_tobs.dat est généré par clim2daily dans ctdiml.
%
    
    

if isempty(varargin{1}==1)    
    Tmat1 = load(Tfname1);
    N = load(Nfname);
    P = 1:size(Tmat1,1);

    for i = 1:length(N)
    
    
        plot(Tmat1(:,i),P);
        hold on
        plot(Tmat1(:,1),P, 'r')
        hold off
        set(gca, 'ydir', 'reverse')
        axis([-2 18 0 300])
        ylabel('depth (m)')
        xlabel('T(^{\circ}C)')
        title(datestr(N(i),6))
        pause(0.1)
    end

    
elseif size(varargin,2)==1
    
    Tmat1 = load(Tfname1);
    Tmat2 = load(varargin{1});
    N = load(Nfname);
    P = 1:size(Tmat1,1);

    for i = 1:length(N)
    
        figure(1)
        clf
        set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 15 10])

        subplot(1,2,1)
        plot(Tmat1(:,i),P);
        hold on
        plot(Tmat1(:,1),P, 'r')
        hold off
        set(gca, 'ydir', 'reverse')
        axis([-2 18 0 300])
        ylabel('depth (m)')
        xlabel('T(^{\circ}C)')
        title(datestr(N(i),6))
        text(7,250,'K=K(z)')
        
        subplot(1,2,2)
        plot(Tmat2(:,i),P);
        hold on
        plot(Tmat1(:,1),P, 'r')
        hold off
        set(gca, 'ydir', 'reverse')
        axis([-2 18 0 300])
        set(gca, 'yticklabel', [])
        xlabel('T(^{\circ}C)')
        title(datestr(N(i),6))
        text(7,250,'K=5\times 10^{-5}')
        %pause(0.1)

    % the name of the figure.png file
        if i<10
            figname = sprintf('figure00%d', i);
        else
            if i<100
                figname = sprintf('figure0%d', i);          
            else %profile>100
                figname = sprintf('figure%d', i);
            end
        end
        print('-dpng', '-r300', figname)
    end

    
   
elseif size(varargin,2)==2
    
    Tmat1 = load(Tfname1);
    Tmat2 = load(varargin{1});
    Tmat3 = load(varargin{2});

    N = load(Nfname);
    P = 1:size(Tmat1,1);
    
    tmax=18;
    tmin=-2;
    cil1 = [tmin, tmax; 50 50];
    cil2 = [tmin, tmax; 150 150];
    y1 = [tmin+0.01 tmax-0.01]';
    y2 = [tmin+0.01 tmax-0.01]'; 
    PB = [50 50]';
    PC = [150 150]';
    PA = [0 0]';
    PD = [300 300]';
    
% $$$     II = find(Tmat3(:,1)<=1);
% $$$     A = min(II); B = max(II);
% $$$     CIL_def = Tmat3(:,1)*0+1;
    
    for i = 1:length(N)-4
    
        figure(1)
        clf
        set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 15 10])

        %% -- Sub 1 -- %%
        subplot(1,3,1)
        plot(Tmat3(:,i),P);
        hold on
        patch([y1; flipdim(y2,1); y1(1)],[PA; flipdim(PB,1); PA(1)],[.8 .8 .8], 'edgecolor', 'none');
        %        patch([y1; flipdim(y2,1); y1(1)],[PC; flipdim(PD,1); PC(1)],[.8 .8 .8], 'edgecolor', 'none'); 
        
        II = find(Tmat3(:,i)<=1);
        A = min(II); B = max(II);
        CIL_def = Tmat3(:,i)*0+1;        
        patch([Tmat3(A:B,i); flipdim(CIL_def(A:B),1); Tmat3(A,i)],[P(A:B)';flipdim(P(A:B)',1); P(A)'],[.4 .4 .7], 'edgecolor', 'none');
        
        plot(Tmat3(:,i),P);
        plot(Tmat3(:,1),P, 'r')
        %        plot(CIL_def(A:B), P(A:B), 'k--', 'LineWidth', 0.25)
        plot(cil1(1,:), cil1(2,:), '--k')
        %        plot(cil2(1,:), cil2(2,:), '--k')
        hold off
        set(gca, 'ydir', 'reverse')
        axis([-2 18 0 300])
        ylabel('depth (m)')
        set(gca, 'yticklabel', [])
        %xlabel('T(^{\circ}C)')
        %title(datestr(N(i),6))
        text(6,250,'T climato', 'BackgroundColor',[1 1 1])
        
        
        %% -- Sub 2 -- %%
        subplot(1,3,2)
        plot(Tmat1(:,i),P);
        hold on
        patch([y1; flipdim(y2,1); y1(1)],[PA; flipdim(PB,1); PA(1)],[.8 .8 .8], 'edgecolor', 'none');
        %        patch([y1; flipdim(y2,1); y1(1)],[PC; flipdim(PD,1); PC(1)],[.8 .8 .8], 'edgecolor', 'none');        

        II = find(Tmat1(:,i)<=1);
        A = min(II); B = max(II);
        CIL_def = Tmat1(:,i)*0+1;        
        patch([Tmat1(A:B,i); flipdim(CIL_def(A:B),1); Tmat1(A,i)],[P(A:B)';flipdim(P(A:B)',1); P(A)'],[.4 .4 .7], 'edgecolor', 'none');        
        
        plot(Tmat1(:,i),P);
        plot(Tmat1(:,1),P, 'r')
        %        plot(CIL_def(A:B), P(A:B), 'k--', 'LineWidth', 0.25)
        plot(cil1(1,:), cil1(2,:), '--k')
        %        plot(cil2(1,:), cil2(2,:), '--k')
        hold off
        set(gca, 'ydir', 'reverse')
        axis([-2 18 0 300])
        xlabel('T(^{\circ}C)')
        title(datestr(N(i),6))
        text(6,250,'K=K(z)', 'BackgroundColor',[1 1 1])
        
        %% -- Sub 3 -- %%
        subplot(1,3,3)
        plot(Tmat2(:,i),P);
        hold on
        patch([y1; flipdim(y2,1); y1(1)],[PA; flipdim(PB,1); PA(1)],[.8 .8 .8], 'edgecolor', 'none');
        %        patch([y1; flipdim(y2,1); y1(1)],[PC; flipdim(PD,1); PC(1)],[.8 .8 .8], 'edgecolor', 'none');        

        II = find(Tmat2(:,i)<=1);
        A = min(II); B = max(II);
        CIL_def = Tmat2(:,i)*0+1;        
        patch([Tmat2(A:B,i); flipdim(CIL_def(A:B),1); Tmat2(A,i)],[P(A:B)';flipdim(P(A:B)',1); P(A)'],[.4 .4 .7], 'edgecolor', 'none');           
        
        plot(Tmat2(:,i),P);
        plot(Tmat1(:,1),P, 'r')
        %        plot(CIL_def(A:B), P(A:B), 'k--', 'LineWidth', 0.25)
        plot(cil1(1,:), cil1(2,:), '--k')
        %        plot(cil2(1,:), cil2(2,:), '--k')
        hold off
        set(gca, 'ydir', 'reverse')
        axis([-2 18 0 300])
        set(gca, 'yticklabel', [])
        %xlabel('T(^{\circ}C)')
        %title(datestr(N(i),6))
        text(6,250,'K=5\times 10^{-5}', 'BackgroundColor',[1 1 1])
        %pause(0.1)

        % the name of the figure.png file
        if i<10
            figname = sprintf('figure100%d', i);
        else
            if i<100
                figname = sprintf('figure10%d', i);          
            else %profile>100
                figname = sprintf('figure1%d', i);
            end
        end
        print('-dpng', '-r300', figname)
    end
   
    
    
  
else
    disp('check input files!')
    return
end
