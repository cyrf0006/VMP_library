function [epsilon,ze] = epsilon_vmp(shear,A,shearFreqHz,T,z,pressFreqHz,varargin)
% function [epsilon,ze] = epsilon_vmp(shear,A,shearFreqHz,T,z,pressFreqHz,varargin)
%
% epsilon_vmp.m--Calculates the dissipation rate [W/kg] by integrating the
% sectioned power spectrum between 1 and 50 Hz. Presently no corrections.
%
% By default, the power spectrum is sectioned in 4-metre bins. Change this
% by specifying a value for the optional input argument "zbin".
%
% The input argument iplt controls the plotting of the results:
%   iplt==0: No plotting done.
%   iplt==[figNum1]: Final summary plot is made in figure window with
%                    handle figNum1.
%   iplt==[figNum1 figNum2 figNum3 spause]: In addition to summary plot, "flash"
%               plots are made during the calculation of epsilon, pausing spause sec.
%   spause: seconds to pause between spectra plots, if negative, do NOT
%   save to file.
%
% Program written by Richard Dewey, June 2005.
%
% Modifications by Kevin Bartlett:
%    2005-06-20: If depth range too small for any bins, return empty
%                variables instead of crashing.

%    2005-06-21: Removed despiking from this routine. Now to be handled by
%                despike.m prior to running epsilon_vmp.m.
%    2006-10-27: New version from Richard. Reformatted (indentation, etc.).
%    2006-10-27: Get rid of global variable pfilename by using
%    parvalpairs.m program to handle input arguments. zbin now included as
%    an optional input argument with the syntax ...'zbin',zbin,..;
%    similarly for pfilename and iplt (see "Syntax", below).
%    2006-10-27: New version makes use of instrument acceleration. Syntax
%    required input of [Ax; Ay; Az]'. For simplicity, changed it to [Ax;
%    Ay; Az] (without the apostrophe to transpose it).
%    2007-10-27: Changed to handle only one shear at a time.
%
% Modifications by F. Cyr:
%   2010-01-08: Modified clipped data for the shear to deal with NaN
%   2010-01-11: Add a break condition if there is flagged data in the bin
%               (the preceeding modification is not needed anymore!
%               I removed it!)
%   2010-01-12: change the length of A to the same as shear when shear is
%               modified at the beginning of the function 
%               see: if ~isempty(shear)...
%   2010-05-27: - remove ffirsti function and move it to
%               clean_shear_spec_rd.m
%               - when no good value for the bin, epsilon=NaN
%               istead of epsilon=0
%  2010-05-28: - Second version of a method to deal with NaNs (see
%                MAIN LOOP). 3 firsts modifs are not needed anymore
%                and have been removed.
%              - Every time this method is applied, the spectrum is
%                plotted to detect errors. Can be removed at the
%                end of plotting section, near the end of main
%                loop.
% 2010-05-31:  - recopy ffirsti function here but left it also in
%                clean_shear_spec_rd.m
%              - keep modif from 2010-05-28, but add an input so
%               that the user can choose to keep or not the spectrum
%               after visual inspection
% 2010-05-31:  - add varargin 'explicit' to explicitly choos eif we
%                keep spectrum. The different between preceeding
%                modifs, is that
%                the question is not asked at each profile but in
%                shear_analysis.m
% 2010-09-??:  - Add a spectrum plot for certain bin. Used for
%                manuscipt 1. Changes are now commented. 
% 2011-02-03:  - Add a condition to skip bin if no falling speed
%                exceeds 0.3m/s (around line 145)
%
%
%
% Syntax:  [epsilon,ze]=epsilon_vmp(shear,A,shearFreqHz,T,z,...
%                         pressFreqHz,<'iplt',iplt>,<'zbin',zbin>,...
%                         <'pfilename',pfilename>)
%
% e.g.,    % Run analyse_vmp.m, open a .p file, crop and make a summary
%          % plot. Then can do the following:
%          vmpData = vmp_data_access('get','cropped'); 
%          shear = vmpData.Sh1.shear;
%          Ax=vmpData.Ax.accel_mpersec_squared;
%          Ay=vmpData.Ay.accel_mpersec_squared;
%          Az=vmpData.Az.accel_mpersec_squared;
%          A=[Ax;Ay;Az]';
%          pressFreqHz = vmpData.P_dP.freqHz;
%          shearFreqHz = vmpData.Sh1.freqHz;
%          smooth = 0.5; thresh = 5; N = 5;
%          [despikedShear,spikeIndex] = despike(4,shear,thresh, smooth,shearFreqHz, N);
%          T = vmpData.SBT1.temperature_C;
%          z = vmpData.P_dP.pressure_db;
%          iplt = 0; zbin=4;
%          [epsilon,ze] = epsilon_vmp(despikedShear,A,shearFreqHz,T,z,...
%             pressFreqHz,'iplt',iplt,'zbin',zbin)

%--------------------------------------------------------------------

% Handle input arguments.
allowNewFields = 0;
isCaseSensitive = 0;
defaultVals.zbin = 4;
defaultVals.iplt = 0;
defaultVals.pfilename = [];
defaultVals.explicit = 0; % explicit on plotting -FC 


parValStruct = parvalpairs(defaultVals,varargin{:},[allowNewFields isCaseSensitive]);
zbin = parValStruct.zbin; % bin size ~4m
pfilename = parValStruct.pfilename;
iplt = parValStruct.iplt;
explicit = parValStruct.explicit; % -FC;

% Transpose acceleration matrix if necessary.
accelMatrixSize = size(A);

if accelMatrixSize(2)~=3 & accelMatrixSize(1)==3
   A = A';
end % if

if size(A,2) ~= 3
   error([mfilename '.m--Acceleration matrix has wrong dimensions.']);
end % if

indx = find(pfilename=='_');

if ~isempty(indx)
   pfilename(indx) = '-';
end

ipplt=0;
% ----- kinematic viscosity & interp1 to remove NaNs ---- %
v=visc(T)/1024;  % kinematic viscosity range of 1-2e-6
I = find(~isnan(v)==1);
XI = 1:length(v); %original length
Y = v(I); % only non NaNs
X = XI(I);
v = interp1(X, Y, XI);
% ------------------------------------------------------- %

zsmth=fix(pressFreqHz*1.5) + mod(pressFreqHz,2) - 1;  % smooth over 1.5 seconds (make odd).
z=flt(z,zsmth);
W = gradient(z) .* pressFreqHz; % the fall speed

% Interpolation of the depth from pressure sensor to the discretisation of
% the shear sensor - FC (is it really needed?)
zs=interp1([1:length(z)],z,[1:length(shear)]*(pressFreqHz/shearFreqHz));

Ws=abs(flt(W,zsmth));  % smooth the fall speed, force positive
indx=find(Ws>0.3);  % only consider data where W>0.2 m/s

if isempty(indx), % no shear data to process if no speed exeeds 0.3m/s
   epsilon = [];
   ze = [];
   return;
end % if

dindx=diff(indx);
if max(dindx)>1, % then there is a gap, only take first block of good fall speeds
    dindx1=min(find(dindx>1));
    indx=indx(1:dindx1);
end
Wm=mean(Ws(indx)); % global average fall speed
z=z(indx);
v=v(indx);
Ws=Ws(indx);
T=T(indx);

% ------ Taking values between z1 & z2 ------- % -FC
% This will remove NaNs at top if there were %
z1=ceil(min(abs(z)));  % top of "falling" (W>20 cm/s) profile
%z1 = 0; -FC (test)
z2=floor(max(abs(z)));  % bottom of falling data

indxs=find(zs>=z1 & zs<=z2);
if ~isempty(shear), %!!! here we reduce the size of shear, we must then reduce the size of A!! - FC
   shear2=shear(indxs);
end
zs2=zs(indxs);
% -------------------------------------------- %


% ------ Param. for FFT (no. pts)------- % -FC
% (zbin/W)*shearFreqHz ~ 3*nfft  power spectral routine needs > 2 nfft of data to overlap 50%
% RKD added longer spectra (2048 and 4096) for large zbin
nnfft=[128 256 512 1024 2048 4096];
nfft=max(nnfft(find(nnfft<((shearFreqHz*zbin)/(3*Wm))))); % nfft is zbin*512/2 - FC 

ipltspec=0; % set to 1 to plot spectra during correction phase
Nfft=min(find([[1:1000]*nfft]>zbin*shearFreqHz/Wm));
N=Nfft*nfft;
N2=fix(N/2);
kc=50; % wavenumber shear probe correction constant
iz=0; %bin counter - FC
pn=0;
% -------------------------------------------- %


% If depth range too small, the following loop will not be entered, and
% output variables will not exist. Return all empty variables if this
% happens kpb 2005-06-20.
zsr=[z1:zbin/2:z2-zbin-1];
if isempty(zsr), % no shear data to process
   epsilon = [];
   ze = [];
   return;
end % if

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---------------- MAIN LOOP !!! ---------------------- %% -FC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for zb=zsr, %loop on each bin, zb is the top of the bin - FC

   zb2=zb+zbin;  %  bottom of bin
   indxs=find(zs >= zb & zs <= zb2); % indices of the shear values in this bin
   
   if ~isempty(indxs) & ~isnan(indxs), % only process if there is data in this depth range (else... next bin! - FC)
       
      
      % ----- compute indexes on which FFT will be applied ----- % 
      indxsmid=fix(mean(indxs)); % mid point of the good data
      iz=iz+1; %count the bin
      ze(iz)=zb+zbin/2; % ze is the mean depth of each bin (returned by the function!) - FC 

      % this is entered if the number of point is too small to
      % perform FFT - FC
      if (indxsmid-N2) < 1,
         indxsmid=N2+1;
      end

      if (indxsmid+N2-1) > length(shear), % should not be entered except very small bins... - FC
         indxsmid=length(shear)-N2;
      end      
      
      %This make sure that the number of pts is always N
      % but may takes few pts out of the range ([zb, zb2])
      indxs=[indxsmid-N2:indxsmid+N2-1];
      % --------------------------------------------------------- %
      
      
      % -------- replacing NaNs (second version, F. cyr) -------- %
      I = find(isnan(shear(indxs))==1);
      iplt = parValStruct.iplt; % reapplies selected values for
                                % plotting
      if ~isempty(I) % only proceed if at least 1 NaN
          if length(I)>round(length(indxs)/2) % more than 50% are NaNs
              disp(['shear not valid in ',num2str(zb+zbin/2),'m bin']);     
              disp(['  -> Returned epsilon = NaN for this bin']);  
              epsilon(iz)=NaN;
              continue
          elseif I(1) == 1% NaNs are at the beginning of bin
              % taking following values, after NaNs
              patching = shear(indxs(I(end))+1:indxs(I(end))+length(I));
              shear(indxs(I))=flipud(patching); %fill NaNs with patching values 
              disp(['Shear extrap. at the begin. of ',num2str(zb+zbin/2),'m bin']);  
              disp([' -> look at spectrum!'])
              iplt=[1 2 3 0.1];
          elseif I(end) == length(indxs) %NaNs are at the end
              % taking preceeding values, before NaNs
              patching = shear(indxs(I(1))-length(I):indxs(I(1))-1);
              shear(indxs(I))=flipud(patching); %fill NaNs with patching values 
              disp(['Shear extrap. at the end of ',num2str(zb+zbin/2),'m bin']);
              disp([' -> look at spectrum!'])
              iplt=[1 2 3 0.1];
          else % NaNs in the middle of the bin
              Y = [shear(indxs(I(1))-1) shear(indxs(I(end))+1)];
              X = [1 length(I)];
              YI = 1:length(I);
              patching = interp1(X, Y, YI);
              shear(indxs(I))=patching; %fill NaNs with patching values 
              disp(['Shear interp. in ',num2str(zb+zbin/2), 'm bin']);
              disp([' -> look at spectrum!'])
              iplt=[1 2 3 0.1];
          end
      end
      
% $$$       %% !!  FC plot for paper on CIL !! %%
% $$$       if zb==42
% $$$         iplt=[1 2 3 0.1]
% $$$       end
      % ----------------------------------------------------------- %

      indxv=find(z >= zb & z <= zb2);
      Wm=mean(Ws(indxv)); % local mean fall speed
      viscosity=mean(v(indxv));  % local viscosity
            
      indxclip=find(abs(shear(indxs))>100); % ?? sort of a despike???? -FC
      if ~isempty(indxclip), 
          shear(indxs(indxclip))=0; 
          disp([' Shear Data Clipped near ',num2str(zb+zbin/2),'m']); 
      end
      
      % ps1 clean spectra has accelerometer signals removed from spectra
      clear ps0 ps1
      
      [ps1,f1]=clean_shear_spec_rd(A(indxs,:),shear(indxs),shearFreqHz,nfft,ipltspec);
      
      ps0=ps1(1,:);
      clear ps1;
      ps1=ps0(:);
      clear ps0;
      [ps0,f0]=ps(shear(indxs),1/shearFreqHz,nfft);  % uncorreced
                                                     % spectra
      % un-comment out the following line if correcting is not working
      if max(ps1)<1e-8 | sum(isnan(ps1))>0 | sum(isinf(ps1))>0,
         ps1=ps0(:);
         f1=f0;
         disp(['Not using clean ps']);
      end % the clean spec is bad

      df=f1(1);
      k1=f1/Wm; % convert frequency to wavenumber
      k1=k1(:)'; % make a wavenumber vector
      k0=f0/Wm; % convert frequency to wavenumber
      k0=k0(:)'; % make a wavenumber vector
      HH=1./(1+(k1/kc).^2);  
      HH=HH(:);% wavenumber correction filter (Macoun&Lueck'04)
        
      ps1=abs(ps1)./HH;
      
      % RKD (0707) a few changes below to correct and improve eps calcualtion
      % first estimate small variance epsilon
      eps0=1e-9; % RKD 07/07
      lost=4.0; % start assuming 75% lost variance
      nasfft=256;
      %
      % RKD 07/07 opened-ip/improved integration limits

      for iloop=1:4, % RKD 07/07 increased the search for ideal integration limits

         % Note: These wavenumber limits KKSL/U are NON-dimensional (k/k_s)
         % RKD 07/07 modified selections below
         if eps0<1e-9,
            kksl=1.5e-2;kksu=3e-2; % lowest epsilon, shink integration band
         elseif eps0>=1e-9 & eps0<1e-8,
            kksl=1.5e-2;kksu=5e-2; % low epsilon, widen integration band
         elseif eps0>=1e-8 & eps0<1e-7,
            kksl=8e-3;kksu=7e-2; % high epsilon, widen integration band
         elseif eps0>=1e-7 & eps0<1e-6
            kksl=3e-3;kksu=9e-2; % higher epsilon, widen integration band
         elseif eps0>=1e-6;
            kksl=1e-3;kksu=1e-1; % highest epsilon, widen integration band
         end

         for ie=1:4,
            
            ks0=(eps0*lost/viscosity^3)^0.25;
            kl0=kksl*ks0;
            ku0=kksu*ks0; % find approximate integration limits 1e-3< k/ks <1e-1
            il=max([1 find(k1 < kl0)]); % index to lower integration limit
            iu=min([find(k1 > ku0) length(k1)]);   % index to upper integration limit
            variance=sum(ps1(il:iu))*df; % integrate power spectrum from kl to ku m^-1
            eps0=7.5*viscosity*variance; % area under spect between limits
            epsn=eps0*lost;
            % RKD 07/07 modified determination of % lost below
            lost=1;
            for in=1:2,
               [uphi,uk]=nasmyth(epsn,viscosity,nasfft);
               uil=max([0 find(uk(1:end-2) < kl0)]) + 1;
               uiu=min([find(uk(3:end) > ku0) length(uk)+1]) - 1;
               uvarb=sum(uphi(uil:uiu));
               uvar=sum(uphi);

               sp=(cumsum(uphi)/uvar);  % WARNING: SOMETIMES DIVIDES BY ZEROS!!! -FC
               if uil==uiu | uil>=nasfft | uiu <=1, keyboard; end % what's happening...
               pll=sp(uil);
               plu=1-sp(uiu);
               lostn=pll + plu; % this is the percentage lost outside integration limits
               lost=1/(1-lostn);
               epsn=eps0*lost;  % updated estimate or epsilon from Nasmyth
            end % in=1:2,
            % RKD 07/07 modified code above
            lost=lost - (lost-2)/2;
            eps0=7.5*viscosity*variance*lost;  % final estimate of Epsilon
         end % for ie=1:4
      end % for iloop=1:4
      
      % check for bad data, skip
      if eps0 < 1e-11 | eps0 > 1e-1 
          disp('Bad FFT, return epsilon=NaN')
          eps0=NaN; 
      end % hold on, bad data, skip

      epsilon(iz)=eps0;
      
      if length(iplt)>2 & ~isnan(eps0),  % if iplt = [# 2 3] then flash-plot each spectra

         pn=pn+1;
         spause=0; % pause between each flash plot

         if length(iplt)>3,
            spause=iplt(4);
         end

         k1=abs(f1/Wm);
         % RKD 07/07 remove k2=abs(f2/Wm); since only 1 shear at a time
         [phi,kn]=nasmyth(epsilon(iz),viscosity,128);
         inp=find(kn > k1(1)/2 & phi > min([min(ps0) min(ps1)]));
         [phil,knl]=nasmyth(epsilon(iz)/10,viscosity,128);
         inl=find(phil==max(phil));
         pnl0=phil(inl(end));
         knl0=knl(inl(end)); 
         [phih,knh]=nasmyth(epsilon(iz)*10,viscosity,128);
         inh=find(phih==max(phih));
         pnh0=phih(inh(end));
         knh0=knh(inh(end));
         f2h=figure(iplt(2));
         if spause<0, set(f2h,'Visible','off'); end
         clf;
         set(gcf,'Units','normalized','Position',[.09 .41 .6 .5]);

         ph=loglog(k1,ps1,'b',kn(inp),phi(inp),'g',k0,ps0,'c');hold on

         loglog([knl0 knh0],[pnl0 pnh0],'-o','Color','g');  % plot a line along peak Nasymth fits
         drawnow
         if isempty(knl0) | isempty(pnl0), keyboard; end
         th(1)=text(knl0*0.8,pnl0*.75,'\times 10^{-1}');
         th(2)=text(knh0*0.85,pnh0*1.45,'\times 10');
         set(th(:),'FontSize',8,'Clipping','On');
         tit=[pfilename,'  Epsilon=',num2str(epsilon(iz)),'  Depth Range ',num2str(zb) ' - ' num2str(zb2),'  W=',num2str(Wm,2)];
         title(tit);
         xlabel('Wavenumber k [m^{-1}]');
         ylabel('Spectral Amplitude');
         smin=min(ps1);smax=max(ps1);
         loglog([k1(il) k1(il)],[smin smax],'r');
         text(k1(il)*0.5,smax*0.8,[num2str(pll*100,3),'%']);
         loglog([k1(iu) k1(iu)],[smin smax],'r');
         text(k1(iu)*1.3,smax*0.8,[num2str(plu*100,3),'%']);
  
         drawnow;
 
% $$$          %% -- FCyr extra spectrum plot for paper on CIL -- %%
% $$$          if zb==42
% $$$          figure(4)
% $$$          clf
% $$$          set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 ...
% $$$                              8 8])
% $$$          ph=loglog(kn(inp),phi(inp),'k--',k0,ps0,'k');hold on
% $$$ 
% $$$          
% $$$          drawnow
% $$$          disp('\epsilon')
% $$$          epsilon(iz)
% $$$          disp('Depth Range') 
% $$$          [zb zb2]
% $$$          
% $$$          xlabel(' k (m^{-1})');
% $$$          ylabel('\Phi (s^{-2} cpm^{-1})');
% $$$          
% $$$          %print('-deps2', 'spectrum_p7_42-43.eps')
% $$$          print('-dpng', '-r300', 'spectrum_p7_42-43.png')
% $$$          end
% $$$           %% -- END of FCyr extra spectrum plot for paper on CIL -- %%
         
        % RKD 07/07 modified plotting to files to handle two shear probes
         if spause~=0 & ~isempty(pfilename),
            pltfn=[pfilename(1:end-2),'-Spectra','.ps'];
            if pn==1 & exist(pltfn,'file')==0, % since we're now doing 1 shear at a time
               disp(pltfn);
               print('-dpsc2',pltfn);
            else
               print('-dpsc2','-append',pltfn);
            end;
         end

         %
         f3h=figure(iplt(3));
         if spause<0, set(f3h,'Visible','off'); end
         clf;
         set(gcf,'Units','normalized','Position',[0.7 .41 .2 .5]);
         ppts=-[1:length(shear(indxs))];
         %plot(shear(indxs)-0.5,ppts,shear2(indxs)+0.5,ppts,A(indxs,1),ppts,A(indxs,2),ppts);
         plot(shear(indxs)-0.5,ppts,A(indxs,1),ppts,A(indxs,2),ppts);
         %title('Shear 1&2, Accel 1&2');
         title('Shear, Accel 1&2');
         drawnow;
         if spause>0, pause(spause); end
         
         % Explicit asking if spectrum OK after visual inspection
         if explicit==1
             answer=0;
             while answer==0
                 R2 = input(' Do you want to keep this bin? (y/n) ', 's');
                 if strcmp('y',R2)==1
                     disp(' Spectrum accepted!')
                     answer=1;
                 elseif strcmp('n',R2)==1  
                     disp(' Spectrum rejected!')
                     answer=1;
                     epsilon(iz)=NaN;
                 else
                     disp(' Bad answer, pleaser type (y or n)')
                 end
             end
         end
         
      end % if plotting

   end % if index is empty
   
end % ------ END OF MAIN LOOP ---- %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------------- Plotting if needed -------------------- %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% replace default value for iplt -FC
iplt = parValStruct.iplt;

if iplt(1)~=0,
   f1h=figure(iplt(1));
   if spause<0, set(f1h,'Visible','off'); end
   clf;
   set(gcf,'Units','normalized','Position',[.09 .41 .6 .5]);
   subplot(1,3,1);
   [sz,iT]=sort(-z);
   plot(T(iT),z(iT));
   grid;
   YLim=[z1 z2];
   set(gca,'YLim',YLim,'YDir','reverse');
   xlabel('Temperature [C]');
   ylabel('Depth [m]');
   
   subplot(1,3,2);
   plot(decimate(shear,10),decimate(zs,10));
   axis([-10 10 z1 z2]);
   set(gca,'YDir','reverse');
   xlabel('Shear [s^{-1}]');
   ylabel('Depth [m]');

   h=subplot(1,3,3);
   semilogx(epsilon,ze);
   set(gca,'YDir','reverse','YLim',[z1 z2]);
   set(h,'LineWidth',1.0);
   grid;
   XTickA=[1e-11 1e-10 1e-9 1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2];
   XTickLabelA=['-11';'-10';' -9';' -8';' -7';' -6';' -5';' -4';' -3';' -2'];
   emin=10^floor(min(log10(epsilon)));
   emax=10^ceil(max(log10(epsilon)));
   ind1=max([find(XTickA<=emin) 1]);
   ind2=min([find(XTickA>=emax) 10]);
   axis([XTickA(ind1) XTickA(ind2) z1 z2]);
   set(gca,'XTick',XTickA(ind1:ind2),'XTickLabel',XTickLabelA(ind1:ind2,:));
   xlabel('Log_{10} Dissipation Rate [W/kg]');
   ylabel('Depth [m]');
   title('Log(Epsilon) [W/kg]');
   %suptitle(pfilename);
   drawnow;
   if spause<0, set(f1h,'Visible','off'); end

% RKD 07/07 modified to handle two separate call for two shear probes
   if spause~=0 & ~isempty(pfilename),
      ipplt=ipplt+1;
      indot=find(pfilename=='.');
      pltfn2=[pfilename(1:indot-1),'-TSE','.ps'];
      if ipplt==1 & exist(pltfn2,'file')==0,
         disp(pltfn2);
         print('-dpsc2',pltfn2);
      else
         print('-dpsc2','-append',pltfn2);
      end
   end
end % plotting

% ------------------------ End of function ------------------------ %


%--------------------------------------------------------------------
function v=visc(T)
% function v=visc(T)
% calculates the kinematic viscosity according to Dewey's thesis
% units are m^2/s
% RKD 4/97
% v=(1.87e-3)*(1.0+0.0323*T+0.00023*T.^2).^(-1);
v=(1.87e-3)*(1.0+0.0323*T+0.00023*T.^2).^(-1);
% fini

%--------------------------------------------------------------------
function xf=flt(x,fpts)
%
% function xf=flt(x,fpts) to low pass filter a time series using filtfilt
%   passed both ways.
%   Since filtfilt passes both ways already, this ties both ends
%   of the filtered time series to the original ends.
%	x = array
%	fpts = number of points for filter length (odd) i.e. 3,5,7
%	xf = filtered array
% RKD 12/94

% make a mask for the data that are NaN's
mask=ones(size(x));
idx=isnan(x);mask(idx)=idx(find(idx)).*NaN;
idx=isinf(x);mask(idx)=idx(find(idx)).*NaN;
x=x.*mask;
% chop off leading and trailing NaN's but keep track so they can be put back
%    later.
npts=length(x);
idx1=min(find(~isnan(x))); idx2=max(find(~isnan(x)));
lnan=[];
tnan=[];

if idx1>1
   lnan=x(1:idx1-1);
end
if idx2<npts
   tnan=x(idx2+1:npts);
end
x=x(idx1:idx2);

% clean up any interior gaps by linearly interpolating over them
[x,tcln]=cleanx(x,1:length(x));

% calculate filter weights
Wn=fix(1.5*fpts);
b=fir1(Wn,(1/(fpts*2)));
m=length(x);
% set weights for forward and reverse sum
winc=1./(m-1);
wt1= 1:-winc :0 ;
wt2= 0: winc :1 ;

% calculate forward filter time series, save
xf1=zeros(size(x));
xf1=filtfilt(b,1,x);

% flip (invert) time series, filter, and unflip
xi=zeros(size(x));
xf2=zeros(size(x));

% flip time series but check for horizontal or vertical matrix
[r,c]=size(x);

if r>c
   xi=flipud(x);
else
   xi=fliplr(x);
end

xi=filtfilt(b,1,xi);

if r>c
   xf2=flipud(xi);
else
   xf2=fliplr(xi);
end

xf=zeros(size(x));
if r>c
   xf=xf1.*wt2' + xf2.*wt1';
else
   xf=xf1.*wt2 + xf2.*wt1;
end

% reattach leading and trailing Nan's and apply mask for interior gaps
if r>c
   xf=[lnan
      xf
      tnan];
else
   xf=[lnan xf tnan];
end

xf=xf.*mask;
% fini

%--------------------------------------------------------------------
function [pspec,f,coh]=ps(x,dt,nfft,y)
%  [pspec,f,coh]=ps(x,dt,nfft,y) calculates power (cross) spectra
%  for the time series x (and optionally y), vectors only(!)
%  divided into 50% overlap segments. f are the frequencies.
%  dt passed for correct scaling and setting frequencies
%     in time units (seconds, hours, days: freq = Hz, cph, cpd)
% Optional input parameters are:
%     y: (then output is cross spectra)
%     nfft: length of fft segments,
%           or pass nfft=1 if you want to routine to estimate nfft
% Optional output parameters are:
%     coh = |Pxy|./sqrt(Pxx.*Pyy);  % the corherency funtion
%
%  Equal to MATLAB's psd*2*dt=psd/f_N  [units^s/freq].
%  Version 1.0 RKD 5/96, 2.0 2/97, 3.0 4/97
x=x(:)';
iy=0;
pspec=[];
if exist('y'),
   iy=1;
   y=y(:)';
end  % only do cross stuff if y exists.

if nargout==4 & (length(y)/1.9) < nfft,
   disp(' You must segment dataset in order to get none-one coherencies.');
end

% chunk and interpolate NaN values in both real and imaginary ts first.
if sum(isnan(x)) > 0,
   x=cleanx(x);
end  % cleanx is an RKD routine

if iy==1,
   if sum(isnan(y)) > 0,
      y=cleanx(y);
   end;
end

% now find new time series length
T=length(x);
if iy == 1,
   Ty=length(y);
   if Ty ~= T,
      disp(['Cleaned X and Y vectors not same length.']);
      if Ty<T,
          dT=T-Tt;
          y([1:dT]+Ty)=mean(y);
      else
          dT=Ty-T;
          x([1:dT]+T)=mean(x);
      end
   end
end

% calculate various spectral parameters
if exist('nfft'),  % check to see if user has set fft length
   M=nfft;         % yes...
else
   M=1;            % no... so do next if loop
end

if M == 1 | T < 3*M,

   n2=2.^(5:13);   % maximum auto spectral length is 8192
   M=n2(max(find(n2<fix(T/3)))); % choose the longest M with at least 3 segments

   if isempty(M),  % if that didn't work, choose less segments
      M=n2(max(find(n2<fix(T)))); % choose the longest M with 1 segments
   end

end

nps=fix(2*T/M) - 1;  % number of 50% overlapping sections/spectra

if nps < 1,
   nps = 1;
end

if M > length(x) | nps == 1,
   M = length(x);
end

if rem(M,2),
   n2=(M+1)/2;    % M is odd
else
   n2=M/2;        % M is even
end

window=hanning(M)';
I=(sum(window.^2));
W=M/I;               % variance lost by windowing
pspec=(0+sqrt(-1)*0)*ones(1,n2); % initialize vector for +f only

if iy==1,
   coh=pspec;
   PXX=pspec;
   PYY=pspec;
end % initialize coherency vectors
%

for jj=1:nps   % loop for segment, fft, sum and average

   nst=1+(jj-1)*M/2;
   nen=nst+M-1;
   indx=fix([nst:nen]);
   X=fft((window.*(detrend(x(indx)))),M); % this is complex

   if iy==1,                    % this loop only if cross spectra/coh
      Y=fft((window.*(detrend(y(indx)))),M); % this is complex
      PS=Y.*conj(X);      % calculate the cross spectra [see B&P (9.3)]
   else
      PS=X.*conj(X);      % or just the auto-spectra, this is not complex
   end

   if iy == 1,  % in addition, calculate Pxx and Pyy for coherency
      Pxx=X.*conj(X);
      Pyy=Y.*conj(Y);
   end

   pspec=pspec + PS(2:n2+1);   % sum spectra not including f=0

   if iy == 1,
      PXX=PXX + Pxx(2:n2+1);   % Must sum the spectra before calculating
      PYY=PYY + Pyy(2:n2+1);   % the coherency function, otherwise coh=1.0
   end

end

%
if iy==1,
   coh = abs(pspec)./sqrt(PXX.*PYY);   % calculate coherency function
end
%
pspec=pspec*(2*dt*W/(M*nps));      % scale power spectra
% if you want c.i. pspec*[plow,phi]=chisqp(95,nps) => a Dewey routine
f=(1:length(pspec))*(1/(M*dt));    % set frequency vector
% fini

%--------------------------------------------------------------------
function [phi,k] = nasmyth(epsilon,nu,N)
% [phi,k] = nasmyth(epsilon,nu,N)
%
% NASMYTH generates a Nasmyth universal shear spectrum
%         for a specified dissipation (epsilon [W/kg]) and
%         viscosity (nu [m^2/s]);
%         phi is the dimensional shear spectrum with N points
%         and k is the wavenumber in cycles per meter.
%         The defult values for nu and N are N = 1000 and
%         nu = 1.0e-6;
%
%  There are 3 forms of the function
%  [phi,k] = nasmyth(epsilon,nu,N)
%  [phi,k] = nasmyth(epsilon,nu)
%  [phi,k] = nasmyth(epsilon)
%
% 1st Version By D. Huang, G2's formula fited by R. Lueck, CEOR, UVic
%
% Feb. 05, 1996
if nargin == 0
   error('No specific dissipation can be used')
elseif nargin == 1
   nu = 1.0e-6;
   N = 1000;
elseif nargin == 2
   N = 1000;
end

x = logspace(-4,0,N);  % x = k / ks;
G2 = 8.05 * x.^(1/3) ./ (1 + (20 * x).^(3.7));
ks = (epsilon ./ nu.^3) .^(1/4);
k = ks * x;
phi = epsilon.^(3/4) * nu.^(-1/4) * G2;
% fini

%--------------------------------------------------------------------
function [a]=avg(x)
% function [mean]=avg(x);
% Calculate the mean of x columns ignoring NaNs
% RKD 02/02
x=x(:)';
[n,m]=size(x);
a=ones(n,1)*NaN;

for i=1:n,
   in=~isnan(x(i,:));
   if sum(in)>0,
      a(i)=mean(x(i,in));
   end
end
%

%--------------------------------------------------------------------
function [tsout,tout]=cleanx(tsin,tin)
%
%  Function [tsout,tout]=cleanx(tsin,tin)
%   or
%  Function [tsout]=cleanx(tsin)
%  Fills in the NaN's in a data set by using linear interpolation,
%
%  4/7/94 RKD - chunk end NaN values first.
% Chunk off NaN values.
if nargin == 1, tin=(1:length(tsin)); end
n=length(tsin);

ist=min(find(isnan(tsin)==0));
iend=max(find(isnan(tsin)==0));
if isempty(ist), disp('No good data in this time series (cleanx.m)'); return; end;
if ist ~= 1 | iend ~= n
   newn=iend-ist+1;
   if newn <= 1, return, end;
   tsnew(1:newn)=tsin(ist:iend);
   tnew(1:newn)=tin(ist:iend);
   tsin=tsnew;
   tin=tnew;
end
%
nn  = 0;                                   % counter for NaNs filled
lt = length(tsin);                         % length of input array
%if lt ~= n
%   fprintf(1,'%5d NaNs truncated from time series\n',(n-lt));
%end
nan1 = isnan(tsin);
nnan = sum(nan1);                          % # of NaN's to be replaced
%
% send info to screen
%
% loop through data filling in NaN's
for ii = 2:1:lt-1
   if nan1(ii) == 1
      
      jj = (ii-1) + ffirsti(tsin(ii:lt)); % find next good point
      kk = jj - ii;                       % number of points to fill
      delx = (tsin(jj)-tsin(ii-1))/(kk+1);
      for i=1:kk
          tsin(ii-1+i)=tsin(ii-1)+(i*delx);
      end
      nan1 = isnan(tsin);     % recalculate the NaN array
      nn = nn+kk;             % increment summary counter
   end
end
%if nn > 0
%  fprintf(1,'%5d NaNs were filled\n',nn);
%end
tsout=tsin;
tout=tin;
% fini

function [kk] = ffirsti(xdata)
%
% Function FFIRSTI.M
% Finds the index of the first good data point in xdata.
%
%       Usage:  index = ffirsti( xdata )
%       Where:  index = output index
%               xdata = input data column

%   MDM  29 Dec 93   Original (FFIRST)
%   MDM  18 Jan 94   Modified FFirst to output index instead of time.
%   JTG  12 Dec 94   Changed to use FIND instead of a loop
%                    Send message if no good data exists and returns 0
%

kk=min(find(~isnan(xdata)));
if isempty(kk)
     kk=0;
     disp('No good data in array')
     return
else
%      fprintf(1,'First good index is %6d \n',kk);
end     


