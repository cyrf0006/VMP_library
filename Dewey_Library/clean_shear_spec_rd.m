function  [clean_UU, F, AA, UU, UA] = clean_shear_spec_rd(A, U, rate, n_fft,iplot)
%
%function [clean_UU, F, AA, UU, UA] = clean_shear_spec_rd(A, U, rate, n_fft, iplot)
%
% function to remove acceleration contamination from all shear channels
% using all available acceleration channels and other signals related shear contamination.
% Uses full and partial coherency. 
% A = [ax ay az ...] are the accelerations (can be more than 3),
% U = [sv sw ...] are the shear probe signals (can be more than 2).
% The time series for each vector must progress down the columns of the
% matrices.
% n_fft length of fft for calculation of auto- and cross-spectra [i.e. 1024]
% rate is the sampling rate of the time series [Hz].
% If included iplot=1 to plot co-spectra, otherwise do not plot
%
% The outputs are:
% the cleaned up cross-spectrum of shear, clean_UU,
% the acceleration cross-spectra, AA
% the original (dirty) shear cross-spectrum, UU,
% the cross-spectrum of dirty shear and acceleration, UA, and
% the frequency of the cross-spectra, F
% If you only need the cleaned-up cross-spectra of shear use
% clean_UU = clean_shear_spec(A,U,n_fft,rate);
%
% Written by Lou Goodman, circa 2003.
% Dicked with by R. Lueck 2004-05-03.
% Uses RKD's ps routine instead of psd_rolf and csd_rolf
% do some checking
if ((length(A(:)) == 1) | (length(U(:)) == 1))
    error('Acceleration and shear matrices must contain vectors')
end
if (size(A,2) > size(A,1))
    error('the vectors do not seem to go down the columns of the acceleration matrix')
end
if (size(U,2) > size(U,1))
    error('The vectors do not seem to go down the columns of the shear matrix')
end
if (size(A,1)~=size(U,1))
    error('Acceleration and shear matrices must have the same number of rows')
end
%if (size(A,1) <= 2*n_fft)
%    error('Your fft length is too long for the length of the vectors')
%end
if (rate <= 0)
    error('Sampling rate is negative')
end
if nargin<5, iplot=0; end
% end of checking
time=1/rate*[1:length(U(:,1))];

% Now the acceleration cross-spectra
%AA = zeros(size(A,2),size(A,2),n_fft/2 + 1);% pre-allocated the matrix
for k = 1:size(A,2),
    [AA(k,k,:), F] = ps(A(:,k),1/rate,n_fft); % use RKD's ps
%    [AA(k,k,:), F] = psd_rolf(A(:,k),n_fft,rate);
    for m= k+1:size(A,2)
      AA(k,m,:) = ps(A(:,k),1/rate,n_fft,A(:,m));  % RKD
%      AA(k,m,:) = csd_rolf(A(:,k),A(:,m),n_fft,rate);
      AA(m,k,:) = conj(AA(k,m,:));
  end
end

% Now the shear probe cross-spectra
% UU = zeros(size(U,2),size(U,2),n_fft/2 + 1);
for k = 1:size(U,2)
     [UU(k,k,:), F] = ps(U(:,k),1/rate,n_fft); % RKD
%    UU(k,k,:) = psd_rolf(U(:,k),n_fft,rate);
    for m = k+1:size(U,2)
        UU(k,m,:) = ps(U(:,k),1/rate,n_fft,U(:,m)); % RKD
%        UU(k,m,:) = csd_rolf(U(:,k),U(:,m),n_fft,rate);
        UU(m,k,:) = conj(UU(k,m,:));        
    end
end

% Now the shear probe and acceleration cross-spectra
% UA = zeros(size(U,2),size(A,2),n_fft/2 + 1);
for k = 1:size(U,2)
    for m = 1:size(A,2)
        UA(k,m,:) = ps(U(:,k),1/rate,n_fft,A(:,m)); % RKD
%        UA(k,m,:) = csd_rolf(U(:,k),A(:,m),n_fft,rate);
    end
end


corrUU=NaN*size(UU);
for ii = 1:length(F)
   corrUU(:,:,ii)=(UA(:,:,ii)/AA(:,:,ii))*(conj(UA(:,:,ii)).');
   clean_UU (:,:,ii) = UU(:,:,ii) - corrUU(:,:,ii);
end
clean_UU = abs(squeeze(clean_UU));
% setiplot = 1 for diagnostic view of spectral peak reduction
if iplot,
    % figure(4);clf;plot(time,U,time,detrend(A));
    figure(5);clf;orient tall;AA1=squeeze(AA(1,1,:));AA2=squeeze(AA(1,2,:));AA3=squeeze(AA(1,3,:));
    subplot(3,1,1);loglog(F,abs(AA1),F,abs(AA2),F,abs(AA3));title('Accelerometers 1,2,3');
    UU1=squeeze(UU(1,1,:));UA1=squeeze(UA(1,1,:));UA2=squeeze(UA(1,2,:));UA3=squeeze(UA(1,3,:));
    subplot(3,1,2);loglog(F,abs(UU1),F,abs(UA1),F,abs(UA2),F,abs(UA3));title('Shear (b) and Cross Spectra Shear x Accel 1,2,3');
    cUUA=squeeze(corrUU);
    subplot(3,1,3);loglog(F,abs(UU1),'b',F,abs(cUUA),'g',F,clean_UU,'r');title('Shear (b), Correction (g), and Cleaned (r)');
    % pause;
end
% fini
%--------------------------------------------------------------------
function [pspec,f,coh]=ps(x,dt,nfft,y)
%  [pspec,f,coh,nof]=ps(x,dt,nfft,y) calculates power (cross) spectra
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
          dT=T-Ty;
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
ndof=2*T/I; % number of degrees of freedom
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
   coh = abs(pspec)./sqrt(PXX.*PYY);   % calculate coherency function  (WARNING!!!! Divide often by zero -FC)
end
%
pspec=pspec*(2*dt*W/(M*nps));      % scale power spectra
% if you want c.i. pspec*[plow,phi]=chisqp(95,nps) => a Dewey routine
f=(1:length(pspec))*(1/(M*dt));    % set frequency vector
% fini

%--------------- added by F Cyr (from epsilon_vmp.m) ----------------------
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

%--------------- added by F Cyr (from epsilon_vmp.m) ----------------------
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

