function viztopo(x,z,H);

[dum imax]=size(x);   % Find size
[dum kmax]=size(z);   % Find size 

ii=0;
for i=1:imax-1
  if H(i+1)==H(i)
    ii=ii+1;
    XV(ii)=x(i);
    ZV(ii)=H(i);
  else
    ii=ii+1;
    XV(ii)=(x(i+1)+x(i))/2;
    ZV(ii)=H(i);
    ii=ii+1;
    XV(ii)=(x(i+1)+x(i))/2;
    ZV(ii)=H(i+1);
  end
end

dH=0.01*max(max(H));
XV(ii+1)=x(imax);
ZV(ii+1)=H(imax);
XV(ii+2)=x(imax);
ZV(ii+2)=max(z);
XV(ii+3)=x(1);
ZV(ii+3)=ZV(ii+2);
XV(ii+4)=XV(1);
ZV(ii+4)=ZV(1);

fill(XV,ZV,[[40 32 33]/255);
%fill(XV,ZV,[0. 0. 0.]);
set(gca,'YDir','reverse');
