function misfit = spectral_fit(eps,nu,K,P)
%function misfit = spectral_fit(cv,K,P)

tiny = 1.e-16;

if eps < 0
    
  misfit=10.d10;
  
else
  fit=PSD_Nasmyth(eps,nu,K);
  %fit=PSD_Panchev_Kesich(eps,nu,K);

  %misfit=sum(((fit-P)./P).^2)
  %misfit=sum(((fit-P)).^2)
  
  % Suggested by Ruddick et. al. The tiny is simply to avoid
  % log(0)
  misfit=sum((log(fit+tiny)-log(P+tiny)).^2);
  
  %misfit=sum((log10(fit+tiny)-log10(P+tiny)).^2);
  %misfit=mean(abs( (P./(fit+tiny))-1));
  %misfit=mad(P./(fit+tiny));
  
end