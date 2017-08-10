function PSD_Nasmyth = f(eps,nu,k)
  
  % Kolmogorov wavenumber
   ks = (eps/(nu^3))^(1/4);
 
   G = (8.05 * (k./ks).^(1/3) ) ./ (1 + (20 * k./ks).^3.7);
   
   PSD_Nasmyth = ((eps.^(3/4)) ./ (nu^(1/4))) .* G;

end