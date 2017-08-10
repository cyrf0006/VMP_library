function y = batchSpectrum(k, chi, kB)

%This function (m-file) takes a given kB values, the observed spectral
%data, the observed k values, and returns a theoretical spectrum based on them.
%The observed spectral data is included to ensure that the theoretical spectral
%values that are calculated are done at the same intervals as the original,
%making it possible (later) to directly compare the theoretical and observed spectra.
%Batch_spec calls no other m-files or functions.
%
%INPUTS - k: 
%       - chi: 
%       - kB: 
%OUTPUTS -y: a n by 1 vector containing the theoreical spectral values, where
%					x is the length of k_obs or S_obs.


Dt=.00000014;
p  = 0.2316419;
b1=  0.319381530;
b2 =-0.356563782;
b3=  1.781477937;
b4= -1.821255978;
b5=  1.330274429 ;


Q=[];
Z=[];
t=[];  
grif=[];

alpha=[];  
alpha=sqrt(2.*(3.2))*(k./kB);  
t = 1.0./(1.0+p.*alpha);
Z = (1.0./sqrt(2*pi)).*exp(-alpha.*alpha/2.0);
Q = Z .* (((((b5.*t + b4).*t + b3).*t + b2).*t + b1).*t);

grif=alpha .* (exp(-alpha.*alpha/2.0) - alpha .* sqrt(2*pi) .* Q);
y= sqrt(3.2/2).*((chi.*grif)./(kB.*Dt));


      
      
      