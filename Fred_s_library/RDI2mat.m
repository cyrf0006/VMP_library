function RDI2mat (fname_in, fname_out)

numav = 1; % No averaging

ADCP = rdradcp(fname_in,numav);

save(fname_out,'ADCP');
