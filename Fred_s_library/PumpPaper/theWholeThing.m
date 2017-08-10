clf

zbin=10;

cd /home/cyrf0006/PhD/Nitrates/

cd ./Stat16_50km/
scatterNO3('dat_profiles','scatter_stat16', zbin, 500, 2000)

cd ../Stat17_50km/
scatterNO3('dat_profiles','scatter_stat17', zbin, 500, 2000)

cd ../Stat18_50km/
scatterNO3('dat_profiles','scatter_stat18', zbin, 500, 2000)

cd ../Stat19_50km/
scatterNO3('dat_profiles','scatter_stat19', zbin, 500, 2000)

cd ../Stat20_50km/                                    
scatterNO3('dat_profiles','scatter_stat20', zbin, 500, 2000)

cd ../Stat21_50km/
scatterNO3('dat_profiles','scatter_stat21', zbin, 500, 2000)

cd ../Stat22_50km/
scatterNO3('dat_profiles','scatter_stat22', zbin, 500, 2000)

cd ../riki        
scatterNO3('dat_profiles','scatter_stat23', zbin, 500, 2000)

cd ../Stat24_20km/                                    
scatterNO3('dat_profiles','scatter_stat24', zbin, 500, 2000)

cd ../Stat25_20km/                                    
scatterNO3('dat_profiles','scatter_stat25', zbin, 500, 2000)

cd ../transect/
NO3_transect