% Script written to execute many times SST_map.m


[thresh1, S1] = SST_map('sstMap_may.list', [47.25 49.5], [-70.5 -67])
[thresh2, S2] = SST_map('sstMap_june.list', [47.25 49.5], [-70.5 -67])
[thresh3, S3] = SST_map('sstMap_july.list', [47.25 49.5], [-70.5 -67])
[thresh4, S4] = SST_map('sstMap_aug.list', [47.25 49.5], [-70.5 -67])
[thresh5, S5] = SST_map('sstMap_sept.list', [47.25 49.5], [-70.5 -67])
[thresh6, S6] = SST_map('sstMap_oct.list', [47.25 49.5], [-70.5 -67])

[thresh1 thresh2 thresh3 thresh4 thresh5 thresh6]'
[S1 S2 S3 S4 S5 S6]'

