%% 单天线，C/A码
clear; clc;
sample_offset = 0*4e6;
acqResults = GPS_L1_CA_acq('E:\GNSS data\data_20190413_191804_ch1.dat', sample_offset, 8000);

%% 双天线，C/A码
clear; clc;
sample_offset = 0*4e6;
acqResults_A = GPS_L1_CA_acq('E:\GNSS data\data_20190414_220821_ch1.dat', sample_offset, 8000);
acqResults_B = GPS_L1_CA_acq('E:\GNSS data\data_20190414_220821_ch2.dat', sample_offset, 8000);

%% 