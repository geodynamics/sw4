%%
% by Ke Chen, July 2, 2018
%%
clc;
clear;
close all;
ngrids = 3;
rhop = 265.0000;
mup = 1.5746e9;
lambdap = 2.5335e9;
fileID = fopen('onep-pr10.bin','w');
fwrite(fileID,ngrids,'int');
fwrite(fileID,rhop,'double');
fwrite(fileID,mup,'double');
fwrite(fileID,lambdap,'double');
fclose(fileID);