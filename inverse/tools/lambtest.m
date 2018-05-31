%
% LAMBTEST
% Evaluate some resolution and runtime numbers for solving Lamb's problem with SW4 and WPP
%
% [p dist scpu wcpu msol]=lambtest()
%
% p: points per wave length
% dist: src-rec distance
% scpu: sw4-cpu times
% wcpu: wpp-cpu times
% msol: max norm of solution on finest grid
%
function [p dist scpu wcpu msol]=lambtest()
% there are 7 different resolutions, here, P = N/2.75, N=10,15,20,30,40,60,80
p=[3.63 5.45 7.27 10.91 14.54 21.81 29.09];
% cpu time computed as CPU*cores
scpu=[947.2 4774.4 13952 70656 225280 1056256 3412992];
wcpu=[384 2048 5888 32000 90880 433152 1297920];

% there were 10 receivers at distances:
dist=[1 2 3 4 5 6 7 8 9 10];
% max norm of solution
msol=[1.0850e-01 7.1686e-02 5.7939e-02 5.0084e-02 4.4756e-02 4.0830e-02 3.7781e-02 3.5333e-02 3.3305e-02 3.1589e-02];
