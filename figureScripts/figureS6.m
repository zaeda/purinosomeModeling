%script to generate Figure S6
%CMJ 20160804
clear variables

M=50;

figure('units','normalized','outerposition',[0 0 1 1])

sweepPerformance('Pout',-5,2,1,2,M,1,1);