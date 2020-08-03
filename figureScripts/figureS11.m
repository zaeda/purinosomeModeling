%script to generate Figure S11
%CMJ 20160804
clear variables

M=50;

figure('units','normalized','outerposition',[0 0 1 1])

sweepPerformance('kcA',-5,2,1,3,M,3,1);

sweepPerformance('kcA',-5,2,0,3,M,3,2);

sweepPerformance('kcP',-5,2,0,3,M,3,3);