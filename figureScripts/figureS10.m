%script to generate Figure S10
%CMJ 20160804
clear variables

figure('units','normalized','outerposition',[0 0 1 1])

M=100;

parameters=PduParams_MCP;
subplot(1,2,1)
sweep_params_ConstantMCP('kcA',-5,5,1,0.1,M,parameters);

parameters=PduParams_MCP;
subplot(1,2,2)
sweep_params_ConstantMCP('kcA',-5,5,1,10,M,parameters);