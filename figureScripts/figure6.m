%script to generate Figure 7
%CMJ 20160804
clear variables

figure('units','normalized','outerposition',[0 0 1 1])

M=100;

parameters=PduParams_MCP;
subplot(2,3,1)
sweep_params_ConstantMCP('kcA',-5,5,1,1,M,parameters);

parameters=PduParams_MCP;
subplot(2,3,2)
sweep_params_ConstantMCP('kcA',-5,5,0,1,M,parameters);

parameters=PduParams_MCP;
subplot(2,3,4)
sweep_params_ConstantMCP('kcP',-5,5,0,1,M,parameters);

parameters=PduParams_MCP;
subplot(2,3,5)
sweep_paramsX2('kcA',-5,5,'kcP',-5,5,0,1,M,parameters);