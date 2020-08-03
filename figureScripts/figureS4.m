%script to generate Figure S4
%CMJ 20160804
add_paths

clear variables

figure('units','normalized','outerposition',[0 0 1 1])

parameters=PduParams_NoMCP;
subplot(2,3,1)
[rCDE, rPQ, fluxA, AcytoMean, AMCPMean] = GradientsAcrossCellNoMCP(parameters,1);

parameters=PduParams_MCP;
parameters.kcA=10e10;
parameters.kcP=10e10;
subplot(2,3,2)
[rCDE, rPQ, fluxA, AcytoMean, AMCPMean] = GradientsAcrossCell(parameters,1);

parameters=PduParams_MCP;
subplot(2,3,4)
[rCDE, rPQ, fluxA, AcytoMean, AMCPMean] = GradientsAcrossCell(parameters,1);

parameters=PduParams_MCP;
parameters.kcA=10e-10;
parameters.kcP=10e-10;
subplot(2,3,5)
[rCDE, rPQ, fluxA, AcytoMean, AMCPMean] = GradientsAcrossCell(parameters,1);