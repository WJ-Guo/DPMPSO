clc;
clear all;
close all;
psoOptions = get_psoOptionsCompositionFunction8;
SwarmSize = psoOptions.Vars.SwarmSize;
Dim = psoOptions.Vars.Dim;
initrandSWARM = rand(SwarmSize,Dim);
initrandDim = rand(SwarmSize,Dim);
tic
format 	shortE     
MaxCrossValid = 1;

sumfxmin1 = [];sumhistory1=[];sumfes1 = [];
sumfxmin2 = [];sumhistory2=[];sumfes2 = [];
sumfxmin3 = [];sumhistory3=[];sumfes3 = [];
sumfxmin4 = [];sumhistory4=[];sumfes4 = [];
sumfxmin5 = [];sumhistory5=[];sumfes5 = [];
sumfxmin6 = [];sumhistory6=[];sumfes6 = [];
sumfxmin7 = [];sumhistory7=[];sumfes7 = [];
sumfxmin8 = [];sumhistory8=[];sumfes8 = [];
sumfxmin9 = [];sumhistory9=[];sumfes9 = [];

FEs1 = [];
FEs2 = [];
FEs3 = [];
FEs4 = [];
FEs5 = [];
FEs6 = [];
FEs7 = [];
FEs8 = [];
FEs9 = [];

Success1 = 0;
Success2 = 0;
Success3 = 0;
Success4 = 0;
Success5 = 0;
Success6 = 0;
Success7 = 0;
Success8 = 0;
Success9 = 0;

for i=1:MaxCrossValid
%         [fxmin1, xmin1,Swarm1,accept_fes1,accept_success1,history1] = PSO_fes (get_psoOptionsCompositionFunction8,i);
%         [fxmin2, xmin2,Swarm2,accept_fes2,accept_success2,history2] = MSPSO_fes(get_psoOptionsCompositionFunction8,i);
%         [fxmin3, xmin3,Swarm3,accept_fes3,accept_success3,history3] = CLPSO_fes(get_psoOptionsCompositionFunction8,i);
%         [fxmin4, xmin4,Swarm4,accept_fes4,accept_success4,history4] = HCLPSO_fes(get_psoOptionsCompositionFunction8,i);
%         [fxmin5, xmin5,Swarm5,accept_fes5,accept_success5,history5] = ALPSO_fes(get_psoOptionsCompositionFunction8,i);
%         [fxmin6, xmin6,Swarm6,accept_fes6,accept_success6,history6] = FIPS_fes(get_psoOptionsCompositionFunction8,i);
%         [fxmin7, xmin7,Swarm7,accept_fes7,accept_success7,history7] = LSHDE_fes(get_psoOptionsCompositionFunction8,i);
%         [fxmin8, xmin8,Swarm8,accept_fes8,accept_success8,history8] = CLPSO_LOD_fes(get_psoOptionsCompositionFunction8,i);
        [fxmin9, xmin9,Swarm9,accept_fes9,accept_success9,history9] = DPMPSO(get_psoOptionsCompositionFunction8,i);
        
%          sumfxmin1 = [sumfxmin1;fxmin1];sumhistory1=[sumhistory1,history1(:,2)];sumfes1 = [sumfes1,history1(:,1)];
%          sumfxmin2 = [sumfxmin2;fxmin2];sumhistory2=[sumhistory2,history2(:,2)];sumfes2 = [sumfes2,history2(:,1)];
%          sumfxmin3 = [sumfxmin3;fxmin3];sumhistory3=[sumhistory3,history3(:,2)];sumfes3 = [sumfes3,history3(:,1)];
%          sumfxmin4 = [sumfxmin4;fxmin4];sumhistory4=[sumhistory4,history4(:,2)];sumfes4 = [sumfes4,history4(:,1)];
%          sumfxmin5 = [sumfxmin5;fxmin5];sumhistory5=[sumhistory5,history5(:,2)];sumfes5 = [sumfes5,history5(:,1)];
%          sumfxmin6 = [sumfxmin6;fxmin6];sumhistory6=[sumhistory6,history6(:,2)];sumfes6 = [sumfes6,history6(:,1)];
%          sumfxmin7 = [sumfxmin7;fxmin7];sumhistory7=[sumhistory7,history7(:,2)];sumfes7 = [sumfes7,history7(:,1)];
%          sumfxmin8 = [sumfxmin8;fxmin8];sumhistory8=[sumhistory8,history8(:,2)];sumfes8 = [sumfes8,history8(:,1)];
         sumfxmin9 = [sumfxmin9;fxmin9];sumhistory9=[sumhistory9,history9(:,2)];sumfes9 = [sumfes9,history9(:,1)];

%          FEs1(i) = accept_fes1;
%          FEs2(i) = accept_fes2;
%          FEs3(i) = accept_fes3; 
%          FEs4(i) = accept_fes4;
%          FEs5(i) = accept_fes5;
%          FEs6(i) = accept_fes6; 
%          FEs7(i) = accept_fes7;
%          FEs8(i) = accept_fes8;
         FEs9(i) = accept_fes9;
% 
%          Success1 = Success1+accept_success1;
%          Success2 = Success2+accept_success2;
%          Success3 = Success3+accept_success3;
%          Success4 = Success4+accept_success4;
%          Success5 = Success5+accept_success5;
%          Success6 = Success6+accept_success6;
%          Success7 = Success7+accept_success7;
%          Success8 = Success8+accept_success8;        
         Success9 = Success9+accept_success9;       
         close all;
end
% 
%     fprintf('PSO 运行结果:\n');
%     max(sumfxmin1)	
%     min(sumfxmin1)
%     mean(sumfxmin1)	
%     std(sumfxmin1)
%     max(FEs1)
% 	  sumhistory1 = mean(sumhistory1,2);
%     sumfes1 = mean(sumfes1,2);
%     sumhistory1 = [sumfes1 sumhistory1];
%     save sumhistory1PSOSchwefels.dat  sumhistory1 -ascii
% 
%     fprintf('MSPSO 运行结果:\n');
%     max(sumfxmin2)	
%     min(sumfxmin2)
%     mean(sumfxmin2)	
%     std(sumfxmin2)
%     max(FEs2)
% 	sumhistory2 = mean(sumhistory2,2);
%     sumfes2 = mean(sumfes2,2);
%     sumhistory2 = [sumfes2 sumhistory2];
%     save sumhistory2MSPSOSchwefels.dat  sumhistory2 -ascii
%      
%     fprintf('CLPSO 运行结果:\n');
%     max(sumfxmin3)	
%     min(sumfxmin3)
%     mean(sumfxmin3)	
%     std(sumfxmin3)
%     max(FEs3)
%     sumhistory3 = mean(sumhistory3,2);
%     sumfes3 = mean(sumfes3,2);
%     sumhistory3 = [sumfes3 sumhistory3];
%     save sumhistory3CLPSOSchwefels.dat  sumhistory3 -ascii
%   
%     fprintf('HCLPSO 运行结果:\n');
%     max(sumfxmin4)	
%     min(sumfxmin4)
%     mean(sumfxmin4)	
%     std(sumfxmin4)
%     max(FEs4)
%     sumhistory4 = mean(sumhistory4,2);
%     sumfes4 = mean(sumfes4,2);
%     sumhistory4 = [sumfes4 sumhistory4];
%     save sumhistory4HCLPSOSchwefels.dat  sumhistory4 -ascii
%    
%     fprintf('MSNSSA 运行结果:\n');
%     max(sumfxmin5)	
%     min(sumfxmin5)
%     mean(sumfxmin5)	
%     std(sumfxmin5)
%     max(FEs5)
%     sumhistory5 = mean(sumhistory5,2);
%     sumfes5 = mean(sumfes5,2);
%     sumhistory5 = [sumfes5 sumhistory5];
%     save sumhistory5MSNSSASchwefels.dat  sumhistory5 -ascii
% 
%     fprintf('HCOAG 运行结果:\n');
%     max(sumfxmin6)	
%     min(sumfxmin6)
%     mean(sumfxmin6)	
%     std(sumfxmin6)
%     max(FEs6)
%     sumhistory6 = mean(sumhistory6,2);
%     sumfes6 = mean(sumfes6,2);
%     sumhistory6 = [sumfes6 sumhistory6];
%     save sumhistory6HCOAGSchwefels.dat  sumhistory6 -ascii
% 
%     fprintf('L-SHADE 运行结果:\n');
%     max(sumfxmin7)	
%     min(sumfxmin7)
%     mean(sumfxmin7)	
%     std(sumfxmin7)
%     max(FEs7)
%     sumhistory7 = mean(sumhistory7,2);
%     sumfes7 = mean(sumfes7,2);
%     sumhistory7 = [sumfes7 sumhistory7];
%     save sumhistory7LSHADESchwefels.dat  sumhistory7 -ascii
% 
%     fprintf('CMAES 运行结果:\n');
%     max(sumfxmin8)	
%     min(sumfxmin8)
%     mean(sumfxmin8)	
%     std(sumfxmin8)
%     max(FEs8)
%     sumhistory8 = mean(sumhistory8,2);
%     sumfes8 = mean(sumfes8,2);
%     sumhistory8 = [sumfes8 sumhistory8];
%     save sumhistory8CMAESSchwefels.dat  sumhistory8 -ascii

    fprintf('DPMPSO 运行结果:\n');
    max(sumfxmin9)	
    min(sumfxmin9)
    mean(sumfxmin9)	
    std(sumfxmin9)
    max(FEs9)
%     sumhistory9 = mean(sumhistory9,2);
%     sumfes9 = mean(sumfes9,2);
%     sumhistory9 = [sumfes9 sumhistory9];
%     save sumhistory9DPMPSOSchwefels.dat  sumhistory9 -ascii
	
     toc
     t=toc