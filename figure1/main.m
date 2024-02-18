clear;clc
addpath 'tool\' 'connect\' 'func_analysis\' 'main_fun\' 'fun\'
%% Parameters 
Naver=1; step=0.01; rng(1); tic
gPoisson=0.5;
NR_node=100;  NS_node=100;
%% Fixed adjacency matrix parameters
KEtoE=20; KEtoI=20; KItoE=20; KItoI=20;
[matrixR] = func_WS_network(NR_node,KEtoE,KEtoI,KItoE,KItoI);
KEtoE1 = 5;
[matrixS] = func_WS_network(NS_node,KEtoE1,KEtoI,KItoE,KItoI);
Per=0.1;
[matrixRtoS] = func_RtoS(NR_node,NS_node,Per);
[matrixStoR] = func_StoR(NS_node,NR_node,Per);

%%
gAMPA_R = 0.0062;   gNMDA_R = 0.0038;   gGABA_R = 0.0037;   gStoR = 0.003; gRtoS = 0.001;
gAMPA_S = gAMPA_R; gNMDA_S = gNMDA_R; gGABA_S = gGABA_R;
tao_AMPA=2; tao_GABA=10;  tao_NMDA_rise=2;  tao_NMDA_decay=100;
%%
[time,outputS,outputR,Nt,break_time,sum_Energy_ions,time_fir_R,n_fir_R,time_fir_S,n_fir_S] = ...
    trial_R(gPoisson,NS_node,NR_node,step,gAMPA_S,gNMDA_S,gGABA_S,gRtoS,gAMPA_R,gNMDA_R,gGABA_R,gStoR,...
    matrixS,matrixR,matrixRtoS,matrixStoR);
%%

[time_firR,n_firR,output_spike_R] = single_text(time,outputR,Nt,NS_node);
[time_firS,n_firS,output_spike_S] = single_text(time,outputS,Nt,NR_node);

%% Visualization
% time_firR=cell2mat(time_fir_R); n_firR=cell2mat(n_fir_R);
% time_firS=cell2mat(time_fir_S); n_firS=cell2mat(n_fir_S);

figure (20)
subplot(2,1,1),plot(time_firS,n_firS,'ro'); axis([0 6000,-inf inf]);  % %%%%%%%%%%%%%%%%%%%%%%
subplot(2,1,2),plot(time_firR,n_firR,'ro'); axis([0 6000,-inf inf]);  % %%%%%%%%%%%%%%%%%%%%%%
title ('Spiking evolution diagram')

%% firing rate
delta_t=round(10/step); t_min=0; t_max=round(6000/step);

[tS,num_spikeS] = fun_Nspike(output_spike_S,delta_t,t_min,t_max);
[tR,num_spikeR] = fun_Nspike(output_spike_R,delta_t,t_min,t_max);
mean_firingrate_S = num_spikeS/(delta_t*step*0.001*(NS_node));
mean_firingrate_R = num_spikeR/(delta_t*step*0.001*(NR_node));

figure (22)
subplot(4,1,1),plot(tS*step,mean_firingrate_S); axis([-inf inf,-inf inf]);  %
subplot(4,1,2),plot(tR*step,mean_firingrate_R); axis([-inf inf,-inf inf]);  %

%% Output time series data
fp1 = fopen('RasterImage_S.dat','w');
for iii=1:length(time_firS)
    fprintf(fp1,'%f %f\n',time_firS(iii),n_firS(iii));
end
fp2 = fopen('RasterImage_R.dat','w');
for iii=1:length(time_firR)
    fprintf(fp2,'%f %f\n',time_firR(iii),n_firR(iii));
end
fp3 = fopen('FiringRate_S&R.dat','w');
for iii=1:length(tR)
    fprintf(fp3,'%f %f %f\n',tR(iii),mean_firingrate_S(iii),mean_firingrate_R(iii));
end
fclose all;