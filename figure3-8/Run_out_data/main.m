clear;clc %%%%*****    计算平均持续时间
addpath 'tool\' 'connect\' 'func_analysis\' 'main_fun\'
%% Parameters
Naver = 1000; Np=1; step=0.01; rng(1);
gPoisson=0.5;
NR_node=100;  NS_node=100;

%%
gAMPA_R = 0.0062;   gNMDA_R = 0.0038;   gGABA_R = 0.0037;   gStoR = 0.003; gRtoS = 0.001;
gAMPA_S = gAMPA_R; gNMDA_S = gNMDA_R; gGABA_S = gGABA_R;
% gAMPA_R = 0.0022;   gNMDA_R = 0.0042;   gGABA_R = 0.0018;   gStoR = 0.003; gRtoS = 0.001;
%% 统计量循环
ticd
for kkk = 1:1
    for ppp = 1:Np
        v_rand=rand(Naver);
%         KEtoE=15+ppp*1;
        KEtoE=20;
        KEtoI=20; KItoE=20; KItoI=20;
        parfor jjj = 1:Naver % 并行运算  % parfor 并行只是为了取平均*****************
            rng(jjj+10000);
            [matrixR] = func_WS_network(NR_node,KEtoE,KEtoI,KItoE,KItoI);
            KEtoE1 = 5;
            [matrixS] = func_WS_network(NS_node,KEtoE1,KEtoI,KItoE,KItoI);            
            Per=0.1;
            [matrixRtoS] = func_RtoS(NR_node,NS_node,Per);
            [matrixStoR] = func_StoR(NS_node,NR_node,Per);
            rng(jjj);
            [time(jjj,:),outputS(jjj,:,:),outputR(jjj,:,:),LFP_S(jjj,:),LFP_R(jjj,:),Nt(jjj),break_time(jjj),sum_Energy_ions(jjj),time_fir(jjj),n_fir(jjj)] = ...
            trial_R(gPoisson,NS_node,NR_node,step,gAMPA_S,gNMDA_S,gGABA_S,gRtoS,gAMPA_R,gNMDA_R,gGABA_R,gStoR,matrixS,matrixR,matrixRtoS,matrixStoR);
            T_persistent(ppp,jjj) = break_time(jjj);
            mean_Energy(ppp,jjj) = sum_Energy_ions(jjj);
            matrix_rec_R(kkk,ppp,jjj,:,:)=matrixR;
            matrix_rec_S(kkk,ppp,jjj,:,:)=matrixS;
            matrix_rec_RtoS(kkk,ppp,jjj,:,:)=matrixRtoS;
            matrix_rec_StoR(kkk,ppp,jjj,:,:)=matrixStoR;
            fprintf('进度 %d/%d 持续时间= %f\n',jjj,Naver,T_persistent(ppp,jjj));
        end

        mean_T_persistent=sum(T_persistent)/Naver;
        for jj = 1:Naver
            fprintf('进度 %d/%d 变量1 = %d 持续时间= %f\n',jj,Naver,ppp,T_persistent(ppp,jj));
            fprintf(fp5,'%f %f %f %f %d %d %f\n',T_persistent(ppp,jj),gAMPA_R,gNMDA_R,gGABA_R,kkk,ppp,mean_Energy(ppp,jj));
        end
        fprintf(fp4,'%f %f\n',gNMDA_R,mean_T_persistent);
        toc
    end
end
%%
[T_p,P_Tp] = fun_histogram_Tdur(T_persistent);
figure (8)
subplot (1,1,1), bar(log10(T_p),P_Tp/Naver); axis([log10(10) log10(5100),-inf inf]);  % semilogy

save matrxi_1000_20_0.006.2_3.8_3.7_nonave_ran10000+fr.mat
