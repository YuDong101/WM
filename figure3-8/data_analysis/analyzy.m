
% clc;clear; load K=18+7_data_spine.mat
clc;clear; load matrxi_1000_20_0.006.2_3.8_3.7_nonave_ran10000+fr
addpath 'tool\' 'connect\' 'func_analysis\' 'main_fun\' 'BCT' 'BCT\data_and_demos'

kkk=1; sss=1; tic
for ppp = 1:Np
    for sss = 1:Naver
        [matrixx] = matrix_trans(matrix_rec_R,sss,ppp);
        matrix=matrixx;
        matrixx_E=matrix(1:round(0.8*NR_node),1:round(0.8*NR_node));
        matrix_E=matrixx_E;
        [Dds_out(kkk,:),Dds_out_to_E(kkk,:),Dds_out_to_I(kkk,:),Dds_in(kkk,:),Dds_in_from_E(kkk,:),Dds_in_from_I(kkk,:),Dds_avg_out(kkk,:),Dds_avg_in(kkk,:)]...
            = func_Degree_Distribution(matrix); % 出入度的计算 不在乎来源
        [Dds_EfromE(kkk,:),Dds_EfromI(kkk,:),Dds_IfromE(kkk,:),Dds_IfromI(kkk,:)]=func_Degree_Distribution_E_from_I(matrix); % 来源为 E
        [Dds_fromE(kkk,:),Dds_fromI(kkk,:)]=func_Degree_Distribution_EI(matrix);
        [Cc,Cc_avg(kkk,:)]          = func_Cluster_Coeff(matrix);
        [Cc_E,Cc_avg_E(kkk,:)]      = func_Cluster_Coeff(matrix_E);
        [Lens_E,Lens_avg_E(kkk,:)]  = func_Path_Length(matrix_E);
        [Lens,Lens_avg(kkk,:)]      = func_Path_Length(matrix);
        [sum_Num_Cycle(kkk,:),sum_Num_Cycle_E(kkk,:)] = Calculate_Basis_Cycle(matrix);
        [sum_Num_All_Cycle(kkk,:),sum_Num_All_Cycle_E(kkk,:)] = Calculate_All_Cycle(matrix);
        [f_motifs(kkk,:),F_motifs(kkk,:,:)]=motif3struct_bin(matrix);
        [f_motifs_E(kkk,:),F_motifs_E(kkk,:,:)]=motif3struct_bin(matrix_E);
        T_persistent_c(kkk) = T_persistent(ppp,sss);
        mean_Energy_c(kkk) = mean_Energy(ppp,sss);
        kkk=kkk+1;
%         kkk
    end
end
%%
u=1;
v2=7;
v1=4;
A=matrix;
s=uint32(sum(10.^(5:-1:0).*[A(v1,u) A(v2,u) A(u,v1)...
                A(v2,v1) A(u,v2) A(v1,v2)]))
[A(v1,u) A(v2,u) A(u,v1) A(v2,v1) A(u,v2) A(v1,v2)]
toc
%% 计算持续时间的分布

fp7 = fopen('直方图分布 gAMPA=0.0 gGABA=0.0 gPoisson=0.0.dat','w');
[T_p,P_Tp] = fun_histogram_Tdur(T_persistent_c);
for jj = 1:49
     fprintf(fp7,'%f %f\n',T_p(jj),P_Tp(jj)/Naver);
end
fclose(fp7);
figure (8)
subplot (1,1,1), bar(log10(T_p),P_Tp/Naver); axis([log10(10) log10(5100),-inf inf]);  % semilogy
title '持续时间分布'
%% motif
arry_persist=find(T_persistent_c<=4000);
P_motifs=f_motifs_E./f_motifs;
d_motifs=f_motifs-f_motifs_E; d_Fmotifs=F_motifs(:,:,1:80)-F_motifs_E;

flag_X=P_motifs;
figure (7)
subplot (4,4,1), plot (flag_X(arry_persist,1),T_persistent_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-1  % semilogy
subplot (4,4,2), plot (flag_X(arry_persist,2),T_persistent_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-2  % semilogy
subplot (4,4,3), plot (flag_X(arry_persist,3),T_persistent_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-3  % semilogy
subplot (4,4,4), plot (flag_X(arry_persist,4),T_persistent_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-4  % semilogy
subplot (4,4,5), plot (flag_X(arry_persist,5),T_persistent_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-5  % semilogy
subplot (4,4,6), plot (flag_X(arry_persist,6),T_persistent_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-6  % semilogy
subplot (4,4,7), plot (flag_X(arry_persist,7),T_persistent_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-7  % semilogy
subplot (4,4,8), plot (flag_X(arry_persist,8),T_persistent_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-8  % semilogy
subplot (4,4,9), plot (flag_X(arry_persist,9),T_persistent_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-9  % semilogy
subplot (4,4,10), plot (flag_X(arry_persist,10),T_persistent_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-10  % semilogy
subplot (4,4,11), plot (flag_X(arry_persist,11),T_persistent_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-11  % semilogy
subplot (4,4,12), plot (flag_X(arry_persist,12),T_persistent_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-12  % semilogy
subplot (4,4,14), plot (flag_X(arry_persist,13),T_persistent_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-13  % semilogy
subplot (4,4,16), plot (mean_Energy_c,T_persistent_c,'ro'); axis([-inf inf,-inf inf]); title T-E  % semilogy
figure (8)
subplot (4,4,1), plot (flag_X(arry_persist,1),mean_Energy_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-1  % semilogy
subplot (4,4,2), plot (flag_X(arry_persist,2),mean_Energy_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-2  % semilogy
subplot (4,4,3), plot (flag_X(arry_persist,3),mean_Energy_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-3  % semilogy
subplot (4,4,4), plot (flag_X(arry_persist,4),mean_Energy_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-4  % semilogy
subplot (4,4,5), plot (flag_X(arry_persist,5),mean_Energy_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-5  % semilogy
subplot (4,4,6), plot (flag_X(arry_persist,6),mean_Energy_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-6  % semilogy
subplot (4,4,7), plot (flag_X(arry_persist,7),mean_Energy_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-7  % semilogy
subplot (4,4,8), plot (flag_X(arry_persist,8),mean_Energy_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-8  % semilogy
subplot (4,4,9), plot (flag_X(arry_persist,9),mean_Energy_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-9  % semilogy
subplot (4,4,10), plot (flag_X(arry_persist,10),mean_Energy_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-10  % semilogy
subplot (4,4,11), plot (flag_X(arry_persist,11),mean_Energy_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-11  % semilogy
subplot (4,4,12), plot (flag_X(arry_persist,12),mean_Energy_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-12  % semilogy
subplot (4,4,14), plot (flag_X(arry_persist,13),mean_Energy_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-13  % semilogy

aaaa=corr(T_persistent_c',mean_Energy_c',"type","Spearman");
corr_Emotifs=corr(f_motifs_E(arry_persist,:),T_persistent_c(arry_persist)',"type","Spearman");
corr_Emotifs_Energy=corr(f_motifs_E(arry_persist,:),mean_Energy_c(arry_persist)',"type","Spearman");

% [~, p_value] = ttest(T_persistent_c' - corr_Emotifs'.* f_motifs_E);

corr_Pmotifs=corr(P_motifs(arry_persist,:),T_persistent_c(arry_persist)',"type","Spearman");
corr_Pmotifs_Energy=corr(P_motifs(arry_persist,:),mean_Energy_c(arry_persist)',"type","Spearman");

corr_motifs=corr(f_motifs(arry_persist,:),T_persistent_c(arry_persist)',"type","Spearman");
corr_motifs_Energy=corr(f_motifs(arry_persist,:),mean_Energy_c(arry_persist)',"type","Spearman");

corr_dmotifs=corr(d_motifs(arry_persist,:),T_persistent_c(arry_persist)',"type","Spearman");
corr_dmotifs_Energy=corr(d_motifs(arry_persist,:),mean_Energy_c(arry_persist)',"type","Spearman");

fp7 = fopen('motifs.dat','w');
fprintf(fp7,'index corr_dmotifs corr_motifs corr_Emotifs corr_Pmotifs\n'); 
for iii = 1:13
    fprintf(fp7,'%d %f %f %f %f\n',iii,corr_dmotifs(iii),corr_motifs(iii),corr_Emotifs(iii),corr_Pmotifs(iii));        
end
fclose(fp7);

fp8 = fopen('motifs_Energy.dat','w');
fprintf(fp8,'index corr_dmotifs_Energy corr_motifs_Energy corr_Emotifs_Energy corr_Pmotifs_Energy\n'); 
for iii = 1:13
    fprintf(fp8,'%d %f %f %f %f\n',iii,corr_dmotifs_Energy(iii),corr_motifs_Energy(iii),corr_Emotifs_Energy(iii),corr_Pmotifs_Energy(iii));        
end
fclose(fp8);

order_R=1;
figure (9)
subplot(211), bar (1:13,corr_Emotifs.^order_R); axis([-inf inf,-1 1]); title ("Emotifs")
subplot(212), bar (1:13,corr_Emotifs_Energy.^order_R); axis([-inf inf,-1 1]);
figure (10)
subplot(211), bar (1:13,corr_Pmotifs.^order_R); axis([-inf inf,-1 1]); title ("Pmotifs")
subplot(212), bar (1:13,corr_Pmotifs_Energy.^order_R); axis([-inf inf,-1 1]);
figure (11)
subplot(211), bar (1:13,corr_motifs.^order_R); axis([-inf inf,-1 1]); title ("motifs")
subplot(212), bar (1:13,corr_motifs_Energy.^order_R); axis([-inf inf,-1 1]);
figure (12)
subplot(211), bar (1:13,corr_dmotifs.^order_R); axis([-inf inf,-1 1]); title ("dmotifs")
subplot(212), bar (1:13,corr_dmotifs_Energy.^order_R); axis([-inf inf,-1 1]);
%%
ppp=[corr_dmotifs,corr_motifs,corr_Emotifs,corr_Pmotifs];
pppp=[corr_dmotifs_Energy,corr_motifs_Energy,corr_Emotifs_Energy,corr_Pmotifs_Energy];
figure (13)
subplot(211), plot (1:4,ppp); axis([-inf inf,-inf inf]);
subplot(212), plot (1:4,pppp); axis([-inf inf,-inf inf]);

%% E 和 I 的 出度比  目前用处不大
Dds_out_E = sum(Dds_out(:,1:(0.8*NR_node)),2)./sum(Dds_out(:,:),2); Dds_out_I = sum(Dds_out(:,(0.8*NR_node+1):NR_node),2)./sum(Dds_out(:,:),2);
ratioDds_out_EI = Dds_out_E./Dds_out_I;
figure (5)
subplot (2,1,1), plot (ratioDds_out_EI,T_persistent_c,'ro'); axis([-inf inf,-inf inf]);  % semilogy
title 'EI出度比'
% 区分源和终   *****是否要区分大的节点？
% P_in_E=Dds_in_from_E./Dds_in; P_out_E=Dds_out_to_E./Dds_out;
% P_in_E=sum(Dds_in_from_E./Dds_in,2)/NR_node; P_out_E=sum(Dds_out_to_E./Dds_out,2)/NR_node;
P_in_E=sum(Dds_in_from_E(:,1:(0.8*NR_node)),2)./sum(Dds_in,2); P_out_E=sum(Dds_out_to_E(:,1:(0.8*NR_node)),2)./sum(Dds_out,2);
P_out_E_and_in_E = P_in_E.*P_out_E;

figure (7)
subplot (3,1,1), plot (P_in_E,T_persistent_c,'ro'); axis([-inf inf,-inf inf]);  % semilogy
subplot (3,1,2), plot (P_out_E,T_persistent_c,'ro'); axis([-inf inf,-inf inf]);  % semilogy
subplot (3,1,3), plot (P_out_E_and_in_E,T_persistent_c,'ro'); axis([-inf inf,-inf inf]);  % semilogy
title 'EI出度比-区分源和终'
%% E 收到来自于 E 的突触电流比
% P_Dds_EfromE=sum(Dds_EfromE(:,:),2);P_Dds_IfromI=sum(Dds_IfromI(:,:),2);P_Dds_EfromI=sum(Dds_EfromI(:,:),2);P_Dds_IfromE=sum(Dds_IfromE(:,:),2);

% P_Dds_EfromE=sum(Dds_EfromE(:,:),2)./(sum(Dds_EfromE(:,:),2)+sum(Dds_EfromI(:,:),2)+sum(Dds_IfromE(:,:),2)+sum(Dds_IfromI(:,:),2));
aa1=sum(Dds_EfromE(:,:),2); aa2=sum(Dds_EfromI(:,:),2); aa3=sum(Dds_IfromE(:,:),2); aa4=sum(Dds_IfromI(:,:),2);
P_Dds_EfromE=aa1./(aa2+aa3);
P_Dds_IfromI=aa4./(aa2+aa3);
P_Dds_EfromI=aa2./(aa4+aa1);
P_Dds_IfromE=aa3./(aa4+aa1);

P_Dds_IE=(aa4+aa1)./(aa2+aa3); % (aa4+aa1)./(aa2+aa3)

nor_P_Dds_EfromE=(P_Dds_EfromE-min(P_Dds_EfromE))/(max(P_Dds_EfromE)-min(P_Dds_EfromE));

figure (11)
subplot (2,2,1), plot (P_Dds_EfromE',T_persistent_c,'ro'); axis([-inf inf,-inf inf]); title 'E to E概率'  % semilogy
subplot (2,2,2), plot (P_Dds_IfromI,T_persistent_c,'ro'); axis([-inf inf,-inf inf]); title 'I to I概率' % semilogy
subplot (2,2,3), plot (P_Dds_EfromI,T_persistent_c,'ro'); axis([-inf inf,-inf inf]); title 'I to E概率'  % semilogy
subplot (2,2,4), plot (P_Dds_IfromE,T_persistent_c,'ro'); axis([-inf inf,-inf inf]); title 'E to I概率' % semilogy
figure (12)
subplot (2,2,1), plot (aa1,T_persistent_c,'ro'); axis([-inf inf,-inf inf]); title 'E to E概率'  % semilogy
subplot (2,2,2), plot (aa4,T_persistent_c,'ro'); axis([-inf inf,-inf inf]); title 'I to I概率' % semilogy
subplot (2,2,3), plot (aa2,T_persistent_c,'ro'); axis([-inf inf,-inf inf]); title 'I to E概率'  % semilogy
subplot (2,2,4), plot (aa3,T_persistent_c,'ro'); axis([-inf inf,-inf inf]); title 'E to I概率' % semilogy
[corr(P_Dds_EfromE,T_persistent_c',"type","Spearman")
corr(P_Dds_IfromI,T_persistent_c',"type","Spearman")
corr(P_Dds_EfromI,T_persistent_c',"type","Spearman")
corr(P_Dds_IfromE,T_persistent_c',"type","Spearman")
corr(P_Dds_IE,T_persistent_c',"type","Spearman")]

[corr(aa1,T_persistent_c',"type","Spearman")
corr(aa4,T_persistent_c',"type","Spearman")
corr(aa2,T_persistent_c',"type","Spearman")
corr(aa3,T_persistent_c',"type","Spearman")]

% x_min=min(P_Dds_EfromE);x_max=max(P_Dds_EfromE);dx=0.003;
% [out_x,out_y,output_histogram] = fun_histogram_Bifurcation(P_Dds_EfromE,T_persistent_c,x_min,x_max,dx);
% figure (110)
% subplot(111),contour(out_x,out_y,output_histogram',10);

fp7 = fopen('EtoE.dat','w');
sss=1;
for ppp = 1:Np
    for jjj = 1:Naver
        fprintf(fp7,'%f %f\n',P_Dds_EfromE(sss),T_persistent_c(sss));
        sss=sss+1;
    end
end
fclose(fp7);

fp7 = fopen('ItoI.dat','w');
sss=1;
for ppp = 1:Np
    for jjj = 1:Naver
        fprintf(fp7,'%f %f\n',P_Dds_IfromI(sss),T_persistent_c(sss));
        sss=sss+1;
    end
end
fclose(fp7);

% fp7 = fopen('EtoE分岔.dat','w');
% sss=1;
% for ppp = 1:length(out_x)
%     for jjj = 1:length(out_y)
%         fprintf(fp7,'%f %f %f\n',out_x(ppp),out_y(jjj),output_histogram(ppp,jjj));
%         sss=sss+1;
%     end
% end
% fclose(fp7);

figure (13)
subplot (1,1,1), semilogx (P_Dds_IE,T_persistent_c,'ro'); axis([-inf inf,-inf inf]); title 'E to E概率'  % semilogy
%% 每个节点的突触电流的比例
NE=round(0.8*NR_node);
P_Dds_EI_single=Dds_fromE./(Dds_fromE+Dds_fromI);

[ai,aj] = size(P_Dds_EI_single); aa(1:ai,1:aj) = 0; aa(P_Dds_EI_single(:,1:NE)>0.8)=1;  bb(1:ai,1:aj) = 0; bb(P_Dds_EI_single(:,NE+1:100)>0.8)=1;
P_Dds_E_large=sum(aa,2)./NE; P_Dds_I_large=sum(bb,2)./20;
Raito_richHub=P_Dds_E_large./P_Dds_I_large;

figure (6)
subplot (311), plot (P_Dds_E_large,T_persistent_c,'ro'); axis([-inf inf,-inf inf]);  % semilogy
subplot (312), plot (P_Dds_I_large,T_persistent_c,'ro'); axis([-inf inf,-inf inf]);  % semilogy
subplot (313), plot (Raito_richHub,T_persistent_c,'ro'); axis([-inf inf,-inf inf]);  % semilogy
title 'E比例大的神经元的占比'

x_min=min(P_Dds_E_large);x_max=max(P_Dds_E_large);dx=0.015;
[out_x,out_y,output_histogram] = fun_histogram_Bifurcation(P_Dds_E_large,T_persistent_c,x_min,x_max,dx);
% figure (110)
% subplot(111),contour(out_x,out_y,output_histogram',10);
  
[corr(P_Dds_E_large,T_persistent_c',"type","Spearman")
corr(P_Dds_I_large,T_persistent_c',"type","Spearman")
corr(Raito_richHub,T_persistent_c',"type","Spearman")]

fp7 = fopen('LargeE.dat','w');
sss=1;
for ppp = 1:Np
    for jjj = 1:Naver
        fprintf(fp7,'%f %f\n',P_Dds_E_large(sss),T_persistent_c(sss));
        sss=sss+1;
    end
end
fclose(fp7);

% ratio_I=zeros(Naver*ppp); ratio_E=zeros(Naver*ppp); ratio_E_I=zeros(Naver*ppp); H_topology=zeros(Naver*ppp);

%% 聚类系数和长度
nor_Cc_avg=(Cc_avg-min(Cc_avg))/(max(Cc_avg)-min(Cc_avg));
nor_Cc_avg_E=(Cc_avg_E-min(Cc_avg_E))/(max(Cc_avg_E)-min(Cc_avg_E));
nor_Lens_avg=(Lens_avg-min(Lens_avg))/(max(Lens_avg)-min(Lens_avg));
nor_Lens_avg_E=(Lens_avg_E-min(Lens_avg_E))/(max(Lens_avg_E)-min(Lens_avg_E));
P_lens=Lens_avg_E./Lens_avg; P_Cc=Cc_avg_E./Cc_avg;

Small_World=P_lens./P_Cc;

figure (13)
subplot (2,3,1), plot (Cc_avg,T_persistent_c,'ro'); axis([-inf inf,-inf inf]);  % semilogy
subplot (2,3,2), plot (Cc_avg_E,T_persistent_c,'ro'); axis([-inf  inf,-inf inf]);  % semilogy
subplot (2,3,3), plot (P_Cc,T_persistent_c,'ro'); axis([-inf inf,-inf inf]);  % semilogy
subplot (2,3,4), plot (Lens_avg,T_persistent_c,'ro'); axis([-inf inf,-inf inf]);  % semilogy
subplot (2,3,5), plot (Lens_avg_E,T_persistent_c,'ro'); axis([-inf inf,-inf inf]);  % semilogy
subplot (2,3,6), plot (P_lens,T_persistent_c,'ro'); axis([-inf inf,-inf inf]);  % semilogy

figure (14)
subplot (1,1,1), plot (Small_World,T_persistent_c,'ro'); axis([-inf inf,-inf inf]);  % semilogy
% corr(Cc_avg,T_persistent_c',"type","Spearman");

fp6 = fopen('小世界性质.dat','w');
fprintf(fp6,'index Cc_avg Cc_avg_E P_Cc Lens_avg Lens_avg_E P_lens T_persistent_c\n');
sss=1;
for ppp = 1:Np
    for jjj = 1:Naver
        fprintf(fp6,'%d %f %f %f %f %f %f %f\n',sss,Cc_avg(sss),Cc_avg_E(sss),P_Cc(sss),Lens_avg(sss),Lens_avg_E(sss),P_lens(sss),T_persistent_c(sss));
        sss=sss+1;
    end
end
fclose (fp6);
%%
Cc_corr_All=corr(Cc_avg,T_persistent_c',"type","Spearman"); Cc_corr_E=corr(Cc_avg_E,T_persistent_c',"type","Spearman"); Cc_corr_P=corr(P_Cc,T_persistent_c',"type","Spearman");
lens_corr_All=corr(Lens_avg,T_persistent_c',"type","Spearman"); lens_corr_E=corr(Lens_avg_E,T_persistent_c',"type","Spearman"); lens_corr_P=corr(P_lens,T_persistent_c',"type","Spearman");
corr_Cc=[Cc_corr_All,Cc_corr_E,Cc_corr_P];
corr_lens=[lens_corr_All,lens_corr_E,lens_corr_P];
figure (14)
subplot(211), plot (1:3,corr_Cc); axis([-inf inf,-inf inf]);
subplot(212), plot (1:3,corr_lens); axis([-inf inf,-inf inf]);

%% 所有循环
% P_Num_Cycle(sum_Num_Cycle~=0)=sum_Num_Cycle_E(sum_Num_Cycle~=0)./sum_Num_Cycle(sum_Num_Cycle~=0);
P_Num_All_Cycle=sum_Num_All_Cycle_E(:,3)./(sum_Num_All_Cycle(:,3)+sum_Num_All_Cycle_E(:,3));
nor_P_Num_All_Cycle=(P_Num_All_Cycle-min(P_Num_All_Cycle))./(max(P_Num_All_Cycle)-min(P_Num_All_Cycle));
ratio_Cycle=sum_Num_All_Cycle_E(:,3)./(sum_Num_All_Cycle(:,3));

figure (14)
subplot (2,1,1), plot (P_Num_All_Cycle,T_persistent_c,'ro'); axis([-inf inf,-inf inf]);  % semilogy
subplot (2,1,2), plot (ratio_Cycle,T_persistent_c,'ro'); axis([-inf inf,-inf inf]);  % semilogy
figure (15)
subplot (3,1,1), plot (sum_Num_All_Cycle_E(:,3),T_persistent_c,'ro'); axis([-inf inf,-inf inf]);  % semilogy
subplot (3,1,2), plot (sum_Num_All_Cycle(:,3),T_persistent_c,'ro'); axis([-inf inf,-inf inf]);  % semilogy
subplot (3,1,3), plot (sum_Num_All_Cycle(:,3)+sum_Num_All_Cycle_E(:,3),T_persistent_c,'ro'); axis([-inf inf,-inf inf]);  % semilogy

Num_All_Cycle_E=sum_Num_All_Cycle_E(:,3);
Num_All_Cycle_d=sum_Num_All_Cycle(:,3);
Num_All_Cycle=sum_Num_All_Cycle(:,3)+sum_Num_All_Cycle_E(:,3);

fp6 = fopen('神经环路.dat','w');
fprintf(fp6,'index cycleE cycleEI cycleALL PcycleE T_persistent_c\n');
sss=1;
for ppp = 1:Np
    for jjj = 1:Naver
        fprintf(fp6,'%d %f %f %f %f %f\n',sss,sum_Num_All_Cycle_E(sss,3),sum_Num_All_Cycle(sss,3),sum_Num_All_Cycle_E(sss,3)+sum_Num_All_Cycle(sss,3),P_Num_All_Cycle(sss),T_persistent_c(sss));
        sss=sss+1;
    end
end
fclose (fp6);

[corr(sum_Num_All_Cycle(:,3),T_persistent_c',"type","Spearman")
corr(sum_Num_All_Cycle_E(:,3),T_persistent_c',"type","Spearman")
corr(sum_Num_All_Cycle(:,3)+sum_Num_All_Cycle_E(:,3),T_persistent_c',"type","Spearman")
corr(P_Num_All_Cycle,T_persistent_c',"type","Spearman")]
%% 数据输出
ap=find(T_persistent_c==4000);  at=find(T_persistent_c<4000);
fp7 = fopen('Energy-transisent Scatter data.dat','w');
fprintf(fp7,'index Cc_avg Cc_avg_E P_Cc Lens_avg Lens_avg_E P_lens P_Dds_IE P_Dds_E_large P_Dds_I_large Raito_richHub Num_All_Cycle Num_All_Cycle_d Num_All_Cycle_E P_Num_All_Cycle aa1 aa2 aa3 aa4 T_persistent_c Energy\n'); 
for jj = 1:length(at)
     fprintf(fp7,'%d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',at(jj),Cc_avg(at(jj)),Cc_avg_E(at(jj)),P_Cc(at(jj)),Lens_avg(at(jj)),Lens_avg_E(at(jj)),P_lens(at(jj)),P_Dds_IE(at(jj)),P_Dds_E_large(at(jj)),P_Dds_I_large(at(jj)),Raito_richHub(at(jj)),Num_All_Cycle(at(jj)),Num_All_Cycle_d(at(jj)),Num_All_Cycle_E(at(jj)),P_Num_All_Cycle(at(jj)),aa1(at(jj)),aa2(at(jj)),aa3(at(jj)),aa4(at(jj)),T_persistent_c(at(jj)),mean_Energy_c(at(jj)));
end
fclose(fp7);

fp7 = fopen('Energy-persistent Scatter data.dat','w');
fprintf(fp7,'index Cc_avg Cc_avg_E P_Cc Lens_avg Lens_avg_E P_lens P_Dds_IE P_Dds_E_large P_Dds_I_large Raito_richHub Num_All_Cycle Num_All_Cycle_d Num_All_Cycle_E P_Num_All_Cycle aa1 aa2 aa3 aa4 T_persistent_c Energy\n'); 
for jj = 1:length(ap)
     fprintf(fp7,'%d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',ap(jj),Cc_avg(ap(jj)),Cc_avg_E(ap(jj)),P_Cc(ap(jj)),Lens_avg(ap(jj)),Lens_avg_E(ap(jj)),P_lens(ap(jj)),P_Dds_IE(ap(jj)),P_Dds_E_large(ap(jj)),P_Dds_I_large(ap(jj)),Raito_richHub(ap(jj)),Num_All_Cycle(ap(jj)),Num_All_Cycle_d(ap(jj)),Num_All_Cycle_E(ap(jj)),P_Num_All_Cycle(ap(jj)),aa1(ap(jj)),aa2(ap(jj)),aa3(ap(jj)),aa4(ap(jj)),T_persistent_c(ap(jj)),mean_Energy_c(ap(jj)));
end
fclose(fp7);

fp7 = fopen('motifs-全网络的表达数量.dat','w');
fprintf(fp7,'index 1 2 3 4 5 6 7 8 9 10 11 12 13 T_persistent_c Energy\n'); 
for jj = 1:length(f_motifs(:,1))
     fprintf(fp7,'%d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',jj,f_motifs(jj,1),f_motifs(jj,2),f_motifs(jj,3),f_motifs(jj,4),f_motifs(jj,5),f_motifs(jj,6),f_motifs(jj,7),f_motifs(jj,8),f_motifs(jj,9),f_motifs(jj,10),f_motifs(jj,11),f_motifs(jj,12),f_motifs(jj,13),T_persistent_c(jj),mean_Energy_c(jj));
end
fclose(fp7);

fp7 = fopen('motifs-E网络的表达数量.dat','w');
fprintf(fp7,'index 1 2 3 4 5 6 7 8 9 10 11 12 13 T_persistent_c Energy\n'); 
for jj = 1:length(f_motifs_E(:,1))
     fprintf(fp7,'%d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',jj,f_motifs_E(jj,1),f_motifs_E(jj,2),f_motifs_E(jj,3),f_motifs_E(jj,4),f_motifs_E(jj,5),f_motifs_E(jj,6),f_motifs_E(jj,7),f_motifs_E(jj,8),f_motifs_E(jj,9),f_motifs_E(jj,10),f_motifs_E(jj,11),f_motifs_E(jj,12),f_motifs_E(jj,13),T_persistent_c(jj),mean_Energy_c(jj));
end
fclose(fp7);

fp7 = fopen('motifs-I涉及的表达数量.dat','w');
fprintf(fp7,'index 1 2 3 4 5 6 7 8 9 10 11 12 13 T_persistent_c Energy\n');
for jj = 1:length(d_motifs(:,1))
     fprintf(fp7,'%d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',jj,d_motifs(jj,1),d_motifs(jj,2),d_motifs(jj,3),d_motifs(jj,4),d_motifs(jj,5),d_motifs(jj,6),d_motifs(jj,7),d_motifs(jj,8),d_motifs(jj,9),d_motifs(jj,10),d_motifs(jj,11),d_motifs(jj,12),d_motifs(jj,13),T_persistent_c(jj),mean_Energy_c(jj));
end
fclose(fp7);

fp7 = fopen('三种motifs表达的平均值和标准差.dat','w');
mean_motifs=mean(f_motifs,1); mean_Emotifs=mean(f_motifs_E,1); mean_dmotifs=mean(d_motifs,1);
std_motifs=std(f_motifs,0,1); std_Emotifs=std(f_motifs_E,0,1); std_dmotifs=std(d_motifs,0,1);
for jj = 1:length(d_motifs(1,:))
     fprintf(fp7,'%d %f %f %f %f %f %f\n',jj,mean_motifs(jj),std_motifs(jj),mean_Emotifs(jj),std_Emotifs(jj),mean_dmotifs(jj),std_dmotifs(jj));
end
fclose(fp7);

%% 找相关联最大值
% corr_max=-1; tic
% % arry_persist=1:length(T_persistent_c);
% clear tex_topology tttt Ecorr_topology F_EfromE F_All_Cycle F_Cc_avg_E F_EI_large F_Lens_avg_E arry_persist
% arry_persist=find(T_persistent_c<=4000); muti=1.0; plus=0.0;
% F_EfromE=nor_P_Dds_EfromE(arry_persist)*muti+plus; F_All_Cycle(:,1)=nor_P_Num_All_Cycle(arry_persist)*muti+plus; % P_motifs(arry_persist,7)
% F_Cc_avg_E=nor_P_Cc(arry_persist)*muti+plus; F_EI_large=nor_P_Dds_EI_large(arry_persist)*muti+plus;
% F_Lens_avg_E=nor_P_lens(arry_persist)*muti+plus;
% 
% ii1=1;
% for a1=-0.1:0.1:1
%     ii2=1;
%     for a2=-0.1:0.1:1
%         ii3=1;
%         for a3=-0.5:0.1:0.5
%             ii4=1;
%             for a4=-1:0.1:0.1
%                 ii5=1;
%                 for a5=-1:0.1:0.1
%                     tex_topology(ii1,ii2,ii3,ii4,ii5,:)=F_EfromE.^a1.*...
%                         F_All_Cycle.^a2.*...
%                             F_Cc_avg_E.^a3.*...
%                             F_EI_large.^a4.*...
%                             F_Lens_avg_E.^a5;
%                     tttt(:)=tex_topology(ii1,ii2,ii3,ii4,ii5,:);
%                     corr_topology(ii1,ii2,ii3,ii4,ii5)=corr(tttt',T_persistent_c(arry_persist)',"type","Pearson");
%                     if corr_topology(ii1,ii2,ii3,ii4,ii5)>corr_max
%                         corr_max=corr_topology(ii1,ii2,ii3,ii4,ii5);
%                         corr_ma1=a1; corr_ma2=a2 ;corr_ma3=a3; corr_ma4=a4; corr_ma5=a5;
%                     end
%                     ii5=ii5+1;
%                 end
%                 ii4=ii4+1;
%             end
%             ii3=ii3+1;
%         end
%         ii2=ii2+1;
%     end
%     ii1=ii1+1;
% end
% %%
% acaca=...
% (nor_P_Dds_EfromE*muti+plus).^corr_ma1.*...
% (nor_P_Num_All_Cycle(:)*muti+plus).^corr_ma2.*...
% (nor_P_Cc*muti+plus).^corr_ma3.*...
% (nor_P_Dds_EI_large*muti+plus).^corr_ma4.*...
% (nor_P_lens*muti+plus).^corr_ma5;
% acaca(acaca==inf)=0;
% toc
% %% 相关系数r和r方
% % R2_acaca=corr(acaca(arry_persist),T_persistent_c(arry_persist)',"type","Pearson")^2;
% % R2_All_Cycle=corr(P_Num_All_Cycle(arry_persist),T_persistent_c(arry_persist)',"type","Pearson")^2;
% % R2_EfromE=corr(P_Dds_EfromE(arry_persist),T_persistent_c(arry_persist)',"type","Pearson")^2;
% % R2_EI_large=corr(P_Dds_EI_large(arry_persist),T_persistent_c(arry_persist)',"type","Pearson")^2;
% % R2_Cc=corr(Cc_avg_E(arry_persist),T_persistent_c(arry_persist)',"type","Pearson")^2;
% % R2_Lens=corr(Lens_avg_E(arry_persist),T_persistent_c(arry_persist)',"type","Pearson")^2;
% % R_acaca=corr(acaca(arry_persist),T_persistent_c(arry_persist)',"type","Pearson");
% % R_All_Cycle=corr(P_Num_All_Cycle(arry_persist),T_persistent_c(arry_persist)',"type","Pearson");
% % R_EfromE=corr(P_Dds_EfromE(arry_persist),T_persistent_c(arry_persist)',"type","Pearson");
% % R_EI_large=corr(P_Dds_EI_large(arry_persist),T_persistent_c(arry_persist)',"type","Pearson");
% % R_Cc=corr(Cc_avg_E(arry_persist),T_persistent_c(arry_persist)',"type","Pearson");
% % R_Lens=corr(Lens_avg_E(arry_persist),T_persistent_c(arry_persist)',"type","Pearson");
% 
% R_acaca=corr(acaca(arry_persist),T_persistent_c(arry_persist)',"type","Spearman");
% R_All_Cycle=corr(P_Num_All_Cycle(arry_persist),T_persistent_c(arry_persist)',"type","Spearman");
% R_EfromE=corr(P_Dds_EfromE(arry_persist),T_persistent_c(arry_persist)',"type","Spearman");
% R_EI_large=corr(P_Dds_EI_large(arry_persist),T_persistent_c(arry_persist)',"type","Spearman");
% R_Cc=corr(P_Cc(arry_persist),T_persistent_c(arry_persist)',"type","Spearman");
% R_Lens=corr(P_lens(arry_persist),T_persistent_c(arry_persist)',"type","Spearman");
% 
% bar_dat=[R_EfromE R_All_Cycle R_Cc R_EI_large R_Lens 0 R_acaca];
% % bar2_dat=[R2_EfromE R2_All_Cycle R2_Cc R2_EI_large R2_Lens 0 R2_acaca];
% bar_weight=[corr_ma1 corr_ma2 corr_ma3 corr_ma4 corr_ma5 0 0];
% figure (18)
% subplot (2,1,1), bar (1:length(bar_dat),bar_dat);  axis([0 8,-1 1]);  % semilog
% subplot (2,1,2),  bar (1:length(bar_dat),bar_weight);  axis([0 8,-1 1]);  % semilog
% 
% figure (15)
% subplot (2,2,1), plot (acaca(arry_persist),T_persistent_c(arry_persist),'ro'); axis([-inf inf,-inf inf]);  % semilogy
% 
% % x_min=min(acaca);x_max=max(acaca);dx=0.1;
% x_min=0;x_max=1.9;dx=0.1;
% [out_x,out_y,output_histogram] = fun_histogram_Bifurcation(acaca,T_persistent_c,x_min,x_max,dx);
% figure (110)
% % subplot(111),contour(out_x,out_y,output_histogram',10);
% subplot(111),heatmap(out_x,out_y,output_histogram');
% 
% fp7 = fopen('F_topology-Duration.dat','w');
% sss=1;
% for ppp = 1:Np
%     for jjj = 1:Naver
%         fprintf(fp7,'%f %f\n',acaca(sss),T_persistent_c(sss));
%         sss=sss+1;
%     end
% end
% fclose(fp7);
% 
% fp7 = fopen('W_topology.dat','w');
% fprintf(fp7,'Factor Correlation Weight\n');
%     fprintf(fp7,'Cc %f %f\n',R_Cc,corr_ma3);
%     fprintf(fp7,'Len %f %f\n',R_Lens,corr_ma5);
%     fprintf(fp7,'EtoE %f %f\n',R_EfromE,corr_ma1);
%     fprintf(fp7,'LargeE %f %f\n',R_EI_large,corr_ma4);
%     fprintf(fp7,'Cycle %f %f\n',R_All_Cycle,corr_ma2);
%     fprintf(fp7,'F_topology %f %f\n',R_acaca,0);
% fclose(fp7);
% 
% fp7 = fopen('Ftopology分岔.dat','w');
% sss=1;
% for ppp = 1:length(out_x)
%     for jjj = 1:length(out_y)
%         fprintf(fp7,'%f %f %f\n',out_x(ppp),out_y(jjj),output_histogram(ppp,jjj));
%         sss=sss+1;
%     end
% end
% fclose(fp7);
% 
% 
% %% 能量相关联最大值
% Emax_corr=-1; tic
% % arry_persist=1:length(T_persistent_c);
% clear Etex_topology Etttt Ecorr_topology F_EfromE F_All_Cycle F_Cc_avg_E F_EI_large F_Lens_avg_E arry_persist
% arry_persist=find(T_persistent_c<=4000); muti=0.9; plus=0.1;
% F_EfromE=nor_P_Dds_EfromE(arry_persist)*muti+plus; F_All_Cycle(:,1)=nor_P_Num_All_Cycle(arry_persist)*muti+plus; % P_motifs(arry_persist,7)
% F_Cc_avg_E=nor_P_Cc(arry_persist)*muti+plus; F_EI_large=nor_P_Dds_EI_large(arry_persist)*muti+plus;
% F_Lens_avg_E=nor_P_lens(arry_persist)*muti+plus;
% % -1:0.1:1
% ii1=1;
% for a1=-0.1:0.1:1
%     ii2=1;
%     for a2=-0.1:0.1:1
%         ii3=1;
%         for a3=-0.5:0.1:0.5
%             ii4=1;
%             for a4=-1:0.1:0.1
%                 ii5=1;
%                 for a5=-1:0.1:0.1
%                     Etex_topology(ii1,ii2,ii3,ii4,ii5,:)=...
%                         F_EfromE.^a1.*...
%                         F_All_Cycle.^a2.*...
%                         F_Cc_avg_E.^a3.*...
%                         F_EI_large.^a4.*...
%                         F_Lens_avg_E.^a5;
%                     Etttt(:)=Etex_topology(ii1,ii2,ii3,ii4,ii5,:);
%                     Etttt(Etttt==inf)=0;
%                     Ecorr_topology(ii1,ii2,ii3,ii4,ii5)=corr(Etttt',mean_Energy_c(arry_persist)',"type","Spearman");
%                     if Ecorr_topology(ii1,ii2,ii3,ii4,ii5)>=Emax_corr
%                         Emax_corr=Ecorr_topology(ii1,ii2,ii3,ii4,ii5);
%                         Ecorr_ma1=a1; Ecorr_ma2=a2 ;Ecorr_ma3=a3; Ecorr_ma4=a4; Ecorr_ma5=a5;
%                     end
%                     ii5=ii5+1;
%                 end
%                 ii4=ii4+1;
%             end
%             ii3=ii3+1;
%         end
%         ii2=ii2+1;
%     end
%     ii1=ii1+1;
% end
% %
% Eacaca=...
% (nor_P_Dds_EfromE*muti+plus).^Ecorr_ma1.*...
% (nor_P_Num_All_Cycle(:)*muti+plus).^Ecorr_ma2.*...
% (nor_P_Cc*muti+plus).^Ecorr_ma3.*...
% (nor_P_Dds_EI_large*muti+plus).^Ecorr_ma4.*...
% (nor_P_lens*muti+plus).^Ecorr_ma5;
% Eacaca(Eacaca==inf)=0;
% toc
% %%  能量
% 
% % arry_persist=1:length(T_persistent_c);
% arry_persist=find(T_persistent_c<=4000);
% E_sta=-inf;
% 
% figure (16)
% subplot (2,2,1), plot (mean_Energy_c,T_persistent_c,'ro'); axis([E_sta inf,-inf inf]);  % semilogy
% subplot (2,2,2), plot (P_Num_All_Cycle(arry_persist),mean_Energy_c(arry_persist),'ro'); axis([-inf inf,E_sta inf]);  % semilogy
% subplot (2,2,3), plot (P_Dds_EfromE(arry_persist),mean_Energy_c(arry_persist),'ro'); axis([-inf inf,E_sta inf]);  % semilogy
% subplot (2,2,4), plot (P_Dds_EI_large(arry_persist),mean_Energy_c(arry_persist),'ro'); axis([-inf inf,E_sta inf]);  % semilogy
% 
% figure (17)
% subplot (2,2,1), plot (Cc_avg_E(arry_persist),mean_Energy_c(arry_persist),'ro'); axis([-inf inf,E_sta inf]);  % semilogy
% subplot (2,2,2), plot (Lens_avg_E(arry_persist),mean_Energy_c(arry_persist),'ro'); axis([-inf inf,E_sta inf]);  % semilogy
% 
% figure (18)
% subplot (2,2,1), plot (Eacaca(arry_persist),mean_Energy_c(arry_persist)','ro'); axis([-inf inf,E_sta inf]);  % semilog
% 
% 
% % R_Energy_Eacaca=corr(Eacaca(arry_persist),mean_Energy_c(arry_persist)',"type","Pearson");
% % R_Energy_All_Cycle=corr(P_Num_All_Cycle(arry_persist),mean_Energy_c(arry_persist)',"type","Pearson");
% % R_Energy_EfromE=corr(P_Dds_EfromE(arry_persist),mean_Energy_c(arry_persist)',"type","Pearson");
% % R_Energy_EI_large=corr(P_Dds_EI_large(arry_persist),mean_Energy_c(arry_persist)',"type","Pearson");
% % R_Energy_Cc=corr(Cc_avg_E(arry_persist),mean_Energy_c(arry_persist)',"type","Pearson");
% % R_Energy_Lens=corr(Lens_avg_E(arry_persist),mean_Energy_c(arry_persist)',"type","Pearson");
% % 
% R_Energy_Eacaca=corr(Eacaca(arry_persist),mean_Energy_c(arry_persist)',"type","Spearman");
% R_Energy_All_Cycle=corr(P_Num_All_Cycle(arry_persist),mean_Energy_c(arry_persist)',"type","Spearman");
% R_Energy_EfromE=corr(P_Dds_EfromE(arry_persist),mean_Energy_c(arry_persist)',"type","Spearman");
% R_Energy_EI_large=corr(P_Dds_EI_large(arry_persist),mean_Energy_c(arry_persist)',"type","Spearman");
% R_Energy_Cc=corr(P_Cc(arry_persist),mean_Energy_c(arry_persist)',"type","Spearman");
% R_Energy_Lens=corr(P_lens(arry_persist),mean_Energy_c(arry_persist)',"type","Spearman");
% 
% fp7 = fopen('F_topology-Energy.dat','w');
% fprintf(fp7,'Cc Lens EtoE LargeE Cycle Eaca Duration\n');
% sss=1;
% for ppp = 1:Np
%     for jjj = 1:Naver
%         fprintf(fp7,'%f %f %f %f %f %f %f\n',P_Cc(sss),P_lens(sss),P_Dds_EfromE(sss),P_Dds_EI_large(sss),P_Num_All_Cycle(sss),Eacaca(sss),mean_Energy_c(sss));
%         sss=sss+1;
%     end
% end
% fclose(fp7);
% 
% fp7 = fopen('W_topology-Energy.dat','w');
% fprintf(fp7,'Factor Correlation Weight\n');
%     fprintf(fp7,'Cc %f %f\n',R_Energy_Cc,Ecorr_ma3);
%     fprintf(fp7,'Len %f %f\n',R_Energy_Lens,Ecorr_ma5);
%     fprintf(fp7,'EtoE %f %f\n',R_Energy_EfromE,Ecorr_ma1);
%     fprintf(fp7,'LargeE %f %f\n',R_Energy_EI_large,Ecorr_ma4);
%     fprintf(fp7,'Cycle %f %f\n',R_Energy_All_Cycle,Ecorr_ma2);
%     fprintf(fp7,'F_topology %f %f\n',R_Energy_Eacaca,0);
% fclose(fp7);
% 
% bar_dat_Energy=[R_Energy_EfromE R_Energy_All_Cycle R_Energy_Cc R_Energy_EI_large R_Energy_Lens 0 R_Energy_Eacaca];
% bar2_dat_Energy=[R2_Energy_EfromE R2_Energy_All_Cycle R2_Energy_Cc R2_Energy_EI_large R2_Energy_Lens 0 R2_Energy_Eacaca];
% bar_weight=[Ecorr_ma1 Ecorr_ma2 Ecorr_ma3 Ecorr_ma4 Ecorr_ma5 0 0];
% figure (19)
% subplot (2,1,1), bar (1:length(bar2_dat_Energy),bar_dat_Energy);  axis([0 8,-1 1]);  % semilog
% subplot (2,1,2),  bar (1:length(bar2_dat_Energy),bar_weight);  axis([0 8,-1 1]);  % semilog
% %% 分岔
% Bif_flag=P_Dds_EfromE;
% x_min=min(Bif_flag);x_max=max(Bif_flag);dx=(x_max-x_min)/40;
% [out_x,out_y,output_histogram] = fun_histogram_Bifurcation(Bif_flag,T_persistent_c,x_min,x_max,dx);
% figure (111)
% subplot(111),contour(out_x,out_y,output_histogram',20); axis([-inf inf,-inf inf]);
% fp6 = fopen('三维分岔.dat','w');
% for ppp = 1:length(out_x)
%     for jjj= 1:length(out_y)
%         fprintf(fp6,'%f %f %f\n',out_x(ppp),out_y(jjj),output_histogram(ppp,jjj));
%     end
% end
% clear out_x out_y output_histogram
% fclose (fp6);
% %% 输出数据
% fp6 = fopen('节点参数.dat','w');
% fprintf(fp6,'index index H_topology P_Num_Cycle3 P_Num_All_Cycle ratioDds_E_I P_Dds_EfromE P_Dds_EI_large acaca T_persistent_c\n');
% sss=1;
% for ppp = 1:Np
%     for jjj = 1:Naver
%         fprintf(fp6,'%d %d %f %f %f %f %f %f %f %f\n',ppp,jjj,H_topology(sss),P_Num_Cycle(sss),P_Num_All_Cycle(sss),ratioDds_out_EI(sss),P_Dds_EfromE(sss),P_Dds_EI_large(sss),acaca(sss,3),T_persistent_c(ppp,jjj));
% %         fprintf(fp6,'%d %d %f %f %f %f %f %f\n',ppp,jjj,H_topology(sss),ratio_E_I(sss),ratio_E(sss),ratio_I(sss),ratioDds_E_I(sss),P_Dds_EfromE(sss));
%         ssT_persistent_c(sss) = T_persistent_c(ppp,jjj);
%         sss=sss+1;
%     end
% end
% fclose (fp6);
% 
% fp7 = fopen('节点参数2.dat','w');
% sss=1;
% for ppp = 1:Np
%     for jjj = 1:Naver
%         fprintf(fp7,'%f %f\n',sum_Num_All_Cycle(sss,3),sum_Num_All_Cycle(sss,4));
% %         fprintf(fp6,'%d %d %f %f %f %f %f %f\n',ppp,jjj,H_topology(sss),ratio_E_I(sss),ratio_E(sss),ratio_I(sss),ratioDds_E_I(sss),P_Dds_EfromE(sss));
% %         ssT_persistent_c(sss) = T_persistent_c(ppp,jjj);
%         sss=sss+1;
%     end
% end
% fclose (fp7);
% 
%%  单个查看 *******************************************************
clear corrEI_r corrEI_p corrCycle_r corrCycle_p
TatolDuration=[];
TatolDds_fromE=[]; TatolDds_fromI=[]; TatolEI_ratio=[];

for iii=1:190
    index_singel=iii+810;
    analyzy_singel
    corrEI_r(iii,:)=EI_r;corrEI_p(iii,:)=EI_p;
    corrCycle_r(iii,:)=Cycle_r;corrCycle_p(iii,:)=Cycle_p;
end
%%
mean_corrEI=mean(corrEI_r,1);
std_corrEI=std(corrEI_r,0,1);

figure(3);
subplot(311),bar(mean_corrEI);  % 绘制柱形
hold on;
errorbar(1:length(mean_corrEI), mean_corrEI, std_corrEI, 'k', 'linestyle', 'none');

mean_corrCycle=mean(corrCycle_r,1);
std_corrCycle=std(corrCycle_r,0,1);
subplot(312),bar(mean_corrCycle);  % 绘制柱形
hold on;
errorbar(1:length(mean_corrCycle), mean_corrCycle, std_corrCycle, 'k', 'linestyle', 'none');

figure(5);
subplot(311),boxplot(corrEI_p);
subplot(312),boxplot(corrCycle_p);

figure (7)
subplot(221); plot(TatolDds_fromE,TatolDuration,'ro');
subplot(222); plot(TatolDds_fromI,TatolDuration,'ro');
subplot(223); plot(TatolEI_ratio,TatolDuration,'ro');
subplot(4,2,6); histogram(TatolDuration, 'Normalization', 'probability', 'NumBins', 200);
subplot(4,2,8); histogram(TatolEI_ratio, 'Normalization', 'probability', 'NumBins', 200);

%% 输出数据
fp10 = fopen('单细胞相关性_EIvs放电率.dat','w');
for jj = 1:length(mean_corrEI)
     fprintf(fp10,'%d %f %f %f %f\n',jj,mean_corrEI(jj),std_corrEI(jj),mean_corrCycle(jj),std_corrCycle(jj));
end
fclose(fp10);
fp10 = fopen('p值_EIvs放电率.dat','w');
for jj = 1:length(corrEI_p)
     fprintf(fp10,'%d %f %f %f %f %f %f\n',jj,corrEI_p(jj,1),corrEI_p(jj,2),corrEI_p(jj,3),corrCycle_p(jj,1),corrCycle_p(jj,2),corrCycle_p(jj,3));
end
fclose(fp10);

[out_Dura,Dura_histogram] = fun_histogram_one_par(TatolDuration,0,0,0);
[out_E,E_histogram] = fun_histogram_one_par(TatolDds_fromE,0,0,0);
[out_I,I_histogram] = fun_histogram_one_par(TatolDds_fromI,0,0,0);
[out_ratio,ratio_histogram] = fun_histogram_one_par(TatolEI_ratio,0,0,0);
figure (25)
subplot(411), bar(out_Dura,Dura_histogram);
subplot(412), bar(out_E,E_histogram);
subplot(413), bar(out_I,I_histogram);
subplot(414), bar(out_ratio,ratio_histogram);


fp10 = fopen('持续的Duration.dat','w');
for jj = 1:length(TatolDuration)
     fprintf(fp10,'%d %f %f %f %f\n',jj,TatolDds_fromE(jj),TatolDds_fromI(jj),TatolEI_ratio(jj),TatolDuration(jj));
end
fclose(fp10);
fp10 = fopen('E分布.dat','w');
for jj = 1:length(out_E)
     fprintf(fp10,'%f %f\n',out_E(jj),E_histogram(jj));
end
fclose(fp10);
fp10 = fopen('I分布.dat','w');
for jj = 1:length(out_I)
     fprintf(fp10,'%f %f\n',out_I(jj),I_histogram(jj));
end
fclose(fp10);
fp10 = fopen('ratio分布.dat','w');
for jj = 1:length(out_ratio)
     fprintf(fp10,'%f %f\n',out_ratio(jj),ratio_histogram(jj));
end
fclose(fp10);
fp10 = fopen('Duration分布.dat','w');
for jj = 1:length(out_Dura)
     fprintf(fp10,'%f %f\n',out_Dura(jj),Dura_histogram(jj));
end
fclose(fp10);
