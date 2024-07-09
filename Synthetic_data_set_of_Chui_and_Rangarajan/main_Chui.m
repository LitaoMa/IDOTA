% Main M-file of Experiments on synthetic data set of Chui and Rangarajan (under rotation,stretch, noise) in paper "Point Clouds Matching Based on Discrete Optimal
% Transport"
% Authors: Litao Ma, Wei Bian, and Xiaoping Xue
% Data: 09/22/2021

% main_Chui.m
% ---------------------------
% Usage: [] = main_Chui()

% Notes: In this demo, MOSEK, an optimization solver, will be utilized. 
% Please ensure that the relevant solver is installed in advance, 
% and the official download website is: https://www.mosek.com
% ---------------------------
% Last modified: 06/11/2024

close all
clear all
%% --------Load data--------------
name1 = './data/';
%% rotation
name2 ='fish_rotateTest_guas0d5_deform0d03_rotate';
name3 = {'30','60','90','120','150','180'};
X0 = load(['./data/templete_fish.mat']);% variable: x1

%% stretch
% name2 ='fish_rotate_guass0d3_stretchxy_';
% name3 = [-4,-2,-1,-0.5,-0.2,0.2,0.5,1,2,4];
% X0 = load(['./data/templete_fish.mat']);% variable: x1

%% noise
% name2 ='fish_stretch2_rotate_noise0d0';
% name3 = [0.01:0.01:0.05];
% X0 = load(['./data/templete_fish.mat']);% variable: x1

X = X0.x1;  % initial data
%% ---------loop--------------
for idfm =1 % 1:length(name3)
    %% rotation 
    field = [name1,name2,name3{idfm},'.mat'];
    %% stretch
%     if name3(idfm)<0
%         str = 'neg';
%     else
%         str = 'pos';
%     end
%     nam = name3(idfm);
%     if abs(nam)<1
%         nam = ['0d',num2str(abs(nam)*10)];
%     else
%         nam = num2str(abs(nam));
%     end
%     field =[name1,name2,str,nam,'.mat'];
    %% noise
%     field =[name1,name2,num2str(name3(idfm)*100),'.mat'];
       
    Y0 = load(field);    %100 random trials. variable:Y size:98*2*100(fish) or 105*2*100(chinese)
    for jtr = 1:100
        Y = Y0.Y(:,:,jtr);% stretch or noise, 
        %% Get density
        [N,dim] = size(X);
        [M,~] = size(Y);
        p = 1/N*ones(N,1);
        q = 1/M*ones(M,1);
        %% --------our method-----------
        tic;
        [~,~,Gam] = IDOTA(X,Y,p,q);
        time = toc;
        TotalTime(idfm,jtr) = time;  
        %% Determine the final outcome
        [~,perm] = max(Gam,[],2);%one to many
        match_err(idfm,jtr) = plot_match_result(X,Y,perm,0);% 2021.8.12，拟合误差，此方法计算误差更低
    end
    disp(['mean affine err is:',num2str(mean(match_err(idfm,:)))])
end
disp('mean time and err are:')
disp(mean(mean(TotalTime)))
disp(mean(match_err,2))
