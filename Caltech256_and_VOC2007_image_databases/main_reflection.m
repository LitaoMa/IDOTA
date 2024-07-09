% Main M-file of Experiments on images (under reflection) in paper "Point Clouds Matching Based on Discrete Optimal
% Transport"
% Authors: Litao Ma, Wei Bian, and Xiaoping Xue
% Data: 09/15/2021

% main_reflection.m
% ---------------------------
% Usage: [] = main_reflection()

% Notes: In this demo, MOSEK, an optimization solver, will be utilized. 
% Please ensure that the relevant software is installed in advance, 
% and the official download website is: https://www.mosek.com
% ---------------------------
% Last modified: 06/11/2024

clear all
close all
%% Initializing parameters
match_err = []; % Matching error
Time = []; % CPU time

%% load data --car
load('./data/car/data_car_match.mat')% variables:X2,X5,X7,X9,X11,X19
nam={'2','5','7','9','11','19'};

%% reflection matrix
M1 = [-1 0;0 1];
M2 = [1 0; 0 -1];
M3 = [-1 0; 0 -1];
M4 = [0 1; 1 0];
M5 = [0 -1; -1 0];

%% Obtain reflected images
for in = 1:length(nam)
    Ydata=[];
    data = eval(['X',nam{in}]);   
    Ydata(:,:,1) = (M1*data')';
    Ydata(:,:,2) = (M2*data')';
    Ydata(:,:,3) = (M3*data')';
    Ydata(:,:,4) = (M4*data')';
    Ydata(:,:,5) = (M5*data')';
    [N,dim] = size(data);
    % Determine the probability distribution
    p = 1/N*ones(N,1);
    q = p;
    for jn = 1:5
        data1 = Ydata(:,:,jn);
        tic;
%% IDOTA
        [~,~,Gam1] = IDOTA(data,data1,p,q);
        Time(in,jn) = toc;
        % Get pointwise correspondence
        [~,perm] = max(Gam1,[],2);
        match_err(in,jn) = plot_match_result(data,data1,perm,0);    
    end
    mean_err = mean(match_err(in,:));
    mean_time = mean(Time(in,:));
end
mean(match_err,2)
mean(Time,2)

