% Computing matching result M-file of Experiments on synthetic data set of Chui and Rangarajan (under rotation,stretch, noise) in paper "Point Clouds Matching Based on Discrete Optimal
% Transport"
% Authors: Litao Ma, Wei Bian, and Xiaoping Xue
% Data: 09/16/2021

%  TPS_match_result.m
% ---------------------------
% Usage: error = TPS_match_result(X,Y,pind)
%     
% Input:
%       X -- N*d;
%       Y -- M*d;
%       pind -- the corresponding index of Y to X;
% Output:
%       error -- the error after Thin-Plate smoothing Spline transform
%
% ---------------------------
% Last modified: 06/11/2024

function err_background = TPS_match_result(X,Y,pind)

ST = tpaps(X',Y(pind,:)');
Xtran = transpose(fnval(ST,X'));
errdif = Xtran-Y(1:length(pind),:);
err_background = mean(sqrt(sum(errdif.^2,2)));
