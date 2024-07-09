% IDOTA M-file of Experiments on images (under reflection) in paper "Point Clouds Matching Based on Discrete Optimal
% Transport"
% Authors: Litao Ma, Wei Bian, and Xiaoping Xue
% Data: 02/08/2020

% IDOTA.m
% ---------------------------
% Usage: --[~,~,Gam] = IDOTA(X,Y,p,q);
%        --[R,S,Gam] = IDOTA(X,Y,p,q,al,br);  
%
% Input:
%       X -- N*d
%       Y -- M*d
%       p -- the density of X
%       q --  the density of Y
%       [al,br] -- the projection domain of S
% Output:
%       R -- Orthogonal transformation matrix
%       S -- Stretch transformation matrix
%       Gam -- transport plan
%
% ---------------------------
% Last modified: 06/11/2024


function [R1,S1,Gam1] = IDOTA(X,Y,p,q,varargin)

%% Determine the number of input parameters
if nargin==4
    al = -1;
    br = 1;
else
    al = varargin{1};
    br = varargin{2};
end

%% initialized parameters 
%S^0, \gamma^0, and tolerate error tol
[N,dim] = size(X);
[M,~]= size(Y);
S0 = eye(dim);
S0 = diag(median([diag(S0)';al*ones(1,dim);br*ones(1,dim)]));
R0 = eye(dim);
Gam0 = diag(p);
tol = 1e-10;

%% solve R
RTmp1 = X'*Gam0*Y*S0;
[U,~,V] = svd(RTmp1);
R1 = V*U';
%% solve S
tildeX = R1*X'; 
for l=1:dim
    tildXl = tildeX(l,:)';
    Yl = Y(:,l)';
    lTmp = tildXl*Yl;
    barX(l) = sum(sum(Gam0.*lTmp));
    Gam_sumi = sum(Gam0,1);
    barY(l) = sum(Gam_sumi.*(Yl.^2));
end
S1 = median([barX./barY;al*ones(size(barX));br*ones(size(barX))]);
S1 = diag(S1);
%% solve Gam
pcost = 2;
options.linprog_niter = 2000;
options.linprog_tol = 1e-12;
options.verbose = 0;
options.pcost = pcost;

%Compute nonregularized OT
[Gam1,~]=solver_ot((R1*X')',(S1*Y')',p,q, options);
iter = 1;
maxiter = 1000;
%% iteration
while (norm(R0-R1)>tol || norm(S0-S1)>tol || norm(Gam0-Gam1)>tol) && iter<maxiter
    S0 = S1;
    Gam0 = Gam1;
    R0 = R1;
    % R
    RTmp1 = X'*Gam0*Y*S0;
    [U,~,V] = svd(RTmp1);
    R1 = V*U';
    % S
    tildeX = R1*X'; 
    for l = 1:dim
        tildXl = tildeX(l,:)';
        Yl = Y(:,l)';
        lTmp = tildXl*Yl;
        barX(l) = sum(sum(Gam1.*lTmp));
        Gam_sumi = sum(Gam1,1);
        barY(l) = sum(Gam_sumi.*(Yl.^2));
    end
    S1 = median([barX./barY;al*ones(size(barX));br*ones(size(barX))]);
    S1 = diag(S1);
    % Gam
    [Gam1,~] = solver_ot((R1*X')',(S1*Y')',p,q, options);
    iter = iter+1;
end


