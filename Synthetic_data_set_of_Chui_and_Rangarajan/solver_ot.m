% solver_ot M-file of Experiments on images (under reflection) in paper "Point Clouds Matching Based on Discrete Optimal
% Transport"
% Authors: Litao Ma, Wei Bian, and Xiaoping Xue
% Data: 02/08/2020

% solver_ot.m
% ---------------------------
% Usage: --[Sigma,~] = solver_ot(X,Y,p,q,options);
%       
% Input:
%       X -- N*d
%       Y -- M*d
%       p -- the density of X
%       q --  the density of Y
%       options -- The structure of the parameters for optimization
%           options.linprog_niter -- The maximum number of iterations;
%           options.linprog_tol -- The allowable error for optimization;
%           options.verbose -- verbosity level, 0=nothing is echoed, 3=all
%                               is echoed;
%           options.pcost -- the norm to construct cost matrix C, i.e. C_{ij}=\|x_i-y_j\|^pcost;
% Output:
%       Sigma -- transport plan
%
% Notice: 
% Solve:
%       min <C,Sigma> 
% subject to
%       Sigma*1 = p 
%       Sigma'*1 = q
%  (0 <= Sigma_i,j <= 1 )
% ---------------------------
% Last modified: 06/11/2024


function [Sigma,err] = solver_ot(X,Y,p,q,options)

N= size(X,1);
M = size(Y,1);

options.null = 0;

Cost = getoptions(options, 'Cost', []);
if isempty(Cost)
    pcost = getoptions(options,'pcost',2);
    % use L^2 Cost  C_{i,j}=|X_i-Y_j|^2
    Cost = (repmat( sum(X'.^2)', [1 M] ) + ...
            repmat( sum(Y'.^2) , [N 1] ) - 2*X*Y').^(pcost/2); 
end

[a,b] = meshgrid(1:N,1:M);
a = a(:);
b = b(:);
[a2,b2] = meshgrid(1:M,1:N);
a2 = a2(:);
b2 = b2(:);
%Sigma(:)
A = [ sparse(a,a+(b-1)*N,ones(N*M,1)); ...    %  Sigma 1
     sparse(a2,b2+(a2-1)*N,ones(N*M,1)); ...    % 1' Sigma 
];    % diag(1 Sigma)Y - Sigma X -W   
 
c = [p;q];
Xmin = zeros(N*M,1);
Xmax = ones(N*M,1);
C = Cost(:);

%% mosek solver
% Setup Mosek variables.

prob.c = C;
prob.a = A;
prob.blc = c;
prob.buc = c;
prob.blx = Xmin;
prob.bux = Xmax;


% Set parameters.
param = [];
% max number of iterations
param.MSK_IPAR_INTPNT_MAX_ITERATIONS = getoptions(options, 'linprog_niter', 100);%算法内默认400
% tolerance, primal
param.MSK_DPAR_INTPNT_TOL_PFEAS = getoptions(options, 'linprog_tol', 1e-12);%算法内默认1e-8
param.MSK_DPAR_INTPNT_TOL_REL_GAP = getoptions(options, 'linprog_tol', 1e-12);%算法内默认1e-8
% verbosity level, 0=nothing is echoed, 3=all is echoed
 verb = getoptions(options, 'verbose', 3);

% Perform the optimization.
[r,res] = mosekopt(['echo(' num2str(verb) ') minimize info'], prob, param);
if r~=0
    warning(['Mosek problem: ' res.rcodestr]);
end
if strcmp(res.sol.bas.prosta,'PRIMAL_AND_DUAL_FEASIBLE')== 0
    warning(['Infeasible problem!:' res.sol.bas.prosta])
end 
err.niter = res.info.MSK_IINF_INTPNT_ITER;
 sol = res.sol.itr; %nterior-point solution.
%sol = res.sol.bas;% x solution in basic solution.
w = sol.xx;
% w=[Sigma(:) Z W];
Sigma = reshape(w(1:N*M), [N M]);

