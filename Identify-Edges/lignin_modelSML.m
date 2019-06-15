% Copyright © 2013 Cai et al; © 2019 North Carolina State University. All rights reserved.
% “lignin_modelSML” by Cai X, Bazerque JA, Giannakis GB and North Carolina State University is licensed under CC BY 3.0 (https://creativecommons.org/licenses/by/3.0/us/legalcode)


% rng_seed=1010;
% rng(rng_seed) %set seed for random number generator, for reproducibility

%%% Load transcript and protein abundances

load ../TranscriptProteinAbundances.mat % load lignin data

Ytable=Ymeanlines_imputed; % Mean of the experimental replicates
Xmask=Xmeanlines_mask; % indicates targeted transcripts


Y=Ytable{:,:}'; % Rows of Y are the transcripts and proteins, columns are the experiments
X=Xmask{:,:}';

[M,N]=size(Y); % M is number of transcripts and proteins, N is number of experiments

K=~logical(X); % binary matrix size MxN, 1 when transcript/protein was not targeted in an experiment, 0 when it was

[Ysc,Asc]=scaleY(Y); %Scale each transcript and protein so that their abundances fall between 0 and 1


%%set parameters
lambda_factors=10.^(0:-.1:-3); %SML regularization parameters
params.Kcv=5; % # of cross-validation folds
params.rho_factors=10.^(-6:0.1:1); %Ridge regularization parameters
params.lambda_factors=lambda_factors;
params.maxiter=1000; % max iterations for SML algorithm
params.cv_its=5; % # of iterations of the kcv-fold cross-validations

%%network inference     
[ilambda_cv]=cross_validation_SML_B(Y,K,params); % Use cross-validation to determine SML regularization parameter
[rho_factor,sigma2]=cross_validation_ridge_B(Ysc,K,params); % Use cross-validation to determine ridge regularization parameter, and estimate sigma2 for SML

[BRhat,mueRhat]=constrained_ridge_B(Ysc,K,rho_factor); % Calculate ridge regression of Ysc using rho_factor from cross_validation

BR=diag(1./Asc)*BRhat*diag(Asc); % Un-scales calculated BRhat and mueRhat to use with Y
mueR=diag(1./Asc)*mueRhat;

W=1./abs(BRhat); % Use BRhat to compute regularization weights to use with Ysc

Q=Inf; % Ignore discarding rules from SML

% Initialize variables for SML
lambda_factor_prev=1;
BLhat=zeros(M);

for ilambda=1:ilambda_cv
    
    [BLhat,mueLhat]=sparse_maximum_likelihood_B(W,BLhat,Ysc,K,Q,lambda_factors(ilambda),lambda_factor_prev,sigma2,params);    
    
    Q=Inf;
    lambda_factor_prev=lambda_factors(ilambda);
end %ilambda

SL=abs(sign(BLhat)); % Detected edges from SML algorithm
[BChat,mueChat]=constrained_ML_B(BLhat,SL,Ysc,K,sigma2,params); % Use constrainedML algorithm on just the identified edges. This solves for the weights of the identifed edges without the l1-norm penalty

BC=diag(1./Asc)*BChat*diag(Asc); % Un-scale BChat and mueChat to use with Y
mueC=diag(1./Asc)*mueChat;

BCsparse=sparse(BC);
