% Copyright © 2013 Cai et al; © 2019 North Carolina State University. All rights reserved.
% “lignin_modelSML” by Cai X, Bazerque JA, Giannakis GB and North Carolina State University is licensed under CC BY 3.0 (https://creativecommons.org/licenses/by/3.0/us/legalcode)

function [rho_factor_m,sigma2ridge]=cross_validation_ridge_B(Yorig,Korig,params)

if(isfield(params,'rho_factors'))
    rho_factors=params.rho_factors;
else
    rho_factors=10.^(-6:0.2:1);
end

if(isfield(params,'Kcv'))
    Kcv=params.Kcv;
else
    Kcv=5;
end

if(isfield(params,'cv_its'))
    cv_its=params.cv_its;
else
    cv_its=5;
end

[M,N]=size(Yorig);
Nrho=length(rho_factors);
Ntest=floor(N/Kcv); %number of experiments in testing set

N_avg=sum(sum(Korig))/M; % Average number of un-targeted experiments per transcript/protein
Ntest_avg=N_avg/Kcv; % Average number of un-targeted experiments in the testing data set

for its=1:cv_its % do this for cv_its random permutations of the original data, leading to different learning and training sets. rho_factor_m and sigma2ridge from the cv_iteration with the smallest error is chosen
    
    % For each cv_iteration re-initialize the variables
    Errs=zeros(Nrho,Kcv);
    min_attained=0;
    irho=0;
    err_mean=0;
    
    % Randomly permute the order of the experiments
    perm=randperm(N);
    Y=Yorig(:,perm);
    K=Korig(:,perm);
    
    while((irho<Nrho)&&(~min_attained))
        irho=irho+1; % rho_factor index
        for cv=1:Kcv
            
            test_indices=(Ntest*(cv-1)+1):(Ntest*cv);
            learn_indices=setdiff(1:N,test_indices);
            
            
            Ylearn=Y(:,learn_indices);
            Klearn=K(:,learn_indices);
            Ytest=Y(:,test_indices);
            Ktest=K(:,test_indices);
            
            Ypred=zeros(M,length(test_indices)); %re-initialize Ypred for each cv
            err=zeros(M,1); %re-initialize testing error
            
            [Ylearn,Alearn]=scaleY(Ylearn);
            
            [BRhat,mueRhat]=constrained_ridge_B(Ylearn,Klearn,rho_factors(irho));
            
            % Un-scale BRhat and mueRhat to use with Ytest
            BR=diag(1./Alearn)*BRhat*diag(Alearn);
            mueR=diag(1./Alearn)*mueRhat;
            
            % Calculate Ypred using the testing set
            for j=1:Ntest
                Ypred(:,j)=diag(Ktest(:,j))*BR*Ytest(:,j)+diag(Ktest(:,j))*mueR+diag(~Ktest(:,j))*Ytest(:,j);
            end
            
            % Calculate the error in the testing set predictions
            for i=1:M
                nontargs=Ktest(i,:); % only calculate error for the experiments where i was not-targeted; in the experiments where i was targeted, Ypred(i,target)=Ytest(i,target)
                err(i,1)=norm(Ytest(i,nontargs)-Ypred(i,nontargs))^2;
            end
            Errs(irho,cv)=sum(err);
                        
        end%cv
        
        %check if the error starts to increase and then stop scanning for irhos
        err_mean_prev=err_mean;
        err_mean=mean(Errs(irho,:),2);
        if (irho>1)
            if (err_mean_prev<=err_mean)
                min_attained=1; %the minimium was attained for the previous irho
                irho=irho-1;
            end
        end
    end %while irho
    
    errs_mean(its,1)=mean(Errs(irho,:),2);
    irho_min(its,1)=irho;
    sigma2hat(its,1)=errs_mean(its,1)/(Ntest_avg*M-1);
end

[~,min_its]=min(errs_mean); % cv_iteration with the minimum testing error over the cv-folds

rho_factor_m=rho_factors(irho_min(min_its))*N_avg/(N_avg-Ntest_avg); % rho_factor associated with that minimum error
sigma2ridge=sigma2hat(min_its); % estimated sigma2

