function ilambda_ms=cross_validation_SML_B(Yorig,Korig,params)
if(nargin==2)
    params=[];
end
if(isfield(params,'Kcv'))
    Kcv=params.Kcv;
else
    Kcv=5;
end

if(isfield(params,'lambda_factors'))
    lambda_factors=params.lambda_factors;
else
    lambda_factors=10.^(0:-.2:-3);
end

if(isfield(params,'cv_its'))
    cv_its=params.cv_its;
else
    cv_its=5;
end

[M,N]=size(Yorig);
Ntest=floor(N/Kcv); % Number of experiments in testing set
Nlambdas=length(lambda_factors);

for its=1:cv_its % do this for cv_its random permutations of the original data, leading to different learning and training sets. ilambda from the cv_iteration with the smallest error is chosen
    
    % For each cv_iteration re-initialize the variables
    Errs=zeros(Nlambdas,Kcv);
    ilambda=0;
    
    % Randomly permute the order of the experiments
    perm=randperm(N);
    Y=Yorig(:,perm);
    K=Korig(:,perm);
    
    lambda_factor_prev=1;
    while((ilambda<Nlambdas))
        ilambda=ilambda+1; % lambda_factor index
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
            
            if(ilambda==1) 
                % if first ilambda, need to estimate rho_factor and sigma2 
                % of the learning set using cross-validation ridge regression
                % and use constrained ridge regression on the learning set
                % to estimate BRidghat to calculate the weights W for the
                % SML algorithm
                [rho_factor,sigma2learnt{cv}]=cross_validation_ridge_B(Ylearn,Klearn,params);
                [BRidgehat{cv},mueRidgehat{cv}]=constrained_ridge_B(Ylearn,Klearn,rho_factor);
                
                W{cv}=1./abs(BRidgehat{cv}); % regularization penalty weights for SML algorithm
                
                Q{cv}=Inf; % No discarding rules
                BLs{cv}=sparse(zeros(M));
            end
                        
            [BLhat,mueLhat]=sparse_maximum_likelihood_B(W{cv},BLs{cv},Ylearn,Klearn,Q{cv},lambda_factors(ilambda),lambda_factor_prev,sigma2learnt{cv},params);
            BLs{cv}=BLhat;
            mueLs{cv}=mueLhat;
            SL=abs(sign(BLhat));
            Q{cv}=Inf;
            [BChat,mueChat]=constrained_ML_B(BLhat,SL,Ylearn,Klearn,sigma2learnt{cv},params);
            
            % Un-scale BChat and mueChat to use with Ytest
            BC=diag(1./Alearn)*BChat*diag(Alearn);
            mueC=diag(1./Alearn)*mueChat;
            
            % Calculate Ypred using the testing set
            for j=1:Ntest
                Ypred(:,j)=diag(Ktest(:,j))*BC*Ytest(:,j)+diag(Ktest(:,j))*mueC+diag(~Ktest(:,j))*Ytest(:,j);                
            end
            
            % Calculate the error in the testing set predictions
            for i=1:M
                nontargs=Ktest(i,:); % only calculate error for the experiments where i was not-targeted; in the experiments where i was targeted, Ypred(i,target)=Ytest(i,target)
                err(i,1)=norm(Ytest(i,nontargs)-Ypred(i,nontargs))^2;
            end
            Errs(ilambda,cv)=sum(err);
           
        end %cv
        
        lambda_factor_prev=lambda_factors(ilambda);
        
    end %while ilambda
    
    errs_mean=sum(Errs,2)/Kcv; % Mean error over cv-folds for each ilambda
    [errs_min(its,1),ilambda_min(its,1)]=min(errs_mean); % minimum of the mean errors and the associated ilambda
end

[~,its_min]=min(errs_min); % cv_iteration with the minimum testing error over the cv-folds
ilambda_ms=ilambda_min(its_min); % ilambda associated with that minimum error




