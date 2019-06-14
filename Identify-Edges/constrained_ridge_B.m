function [B,me]=constrained_ridge_B(Y,K,rho_factor)

[M,~]=size(Y);
[Y,meanY]=centerY(Y,K); % Center Y around 0. % returns cells Y and meanY

B=zeros(M);
IM=speye(M);

    for i=1:M
        if i<=M/2
            rows_rem=[i i+M/2]; % if transcript, the associated protein can not regulate
            I=speye(M-2);
        else
            rows_rem=i;
            I=speye(M-1);
        end
        
        nontargs=K(i,:); % Experiments where i is not targeted
        Nhat=sum(nontargs);
        yi=Y{i}(i,nontargs)';
        Yi=Y{i}(:,nontargs)';

        Yi(:,rows_rem)=[];
        
        Pi=eye(Nhat);
        rho=rho_factor*norm(Yi'*Pi)^2;
        bi=(Yi'*Pi*Yi+rho*I)\(Yi'*(Pi*yi));

        indexes=1:M;
        indexes(rows_rem)=[];
        B(i,indexes)=bi';
        
        me(i,1)=(IM(i,:)-B(i,:))*meanY{i}; % Solve for mue

    end
   