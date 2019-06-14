function [B,me]=constrained_ML_B(B,S,Y,K,sigma2,params)

if(nargin==5)
    params=[];
end

if(isfield(params,'maxiter'))
    maxiter=params.maxiter;
else
    maxiter=1000;
end

[M,~]=size(Y);
IM=speye(M);

[Y,meanY]=centerY(Y,K);

%the corresponding one in entry i of ei will be added in the loop
ei=zeros(M,1);
s=sum(S,2);
I=speye(M);

for k=1:(maxiter)
    Bold=B;
    for i=1:M
        nontargs=K(i,:); % Un-targeted experiments for gene product i
        Nhat=sum(nontargs); % Number of un-targeted experiments
        if (s(i))
            ei(i)=1; % ei is the unit vector
            
            %zi is used and updated below in order to compute the m_ij=cofactor_ij/det(I-B)
            zi=(I-B)\ei;
            js_i=find(S(i,:)>0);
            for j=js_i
                m_ij=zi(j);
                B_old=B(i,j);
                y_j=Y{i}(j,nontargs)';
                
                BiT=B(i,:);
                BiT(j)=0; %to remove B(i,j) and obtain the ith row of B\ij;
                a_iT=(ei'-BiT)*Y{i}(:,nontargs);
                
                r_ij=y_j'*y_j;
                beta_ij=a_iT*y_j;
                if (m_ij==0) %go to the linear equation
                    B(i,j)=beta_ij/r_ij;
                else%if m_ij ~=0 go to the quadratic equation
                    d_ij=1/m_ij+B(i,j);
                    theta_ij=r_ij*d_ij+beta_ij;
                    k_ij=d_ij*beta_ij-Nhat*sigma2;
                    q_ij=theta_ij^2-4*r_ij*k_ij;
                    Bijp=(1/(2*r_ij))*(theta_ij+sqrt(q_ij));
                    Bijm=(1/(2*r_ij))*(theta_ij-sqrt(q_ij));
                    candsBij=[Bijm, Bijp]';
                    Lss=sigma2*Nhat*log((d_ij-candsBij).^2)/2-0.5*sum((ones(length(candsBij),1)*a_iT-candsBij*Y{i}(j,nontargs)).^2,2);
                    [~,l_max]= max(Lss);
                    B(i,j)=candsBij(l_max);
                end %m_ij
                dB=B_old-B(i,j);
                zi=zi/(1+dB*m_ij);
            end%j
            %now coordinate ascent on matrix of perturbation-gene links F
            ei(i)=0;
        end%s(i)
    end%i
    
    delta_B=norm(Bold-B,'fro')/(norm(Bold,'fro')+1e-10); % Change in B
    if delta_B<1e-2 % If change in B is < 1e-2, then stop solving for B
        for i=1:M
            me(i,1)=(IM(i,:)-B(i,:))*meanY{i}; % Compute mue(i) associated with B(i,:)
        end
        return
    end
end %k
for i=1:M
    me(i,1)=(IM(i,:)-B(i,:))*meanY{i}; % Compute mue(i) associated with B(i,:)
end
