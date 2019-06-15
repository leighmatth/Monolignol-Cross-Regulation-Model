% Copyright © 2013 Cai et al; © 2019 North Carolina State University. All rights reserved.
% “lignin_modelSML” by Cai X, Bazerque JA, Giannakis GB and North Carolina State University is licensed under CC BY 3.0 (https://creativecommons.org/licenses/by/3.0/us/legalcode)

function [B,me]=sparse_maximum_likelihood_B(W,B,Y,K,Q,lambda_factor,lambda_factor_prev,sigma2,params)
if(nargin==8)
    params=[];
end
if(isfield(params,'maxiter'))
    maxiter=params.maxiter;
else
    maxiter=1000;
end

[M,N]=size(Y);
[Y,meanY]=centerY(Y,K);

IM=speye(M);

Ylamb=zeros(M);
for i=1:M
    nontargs=K(i,:);
    Yc=Y{i}(:,nontargs)*Y{i}(:,nontargs)';
    Ylamb(i,:)=Yc(i,:);
end

lambda_max=max(max(abs((-sigma2*sum(K,2).*IM+Ylamb)./W)));

lambda=lambda_max*lambda_factor;

Dlambda=(2*lambda_factor-lambda_factor_prev)*lambda_max*W;
S=(abs(Q)>=Dlambda);

for i=1:M/2 % ti cannot be influenced by its associated protein
    S(i,i+M/2)=0;
end

s=sum(S,2);

ei=zeros(M,1);

B=B.*S;
for k=1:maxiter
    Bold=B;
    for i=1:M
        nontargs=K(i,:); % Un-targeted experiments for gene product i
        Nhat=sum(nontargs); % Number of un-targeted experiments
        if (s(i))
            ei(i)=1; % ei is the unit vector
            
            %zi is used and updated below in order to compute the m_ij=cofactor_ij/det(I-B)
            zi=(IM-B)\ei;
            js_i=find(S(i,:)>0);
            for j=js_i
                m_ij=zi(j);
                B_old=B(i,j);
                if(j~=i)
                    lambdaW=lambda*W(i,j);
                    %compute coefs r_ij(=alpha_ij) beta_ij and d_ij:=[det(I-B)+c_ijB(i,j)]/c_ij
                    y_j=Y{i}(j,nontargs)'; % Only want to solve for B that explains the un-targeted experiments
                    
                    BiT=B(i,:);
                    BiT(j)=0; %to remove B(i,j) and obtain the ith row of B\ij;
                    a_iT=(ei'-BiT)*Y{i}(:,nontargs);
                    
                    r_ij=y_j'*y_j;
                    beta_ij=a_iT*y_j;
                    if(m_ij==0) %go to the linear equation
                        Bij=(beta_ij-lambdaW)/r_ij; %linear equation, assume B_ij>0
                        if Bij>0
                            B(i,j)=Bij;
                        else
                            Bij=(beta_ij+lambdaW)/r_ij; %linear equation, assume B_ij<0
                            if Bij<0
                                B(i,j)=Bij;
                            else
                                B(i,j)=0;
                            end %Bij<0
                        end %B_ij>0
                        
                    else%if m_ij ~=0 go to the quadratic equation
                        
                        %now let see which roots are candidates; 0 is always a candidate
                        %quadratic equation, assume Bij>0
                        d_ij=1/m_ij+B(i,j);
                        theta_ijp=r_ij*d_ij+beta_ij-lambdaW;
                        k_ijp=d_ij*(beta_ij-lambdaW)-Nhat*sigma2;
                        q_ijp=theta_ijp^2-4*r_ij*k_ijp;
                        Bijpp=(1/(2*r_ij))*(theta_ijp+sqrt(q_ijp));
                        Bijpm=(1/(2*r_ij))*(theta_ijp-sqrt(q_ijp));
                        
                        %quadratic equation, assume Bij<0
                        %do not compute theta and k again,
                        %use instead (q_ij-)=(q_ij+) -4lambdaW(beta_ij-r_ijd_ij)
                        %theta_ij=r_ij*d_ij+beta_ij+lambdaW;
                        %k_ij=d_ij*(beta_ij+lambdaW)-N*sigma2;
                        q_ijm=q_ijp +4*lambdaW*(beta_ij-r_ij*d_ij);
                        theta_ijm=theta_ijp+2*lambdaW;
                        Bijmm=(1/(2*r_ij))*(theta_ijm-sqrt(q_ijm));
                        Bijmp=(1/(2*r_ij))*(theta_ijm+sqrt(q_ijm));
                        candsBij=[0, Bijmm, Bijmp, Bijpm, Bijpp]';
                        
                        %discard roots that do not satisfy the sign condition
                        candsBij=candsBij.*(([1 -1 -1 1 1]'.*candsBij)>0);
                        Lss=sigma2*Nhat*log((d_ij-candsBij).^2+1e-16   )/2-0.5*sum((ones(length(candsBij),1)*a_iT-candsBij*Y{i}(j,nontargs)).^2,2)-lambdaW*abs(candsBij);
                        [~,l_max]= max(Lss);
                        B(i,j)=candsBij(l_max);
                    end% m_ij==0
                    
                end%if(j~=i)
                dB=B_old-B(i,j);
                zi=zi/(1+dB*m_ij);
                
            end%j
            
            ei(i)=0; % re-set ei back to zero vector
            
        else%s(i)
            B(i,:)=zeros(1,M);
        end%s(i)
    end%i
    
    delta_B=norm(Bold-B,'fro')/(norm(Bold,'fro')+1e-10); % Change in B
    
    if delta_B<1e-4 % If change in B is < 1e-4, then stop solving for B
        for i=1:M
            me(i,1)=(IM(i,:)-B(i,:))*meanY{i}; % Compute mue(i) associated with B(i,:)
        end
        return
    end %delta_B
end %k

for i=1:M
    me(i,1)=(IM(i,:)-B(i,:))*meanY{i}; % Compute mue(i) associated with B(i,:)
end

