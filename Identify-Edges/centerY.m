function [Yc,Y_bar]=centerY(Y,K)

[M,~]=size(Y);

Y_bar=cell(M,1);
Yc=cell(M,1);
for i=1:M
    nontargs=K(i,:);
    Nhat=sum(nontargs);
    Y_bar{i}=sum(Y(:,nontargs),2)/Nhat;
    Yc{i}=Y-Y_bar{i};
end

end