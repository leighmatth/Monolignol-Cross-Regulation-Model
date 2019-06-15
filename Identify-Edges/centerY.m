% Copyright © 2013 Cai et al; © 2019 North Carolina State University. All rights reserved.
% “lignin_modelSML” by Cai X, Bazerque JA, Giannakis GB and North Carolina State University is licensed under CC BY 3.0 (https://creativecommons.org/licenses/by/3.0/us/legalcode)

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
