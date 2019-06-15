% Copyright © 2013 Cai et al; © 2019 North Carolina State University. All rights reserved.
% “lignin_modelSML” by Cai X, Bazerque JA, Giannakis GB and North Carolina State University is licensed under CC BY 3.0 (https://creativecommons.org/licenses/by/3.0/us/legalcode)

function [Ysc, Ascaling]=scaleY(Y)

%Y is an M x N matrix with non-negative elements

%Scale Data so each row has the range [0 1]
    Ymax=max(Y,[],2);
    Ascaling=1./Ymax;
    Ysc=diag(Ascaling)*Y;
    
end
