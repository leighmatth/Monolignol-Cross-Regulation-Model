function [Ysc, Ascaling]=scaleY(Y)

%Y is an M x N matrix with non-negative elements

%Scale Data so each row has the range [0 1]
    Ymax=max(Y,[],2);
    Ascaling=1./Ymax;
    Ysc=diag(Ascaling)*Y;
    
end