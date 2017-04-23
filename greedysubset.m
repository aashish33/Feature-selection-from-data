% Input: number of features F
%       matrix X of features, with n rows (samples), d columns (features)
%       X(i,j) is the j-th feature of the i-th sample
%       vector y of scalar values, with n rows (samples), 1 column
%       y(i) is the scalar value of the i-th sample
% Output: vector of selected features S, with F rows, 1 column
%       vector thetaS, with F rows, 1 column
%       thetaS(1) corresponds to the weight of feature S(1)
%       thetaS(2) corresponds to the weight of feature S(2)
%       and so on and so forth
function [S thetaS] = greedysubset(F,X,y)
[n d] = size(X);
 l = [];
for p=1:d
    l = [l p];
end
S = [];
thetaS = [];
for f = 1:F
    A = setdiff(l,S);
    x = [];
    J = [];
    xprev=[];
    T = [];
     for e = S
         xprev = [xprev  X(:,[e])];
     end
     for j = A
        x = xprev;
        x = [x  X(:,[j])];
        thetasJ = pinv(x)*y;

        w = 0;
        for t=1:n
            xts = x([t],:);
            w = w + ((y(t) - dot(thetasJ,xts))^2);
        end
        J = [J 1/2*w];
        T = [T thetasJ];
     end
    m = min(J);
    index = find(J==m);
    jhat = A(index);
    S = [S jhat];
    thetaS = [ thetaS T(index)];
end

S = transpose(S);
thetaS = transpose(thetaS);
        
        