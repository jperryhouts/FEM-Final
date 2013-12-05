function [ res ] = residual( f, K, u, coordinates, els )
%RESIDUAL Summary of this function goes here
%   Detailed explanation goes here


ORD=3;



sizes=size(els);
numelt=sizes(1);

res=ones(size(els,1),1)*Inf;


for j=1:numelt
    
    ldof=els(j,:);
    lcoord=coordinates(ldof',:);
    [wq,xq,yq]=tri_quadcofs(lcoord,ORD);
    h=-Inf;
    
    for k=1:3
        for l=k+1:3
            d=norm(lcoord(k,:)-lcoord(l,:));
            if d > h
                h=d;
            end
        end
    end
              

coefs=([ones(3,1), lcoord]\eye(3))*diag(u(ldof));
    linter=[ones(3,1), xq', yq']*coefs;   
    grad=sum(linter(:,2:3).')';
    div=(grad(1)*(K(xq+h/2,yq)-K(xq-h/2,yq)) + grad(2)*(K(xq, yq+h/2)-K(xq,yq-h/2)))/h;
    
    res(j)=h*sqrt(sum(((div+f(xq,yq)).^2)*diag(wq)));
end



end

