%% Authors: Timothy Costa <costa@math.oregonstate.edu>
%%          Duncan McGregor <mcgregod@onid.oregonstate.edu>
%% Created: 2013-12-11


function [ res ] = residual( f, K, u, coordinates, els )
%RESIDUAL Summary of this function goes here
%   Detailed explanation goes here

%coordinates(i,j) maps global dof i and coordinate j (i.e. j=1 ==> x coord)
% els(i,j) maps element i and local dof j to global dof
%f is the forcing funciton handle
%K is the diffusivity function handle





ORD=3;



sizes=size(els);
numelt=sizes(1);

res=ones(size(els,1),1)*Inf;

% this call calculates || [k grad u \cdot n ] || along each edge and creates
% a vector faceres containing residual contributions from edges
[ edges ] = edgemap( els );
[ faceres ] = edgeres( K, u, coordinates, edges, els );

for j=1:numelt

    ldof=els(j,:); %gets local dof
    lcoord=coordinates(ldof',:); %gets local coordinates of dof
    [wq,xq,yq]=tri_quadcofs(lcoord,ORD); %computes quadrature points
    h=-Inf;

    % calculates h
    for k=1:3
        for l=k+1:3
            d=norm(lcoord(k,:)-lcoord(l,:));
            if d > h
                h=d;
            end
        end
    end


    coefs=([ones(3,1), lcoord]\eye(3))*diag(u(ldof)); %calculates columns c1 + c2 x + c3 y for local basis functions

    grad=sum(coefs(:,2:3).')'; %evaluates total gradient
    div=(grad(1)*(K(xq+h/2,yq)-K(xq-h/2,yq)) + grad(2)*(K(xq, yq+h/2)-K(xq,yq-h/2)))/h; %calculates divergence of gradient

    res(j)=h*sqrt(sum(((div+f(xq,yq)).^2)*diag(wq))) + sqrt(h) * faceres(j); %adds residual

end



end
