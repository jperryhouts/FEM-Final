%% Authors:   Timothy Costa <costa@math.oregonstate.edu>
%%            Duncan McGregor <mcgregod@math.oregonstate.edu>
%% Created: 2013-12-11


function [ faceres ] = edgeres( K, u,coordinates,els, edges )

% edgres calculates the contribution to the a
% posteriori error estimation from edges..

% input: K - coefficient vector
% u: fem solution
% edges - a #edges x 4 array containing element and point
%   information for edges
%   edges(i,1) and edges(i,2) give global element number
%     for elements stradling edge i
%   edges(i,3) and edges(i,4) give global numbers for points
%     on the edge

% output: faceres(j) is the

numedge = size(edges(:,1));
faceres=zeros(size(els,1),1);

% loop over edges
for j = 1:numedge
  % calculate k grad u \cdot n on elements edges(j,1) and edges(j,2)
  ldof1=els(edges(j,1),:);
  ldof2=els(edges(j,2),:);
  p1=coordinates(ldof1,:);
  p2=coordinates(ldof2,:);

  ind1(1)=find(ldof1==edges(j,3));
  ind1(2)=find(ldof1==edges(j,4));
  ind1(3)=find(ldof1~=edges(j,3)&&ldof1~=edges(j,4));

  [Q,R]=qr([p1(ind1(2),:).'-p1(ind1(1),:).', p1(ind1(3),:).'-p1(ind1(1),:).']);
  n1=-Q(:,2);
  n2=-n1;

  coefs1 = ([ones(3,1), p1]\eye(3))*diag(u(ldof1)); %get local basis functions
  coefs2 = ([ones(3,1), p2]\eye(3))*diag(u(ldof2)); %get local basis functions
  
  grad1=sum(coefs1(:,2:3).'); %Compute local gradient on element 1
  grad2=sum(coefs2(:,2:3).'); %Compute local gradient on element 2

  m1=.5*(p1(ind(1),:)+p1(ind(2),:)); %compute midpoint of edge
  s=norm(p1(ind(1),:)-p2(ind(2),:)); %computes side length for L2 error
  

  kdivun1 = grad1*n1; %compute normal derivatives
  kdivun2 = grad2*n2; 
  % calculate || kdivu2 - kdivu1 ||
  facediff = sqrt(K(m1(1),m1(2))^2*(kdivun2 - kdivun1)^2*s); %midpoint approximation of L2 error
  % add to element sized edge residual vector
  faceres(edges(j,1)) = faceres(edges(j,1)) + 0.5 * facediff;
  faceres(edges(j,2)) = faceres(edges(j,2)) + 0.5 * facediff;

end

end