%% Authors: Duncan McGregor <mcgregod@onid.oregonstate.edu>
%% Created: 2013-12-11


function edges = edgemap(t)
nel=size(t,1);
k=1;
for i=1:nel
    for j=i+1:nel
        int=intersect(t(i,:),t(j,:));
        if size(int,2)==2
            edges(k,1)=i;
            edges(k,2)=j;
            edges(k,3)=int(1);
            edges(k,4)=int(2);
            k=k+1;
        end
    end
end