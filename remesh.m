%% Authors: Jonathan Perry-Houts <jperryh2@uoregon.edu>,
%%          Duncan McGregor <mcgregod@onid.oregonstate.edu>
%%          Timothy Costa <costa@math.oregonstate.edu>
%% Created: 2013-12-11

% p is the new dof locations
% t is the new set of elements
% err is an a-posteriori error estimate of the initial mesh

%f is a funciton handle for the forcing fucntion of the form f=@(x,y)
%K is a function handle for the local diffusivity of teh form K= @(x,y)
%u is a computed finite element solution
%p is an array dof X 2 where p(i,j) is the jth coordinate of global dof i
%t is a local to global map elements X 3 where t(i,j) is global dof of
%element i's local dof j

function [p,t,err] = remesh(f, K, u, p, t)
    % Calculate residual
    res = residual(f, K, u, p, t);
    err=sum(res);

    % Get locations of elements in the highest 10% residual
    [B,IX] = sort(res, 'descend');
    ref = false(size(res));
    top_ten = IX(1:ceil(length(IX)*0.1));
    ref(top_ten) = true;

    % Refine those elements
    [p,t] = refine(p, t, ref);
end
