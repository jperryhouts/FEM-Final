%% Authors: Jonathan Perry-Houts <jperryh2@uoregon.edu>,
%%          Duncan McGregor <mcgregod@onid.oregonstate.edu>
%%          Timothy Costa <costa@math.oregonstate.edu>
%% Created: 2013-12-11

function [p,t] = remesh(f, K, u, p, t)
    res = residual(f, K, u, p, t);
    [p,t] = refine(p, t, res);
end
