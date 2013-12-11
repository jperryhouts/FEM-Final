%% Authors: Jonathan Perry-Houts <jperryh2@uoregon.edu>,
%%          Duncan McGregor <mcgregod@onid.oregonstate.edu>
%%          Timothy Costa <costa@math.oregonstate.edu>
%% Created: 2013-12-11

function [p,t] = remesh(f, K, u, p, t)
    % Calculate residual
    res = residual(f, K, u, p, t);

    % Get locations of elements in the highest 10% residual
    [B,IX] = sort(res, 'descend');
    ref = false(size(res));
    top_ten = IX(1:ceil(length(IX)*0.1));
    ref(top_ten) = true;

    % Refine those elements
    [p,t] = refine(p, t, ref);
end
