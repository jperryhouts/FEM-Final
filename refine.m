function [p,t] = refine(p,t,tri)
    %  p   : Nx2 array of nodal XY coordinates, [x1,y1; x2,y2; etc]
    %  t   : Mx3 array of triangles as indices, [n11,n12,n13; n21,n22,n23; etc]
    %  tri : Mx1 logical array, with tri(k) = TRUE if kth triangle is to be refined

    DEBUG = false;

    if DEBUG
        %% Test case on unit square.
        p = [0 0;   .5 0;   1 0; ...
            0 .5;           1 .5; ...
            0 1;    .5 1;   1 1];
        t = [1 2 4; ...
            4 2 7; ...
            6 4 7; ...
            2 3 5; ...
            2 5 7; ...
            7 5 8];
        tri = [0; 1; 0; 0; 0; 1];
    else
        if nargin ~= 3
            error('Wrong number of inputs');
        end
        if numel(tri) ~= size(t,1)
            error('tri must be an Mx1 array');
        end
    end

    tri = logical(tri);

    %% find edges which are to be split
    numt = size(t,1);
    vect = 1:numt;

    e = [t(:,[1,2]); t(:,[2,3]); t(:,[3,1])];
    [e,j,k] = unique(sort(e,2),'rows'); %% set of edges
    te = [k(vect), k(vect+numt),k(vect+2*numt)]; %% Edges in each triangle

    split = false(size(e,1),1);
    split(te(tri,:)) = true;

    %% For triangles with only two split edges, split all 3 edges instead.
    while any(sum(double(split(te)),2) == 2)
        newsplits = sum(double(split(te)),2)==2;
        split(te(newsplits,:)) = true;
    end

    %% Get set of triangles which will be split
    split1 = sum(double(split(te)),2) == 1;
    split3 = sum(double(split(te)),2) == 3;

    %% Add new nodes
    nnodes0 = size(p,1);
    newnodes = 0.5 * (p(e(split,1),:) + p(e(split,2),:));
    p = [p; newnodes];

    %% Map indices of split edges to indices of new nodes.
    i = zeros(1,size(e,1));
    i(split) = (1:length(find(split))) + nnodes0;

    tret = t(~(split1|split3),:); % temporary variable for new triangles

    %% Split refined triangles in to 4 triangles
    if any(split3)
        tnodes = [t(split3,:) i(te(split3,:))];
        tret = [tret;tnodes(:,1),tnodes(:,4),tnodes(:,6)];
        tret = [tret;tnodes(:,4),tnodes(:,2),tnodes(:,5)];
        tret = [tret;tnodes(:,5),tnodes(:,3),tnodes(:,6)];
        tret = [tret;tnodes(:,4),tnodes(:,5),tnodes(:,6)];
    end

    %% Split triangles with 1 hanging node in to 2 triangles
    newnodes = sum(i(te(split1,:))');
    othernodes = t(split1,:);
    cnodes = e(split,:);
    for j=1:length(find(split1)) % for each triangle being split
        % find the node opposite the new hanging node
        newnode = newnodes(j);
        connected_nodes = cnodes(newnode-nnodes0,:);
        othernode = setdiff(othernodes(j,:), connected_nodes);
        tret = [tret; newnode, othernode, connected_nodes(1)];
        tret = [tret; newnode, connected_nodes(2), othernode];
    end

    t = tret;

    if DEBUG
        showmesh(p,t)
    end
end

function showmesh(p,t)
    clf
    plot(p(:,1),p(:,2),'b.')
    hold on
    e = [t(:,[1,2]); t(:,[2,3]); t(:,[3,1])];
    e = unique(sort(e,2),'rows');
    for j=1:size(e,1)
        X = [p(e(j,1),1) p(e(j,2),1)];
        Y = [p(e(j,1),2) p(e(j,2),2)];
        plot(X, Y)
        %text(sum(X)/2, sum(Y)/2, ['e',num2str(j)])
    end
    for i=1:size(p,1)
        text(p(i,1), p(i,2), ['n',num2str(i)])
    end
    axis equal off
end
