function [p,t] = refine(p,t,ref)
    % p   : Nx2 array of nodes' XY coordinates
    % t   : Mx3 array of elements as node indices
    % ref : Mx1 logical array; ref(x) = 1 if xth element is to be refined

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
        ref = [0; 1; 0; 0; 0; 1];
    else
        if nargin == 0
            load('coordinates.dat');
            load('elements3.dat');
            load('refine.dat');
            p = coordinates(:,2:3);
            t = elements3(:,2:4);
            ref = false(size(t,1),1);
            ref(refine) = true;
        end
        if size(ref,1) ~= size(t,1)
            error('ref must be an Mx1 array');
        end
    end

    ref = logical(ref);

    %% find edges which are to be split
    nume = size(t,1);
    vect = 1:nume;

    e = [t(:,[1,2]); t(:,[2,3]); t(:,[3,1])];
    [e,j,k] = unique(sort(e,2),'rows'); %% Set of edges
    te = [k(vect), k(vect+nume),k(vect+2*nume)]; %% Edges organized by element

    split = false(size(e,1),1); % Tracks which edges will be split.
    split(te(ref,:)) = true;

    %% For triangles with only two split edges, split all 3 edges instead.
    while any(sum(double(split(te)),2) == 2)
        newsplits = sum(double(split(te)),2)==2;
        split(te(newsplits,:)) = true;
    end

    %% Get logical arrays of elements which will be split
    split2 = sum(double(split(te)),2) == 1; % Elements to be cut in 2
    split4 = sum(double(split(te)),2) == 3; % Elements to be cut in 4

    %% Add new nodes at midpoints of split edges
    nnodes0 = size(p,1);
    newnodes = 0.5 * (p(e(split,1),:) + p(e(split,2),:));
    p = [p; newnodes];

    %% Map indices of split edges to indices of new nodes.
    i = zeros(1,size(e,1));
    i(split) = (1:length(find(split))) + nnodes0;

    tret = t(~(split2|split4),:); % temporary variable for new elements

    %% Split refined triangles in to 4 triangles
    if any(split4)
        tnodes = [t(split4,:) i(te(split4,:))];
        tret = [tret;tnodes(:,1),tnodes(:,4),tnodes(:,6)];
        tret = [tret;tnodes(:,4),tnodes(:,2),tnodes(:,5)];
        tret = [tret;tnodes(:,5),tnodes(:,3),tnodes(:,6)];
        tret = [tret;tnodes(:,4),tnodes(:,5),tnodes(:,6)];
    end

    %% Split triangles with 1 hanging node in to 2 triangles
    newnodes = sum(i(te(split2,:))');
    originalnodes = t(split2,:);
    cnodes = e(split,:);
    for j=1:length(find(split2)) % for each triangle being split
        % find the node opposite the hanging node
        newnode = newnodes(j);
        connected_nodes = cnodes(newnode-nnodes0,:);
        othernode = setdiff(originalnodes(j,:), connected_nodes);
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
