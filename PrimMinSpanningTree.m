function [T] = PrimMinSpanningTree(G,v_id)
% check if vertices have names
if (~sum(ismember(G.Nodes.Properties.VariableNames,'Name')))
    % if not, give names using its indices
    Vnames = int2str(1:numnodes(G));
    G.Nodes.Name = split(Vnames);
end

% check if edges have names
if (~sum(ismember(G.Edges.Properties.VariableNames,'Name')))
    % if not, give names using its indices
    Enames = int2str(1:numedges(G));
    G.Edges.Name = split(Enames);
end


% set the dfnumber of all vertices to -inf
G.Nodes.discN = -inf(numnodes(G),1);

% Let T = the vertex with id v_id
T = graph;
T = addnode(T,1);

% record the original id for the vertex in G in the origId attribute of the
% nodes of T
T.Nodes.origId(1) = v_id;
T.Nodes.Name(1) = G.Nodes.Name(v_id);


% initiate the counting of dfnumber
currentDf = 0;

% set the dfnumber for the starting vertex in G and in T
G.Nodes.discN(v_id) = currentDf;
T.Nodes.discN(1) = currentDf;

% the first set of frontier edges are the edges from the vertex with id
% v_id
[S,nV] = outedges(G,v_id);

S = S(nV~=v_id);


% when the set of frontier edges is not empty
while ~isempty(S)
    % advance the dfnumber
    currentDf = currentDf+1;
    
    % edge idx (in G) of the next edge
    eidx = PrimNextEdge(G,S);
    
    % endpints of the chosen next edge
    endpts = G.Edges.EndNodes(eidx,:);
    endpts = findnode(G,{endpts{1} endpts{2}});
    
    % next vertex and its tree node
    if (~isinf(G.Nodes.discN(endpts(1)))) % endpts(1) is a tree node
        w_id = endpts(2);
        pre_id = endpts(1);
    else % endpts(2) is a tree node
        w_id = endpts(1);
        pre_id = endpts(2);
    end
    % add the new node to the tree
    newNode = table(G.Nodes.Name(w_id), w_id, currentDf,'VariableNames', {'Name','origId', 'discN'});
    T = addnode(T,newNode);
    
    % create the edge and its attributes (endpts and original id in G) to be added in T
    newEdge = table([G.Nodes.discN(pre_id)+1,currentDf+1],G.Edges.Name(eidx),eidx,G.Edges.Weight(eidx), 'VariableNames', {'EndNodes','Name', 'origId', 'Weight'});
    T = addedge(T,newEdge);
    
    % record the dfnumber in G and T for the next discovered vertex
    G.Nodes.discN(w_id) = currentDf;
    
    % update the set of frontier edges
    S = PrimFrontierEdge(G,S);
    
end
end % end function PrimMinSpanningTree


%%Next Edge

function [eidx] = PrimNextEdge(G,S)

    %First, sort S from small to big
    sort_S = sort(S,'Ascend');

    %run a loop across the length of S
    for i=1:length(sort_S)

        %extract the weight and position of the Edge with the minimum weight in S
        [weight, pos] = min(G.Edges.Weight(S));
        
        %set next edge to be the edge at that position.
        eidx = S(pos);
    end

end %end function dfsNextEdge+

%% Frontier Edge

function S_new = PrimFrontierEdge(G,S)
    %First, sort S from small to big
    sort_S = sort(S,'Ascend');
    S_new = [];
    %Start a loop across the length of Nodes of G.
    for i = 1:length(G.Nodes.Name)
        
        % find all edges spanning from the discovered nodes.
        %determine the discovered nodes by its dfN, if it is not -Inf, then it is discovered.
        if (~ isinf(G.Nodes.discN(i)))

            % Set S_new as the frontier edges from that node. However, this
            % S includes both discovred and undiscovered edges
            S1 = outedges(G,G.Nodes.Name(i));
            S_new = cat(2,S_new,S1');
        end
    end

    %Now filter out the discovered edges
    %set an empty discoveredS for later
    discoveredS = [];

    %start a loop across the length of the S_new found above
    for k = 1:length(S_new)
        %locate the endpoints of each edge at each k
        endpoints = G.Edges.EndNodes(S_new(k),:);
        endpoints = findnode(G,{endpoints{1} endpoints{2}});

        %if the dfN of both of the endpoints of that k-th edge is greater
        %or equal to 0, then that edges is discovred.
        if (G.Nodes.discN(endpoints(1))>=0) && (G.Nodes.discN(endpoints(2))>=0)

            %call the k-th discovered edge is newS
            newS = (S_new(k));

            %append all the discovered edges to discoveredS
            discoveredS(end+1) = [newS];

        end
    end

    %use setdiff to extract the undiscovered from S_new vs the discovredS.
    S_new = setdiff(S_new, discoveredS);
end


