function [T] = bfsSpanningTree(G,v_id)

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
eidx = bfsNextEdge(G,S);
    
    endpts = G.Edges.EndNodes(eidx,:);
    %endpts = cellfun(@str2num, endpts);
    endpts = findnode(G,endpts);
    w_id = [];
    for i = 1:(size(endpts,1)*size(endpts,2))
        if (isinf(G.Nodes.discN(endpts(i))))
            newN = endpts(i);
            w_id(end+1) = [newN];
            
        else
            pre_id = endpts(i);
        end
        
    end
    w_id = unique(w_id);
    for i = 1:length(w_id)
        G.Nodes.discN(w_id(i)) = currentDf +i;
    end
    
    

    % add the new node to the tree
    newNode = table(G.Nodes.Name(w_id), w_id', G.Nodes.discN(w_id),'VariableNames', {'Name','origId', 'discN'});
    T = addnode(T,newNode);

    % create the edge and its attributes (endpts and original id in G) to be added in T
    newEdge = table(G.Edges.EndNodes(eidx,:),G.Edges.Name(eidx),eidx','VariableNames', {'EndNodes','Name','origId'});
    T = addedge(T,newEdge);
     S = bfsupdateFrontierEdge(G,S);
    currentDf = max(G.Nodes.discN);
while ~isempty(S)
%     currentDf = currentDf+1;
    eidx = bfsNextEdge(G,S);
    S_test = [];
for i = 1:length(eidx)
    endpoints = G.Edges.EndNodes(eidx(i),:);
    endpoints = findnode(G,{endpoints{1} endpoints{2}});
    if (G.Nodes.discN(endpoints(1)) == -Inf)
        S_test1 = outedges(G, endpoints(1));
    end
    if (G.Nodes.discN(endpoints(2)) == -Inf)
        S_test1 = outedges(G, endpoints(2));
    end

    S_test1 = min(S_test1);
    S_test = cat(2,S_test,S_test1');
    S_test = unique(S_test);
end

eidx = S_test;
    endpts = G.Edges.EndNodes(eidx,:);
    endpts = findnode(G,endpts);
   
    w_id = [];
    for i = 1:(size(endpts,1)*size(endpts,2))
        if (isinf(G.Nodes.discN(endpts(i))))
            newN = endpts(i);
            w_id(end+1) = [newN];
            
        else
            pre_id = endpts(i);
        end
        
    end
    w_id = unique(w_id);
    for i = 1:length(w_id)
        G.Nodes.discN(w_id(i)) = currentDf +i;
    end
    
    

    % add the new node to the tree
    newNode = table(G.Nodes.Name(w_id), w_id', G.Nodes.discN(w_id),'VariableNames', {'Name','origId', 'discN'});
    T = addnode(T,newNode);

    % create the edge and its attributes (endpts and original id in G) to be added in T
    newEdge = table(G.Edges.EndNodes(eidx,:),G.Edges.Name(eidx),eidx','VariableNames', {'EndNodes','Name','origId'});
    T = addedge(T,newEdge);

     S = bfsupdateFrontierEdge(G,S);
currentDf = max(G.Nodes.discN);
    
    

end

end % end function bfsSpanningTree