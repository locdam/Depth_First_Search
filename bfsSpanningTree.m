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

%From line 43 to line , we run the loop one time manually before going into the while loop.
    
    % generate the next edges using bfsNextEdge, however, this next edge does not 
    % have to be only 1 edge, it could also a array of edges.
    eidx = bfsNextEdge(G,S);
    
    %determine the nodes of those next edges
    endpts = G.Edges.EndNodes(eidx,:);
    endpts = findnode(G,endpts);;
   
    %since the next edges is an array, we have to use the while loop through
    %each of the nodes to determine either they are discovered not not-discovered
    %if they are new nodes, then append them to w_id
    
    w_id = [];
    for i = 1:(size(endpts,1)*size(endpts,2))
        if (isinf(G.Nodes.discN(endpts(i))))
            newN = endpts(i);
            w_id(end+1) = [newN];
            
        else
            pre_id = endpts(i);
        end 
    end
    
    % run a loop through w_id to update their dfN
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
    
    % use bfsupdateFrontierEdge function to determine the next set of S.
    S = bfsupdateFrontierEdge(G,S);
    
    %update currentDf, since after each loop, we can have a new set of discovered nodes, we update currentDf as the max of the most 
    % recently discovered nodes.
    currentDf = max(G.Nodes.discN);
    
    % Now, we go through while loop as long as S is not empty
while ~isempty(S)
    
    % generate the next edges using bfsNextEdge, however, this next edge does not 
    % have to be only 1 edge, it could also a array of edges.
    eidx = bfsNextEdge(G,S);
    
    % This S_test will fillter out the edges where there only edges derive from
    % the nodes with dfN = -Inf
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
    % There could be multiple edges derive from the -Inf nodes, this S_test1 will 
    % filter out the only edge with smallest dfN
    S_test1 = min(S_test1);
    
    % Append all the edges and into 1 set and filter out the unique elements in them 
    S_test = cat(2,S_test,S_test1');
    S_test = unique(S_test);
end

    % Set the next edges eidx as S_test
    eidx = S_test;
    endpts = G.Edges.EndNodes(eidx,:);
    endpts = findnode(G,endpts);
   
    %Looping through the nodes of new edges to determine if they discovered or un-discovered
    w_id = [];
    for i = 1:(size(endpts,1)*size(endpts,2))
        if (isinf(G.Nodes.discN(endpts(i))))
            newN = endpts(i);
            w_id(end+1) = [newN];
            
        else
            pre_id = endpts(i);
        end
        
    end
    
    % update the dfN of discovered nodes
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


%%bfnextEdge
function [eidx] = bfsNextEdge(G,S)

%First, sort S from small to big
sort_S = sort(S);

%run a loop across the lenght of S

eidx = [];
for i=1:length(sort_S)
        
        % find the endpoints of S at each ith position
        endpoints = G.Edges.EndNodes(sort_S(i),:);
        endpoints = findnode(G,{endpoints{1} endpoints{2}});
        
        %determine the edges that is discovered or not by the dfN of their nodes.
        if (G.Nodes.discN(endpoints(1)) == -Inf) || (G.Nodes.discN(endpoints(2)) == -Inf)
            newE = sort_S(i);
            eidx(end+1) = [newE];
                    
        end
end
end

%%bfsupdateFrontierEdge
function S_new = bfsupdateFrontierEdge(G,S)
    %First, sort S from small to big
sort_S = sort(S);

S_new = [];
%Start a loop across the length of Nodes of G.
  for i = 1:length(sort_S)
 
      % Look for the node N, which is also the newest
      % discovered node
      endpoints = G.Edges.EndNodes(sort_S(i),:);
      endpoints = findnode(G,{endpoints{1} endpoints{2}});

      if (G.Nodes.discN(endpoints(1)) >0) 
      
        % Set S_new as the frontier edges from that node. However, this
        % S includes both discovred and undiscovered edges
        [S_new1,nV] = outedges(G,G.Nodes.Name(endpoints(1)));
      end
      if (G.Nodes.discN(endpoints(2)) >0) 
          [S_new1,nV] = outedges(G,G.Nodes.Name(endpoints(2)));
        %S_new1 = S_new1'
      end
      S_new = cat(2,S_new,S_new1');
      S_new = unique(S_new);
  end

%Now filter out the discovered edges
%set an empty discovered S for later
discoveredS = [];

%start a loop across the length of the S_new found above
for k = 1:length(S_new)
        
        %locate the endpoints of each edge at each k-th
        endpoints = G.Edges.EndNodes(S_new(k),:);
        endpoints = findnode(G,{endpoints{1} endpoints{2}});
        
        %if the dfN of both of the endpoints of that k-th edge is greater
        %or equal to 0, then that edges is discovred.
        if (G.Nodes.discN(endpoints(1))== -Inf) || (G.Nodes.discN(endpoints(2)) == -Inf)
            
            %call the k-th discovered edge is newS
            newS = S_new(k);
            
            %append all the discovered edges to discoveredS
            discoveredS(end+1) = [newS];
    
        end
        

end

%use setdiff to extract the undiscovered from S_new vs the discovredS.

%S_new = setdiff(S_new, discoveredS);
S_new = discoveredS;
end