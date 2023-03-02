function [eidx] = bfsNextEdge(G,S)

%First, sort S from small to big
sort_S = sort(S);

%run a loop across the lenght of S

eidx = [];
for i=1:length(sort_S)
        
        % find the endpoints of S at each ith position
        endpoints = G.Edges.EndNodes(sort_S(i),:);
        endpoints = findnode(G,{endpoints{1} endpoints{2}});

        if (G.Nodes.discN(endpoints(1)) == -Inf) || (G.Nodes.discN(endpoints(2)) == -Inf)
            newE = sort_S(i);
            eidx(end+1) = [newE];
                    
        end
end
% S_test = [];
% for i = 1:length(eidx)
%     endpoints = G.Edges.EndNodes(eidx(i),:);
%     endpoints = findnode(G,{endpoints{1} endpoints{2}});
%     if (G.Nodes.dfN(endpoints(1)) == -Inf)
%         S_test1 = outedges(G, endpoints(1));
%     end
%     if (G.Nodes.dfN(endpoints(2)) == -Inf)
%         S_test1 = outedges(G, endpoints(2));
%     end
% 
%     S_test1 = min(S_test1);
%     S_test = cat(2,S_test,S_test1');
%     S_test = unique(S_test);
% end
% 
% eidx = S_test;
end