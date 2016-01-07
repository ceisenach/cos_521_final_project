function [ A ] = GraphGen( n, p )

    %generate an Erdos-Renyi Graph
    %inputs:
    % - n: size of the graph
    % - p: probability of the ER


    A = zeros(n); %adjacency matrix of the graph
    for i=1:n-1
        for j=i+1:n            
            A(i,j) = (rand < p);
        end
    end
    
    
    A = A + A';
    
    
    
end



