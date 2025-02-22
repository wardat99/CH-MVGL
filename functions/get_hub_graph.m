function [Adj,A,hn] = get_hub_graph(n,num_nodes,num_views,p)


rng('shuffle');
hn=[];
x = setdiff(1:n,hn);
percentage_rng=0.88:0.2:1;
percentage=percentage_rng(randi(numel(percentage_rng)));
for k=1:num_views
    Adj(:,:,k) = create_ER_Graph(n,p);

end
A=Adj;
for i=1:num_nodes
    hn(i)=x(randi(numel(x)));
    %rv = int32(randi([0, 1], [1, n]));
    rv =get_vector(n,percentage);
    for  j=1:num_views 
        %rn(i)=x(randi(numel(x)));
        
        A(hn(i),:,j)=rv;
        A(:,hn(i),j)=rv';
    end
    x=setdiff(1:n,hn);
end


    function vector=get_vector(n,percentage)
    num_ones = round(percentage * n); % number of 1's 
    num_zeros = n - num_ones; % number of 0's 
    vector = [ones(1, num_ones), zeros(1, num_zeros)];
    vector = vector(randperm(n));
    end


end