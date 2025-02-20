function [S] = l21_operator(L,lambda1,mu)


[n,~,M]=size(L);

F=double(tenmat(L,2,'t'));      
[Sconcat] = solve_l1l2(F,lambda1/mu);
for j=1:M
     S(:,:,j)=Sconcat((j-1)*n+1:j*n,:);
end



end