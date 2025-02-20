function [X,e] = gen_samples_new(A,n_signals,noise_amount,filter)
%
%   

[n,~]=size(A);
mu=zeros(1,n);
h=zeros(n,n);
alpha=6;

if filter == "gaussian"
% Get the graph Laplacian spectrum
L =diag(sum(A)) - A;
[V,e] = eig(L);

% Normalize the L matrix

 e(e< 1e-8)=0;
 e=e./max(max(e));

% generate samples
Sigma=pinv(L);
X0 = mvnrnd(mu,Sigma,n_signals);
X0_hat = V'*X0'; 


    for i=1:n
        for j=1:n
            if e(i,j)>0
                h(i,j)=1./sqrt(e(i,j));
            else
                h(i,j)=0;
            end
        end
    end
    X=V*diag(diag(h))*X0_hat;
elseif  filter == "tikhonov"
    L =diag(sum(A)) - A;
[V,e] = eig(L);

% Normalize the L matrix

 e(e< 1e-8)=0;
 e=e./max(max(e));

% generate samples
Sigma=L;
X0 = mvnrnd(mu,Sigma,n_signals);
X0_hat = V'*X0'; 

     for i=1:n
        for j=1:n
            h(i,j)=1./(1+alpha*e(i,j));
        end
     end 
     X=V*diag(diag(h))*X0_hat;
elseif filter == "heat"
    L =diag(sum(A)) - A;
[V,e] = eig(L);

% Normalize the L matrix

 e(e< 1e-8)=0;
 e=e./max(max(e));

% generate samples
Sigma=L;
X0 = mvnrnd(mu,Sigma,n_signals);
X0_hat = V'*X0'; 

    for i=1:n
        for j=1:n
            h(i,j)=exp(-alpha*e(i,j));
        end
    end 
     X=V*diag(diag(h))*X0_hat;
end



% if filter == "Gaussian"
%     for i=1:n
%         for j=1:n
%             if e(i,j)>0
%                 h(i,j)=1./sqrt(e(i,j));
%             else
%                         h(i,j)=0;
%             end
%         end
%     end
% 
% elseif  filter == "tikhonov"
%     h=1./(1+alpha*e);
% elseif filter == "heat"
%     h=exp(-alpha*e);
% end

% add noise
normX=norm(X);
E = normrnd(0,1,size(X));
normE=norm(E);
X=(X+E*(noise_amount*normX/normE))';
end