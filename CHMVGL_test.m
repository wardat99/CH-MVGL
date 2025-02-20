clear all
clc
close all
warning off;

%% This is a test for CH-MVGL method

addpath('functions')
max_num_run=10;       % The max rum of runs 
for num_run =1:1:max_num_run
n=2^8;              % Number of Nodes
p=0.1;               % Probability of connection for the graph
hn_perc=0.02;        % Percentage of co-hub nodes
num_views =6;        % Number of Views
num_coHub_nodes=int8(hn_perc*n); % Number of co-hub nodes 
[Adj, A,hn] = get_hub_graph(n,num_coHub_nodes,num_views,p); % generate co-hub graphs 
 for jj=1:num_views
     L(:,:,jj)=diag(sum(A(:,:,jj)))-A(:,:,jj); % compute the Laplacian matrix 
 end
n_signals=700;     % Number of signals
noise_amount=0.1;  % Noise level in each signal
for v=1:num_views
X{v} = gen_samples_new(A(:,:,v),n_signals,noise_amount,'heat')'; % generate the samples
end

[P] = generate_P(n)'; % generate matrix P
alpha=1;  
delta1=15; delta2=20; delta3=35;
gamma_1=delta1/p; gamma_2=delta2*n; gamma_3=delta3/hn_perc; gamma_4=delta1/p;

fprintf('CH-MVGL attempt number start: %.2f\n' ,num_run);
[Ak,Ck,G,time_cost(num_run)] = CHMVGL(X,P,gamma_1,gamma_2,gamma_3,gamma_4,alpha);
fprintf('CH-MVGL attempt number end: %.2f\n' ,num_run);

%% compute the F-score

for jj=1:num_views
    f(jj,num_run) = compute_f(Ak(:,:,jj),A(:,:,jj));
end
f_avg_each_run(num_run)=mean(f(:,num_run));
end
fprintf('The avg. F-score across views is: %.2f\n' , mean(f_avg_each_run));
fprintf('The avg. F-score across views is: %.2f\n' , mean(time_cost));
