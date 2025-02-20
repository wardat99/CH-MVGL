function [Ak,Ck,G,time_cost] = CHMVGL(Xk,P,gamma_1,gamma_2,gamma_3,gamma_4,alpha)



n=size(Xk{1},1);
[~,K]=size(Xk);
for k=1:K
    Bk(:,:,k)=P'*Xk{k}*Xk{k}'*P;
end

%% Initialization
Ck=zeros(n,n,K);
V=zeros; W=V; G=V;  Zk=Ck; 
Yk=zeros(n,n,K); Jk=Yk; Mk=Yk; N=zeros(n,n); Q=N; Ei=zeros(n-1,n-1,k); Tk=Yk;
I=eye(n); Ci=Ck;
max_iter = 500;
tol=1e-3;


tic;
for iter=1:max_iter
    
    %% Update Ek, Zk
    for k=1:K
        Ek(:,:,k)=(2*gamma_1*P'*Zk(:,:,k)*P+alpha*P'*Ck(:,:,k)*P+P'*Yk(:,:,k)*P-Bk(:,:,k)')./(2*gamma_1+alpha);
    end

    for k=1:K
        Uk(:,:,k)=2*gamma_1*P*Ek(:,:,k)*P'+alpha*I.*Ck(:,:,k)-Jk(:,:,k);
        Zk(:,:,k)=(Uk(:,:,k)+(Uk(:,:,k).^2+4*(2*gamma_1+alpha)*gamma_2*I).^0.5)./(4*gamma_1+2*alpha);
        Zk(:,:,k)=I.*Zk(:,:,k);
    end

    %% Update Ei, Ck, Ci
    for k=1:K
        Ei(:,:,k)=(alpha*P'*Ci(:,:,k)*P+P'*Tk(:,:,k)*P)./(2*gamma_4+alpha);
    end

   

    for k=1:K
        Ck_d=I.*((alpha*P*Ek(:,:,k)*P'-Yk(:,:,k)+alpha*Zk(:,:,k)+Jk(:,:,k)+alpha*Ci(:,:,k)+alpha*V+alpha*W-Mk(:,:,k))./(3*alpha));
        Ck_d(Ck_d<0)=0;
        Ck_off=(alpha*P*Ek(:,:,k)*P'-Yk(:,:,k)+alpha*Ci(:,:,k)+alpha*V+alpha*W-Mk(:,:,k))./(2*alpha);
        Ck_off(Ck_off>0)=0;
        Ck(:,:,k)=Ck_d+Ck_off; Ck(:,:,k)=0.5*(Ck(:,:,k)+Ck(:,:,k)');
    end
    for k=1:K
        Cii = (alpha*Ck(:,:,k)-alpha*V-alpha*W+Mk(:,:,k)+alpha*P*Ei(:,:,k)*P'-Tk(:,:,k))./(2*alpha);
        Ci_d=Cii.*I; Ci_d(Ci_d<0)=0; 
        Ci_off=Cii-Cii.*I; Ci_off(Ci_off>0)=0; 
        Ci(:,:,k)=Ci_d+Ci_off; Ci(:,:,k)=0.5*(Ci(:,:,k)+Ci(:,:,k)'); clear Cii Ci_d Ci_off
    end

    %% Update W, V, G

    sumW=0;
    for k=1:K
       sumW=sumW+alpha*Ck(:,:,k)-alpha*Ci(:,:,k)-alpha*V+Mk(:,:,k);
    end
    W=(sumW+alpha*V'+N')./(alpha*(K+1));

    sumV=0;
    for k=1:K
        sumV=sumV+alpha*Ck(:,:,k)-alpha*Ci(:,:,k)-alpha*W+Mk(:,:,k);
    end
    V=(sumV+alpha*W'-N+alpha*G+Q)./(alpha*(K+2));
    
    
    G = solve_l1l2(V-Q./alpha,gamma_3./(2*alpha));
    % GG=V-Q./alpha; GGd=diag(GG);
    % 
    %     G_nd=nodiag_construction(GG);
    % 
    % 
    %     H = soft_scal(G_nd,gamma_2/(2*alpha)); 
    %     G=diag_construction(H,GGd);
        



    %% Update Lag.
    for k=1:K
        dYk(:,:,k)=Ck(:,:,k)-P*Ek(:,:,k)*P';
        Yk(:,:,k)=Yk(:,:,k)+alpha*dYk(:,:,k);

        dJk(:,:,k)=Zk(:,:,k)-I.*Ck(:,:,k);
        Jk(:,:,k)=Jk(:,:,k)+alpha*dJk(:,:,k);

        dMk(:,:,k)=Ck(:,:,k)-Ci(:,:,k)-V-W;
        Mk(:,:,k)=Mk(:,:,k)+alpha*dMk(:,:,k);

        dTk(:,:,k)=Ci(:,:,k)-P*Ei(:,:,k)*P';
        Tk(:,:,k)=Tk(:,:,k)+alpha*dTk(:,:,k);
    end
    
    dN=V-W'; N=N+alpha*dN;
    dQ=G-V;  Q=Q+alpha*dQ;

    alpha=1.1*alpha;



    %% stopping criteria
    
    err(1)=norm(dYk(:),'fro'); err(2)=norm(dJk(:),'fro');
    err(3)=norm(dMk(:),'fro'); err(4)=norm(dN(:),'fro');
    err(5)=norm(dQ(:),'fro'); err(6)=norm(dTk,'fro');

    %if mod(iter,10) == 0
    %    disp(['iter ' num2str(iter) ', err=' num2str(err)])
    %end
    err_max=max(err);

    if err_max < tol
       break;
    end
end
time_cost=toc;


ThValue=0.01;
        for k=1:K
        Akk=diag(diag(Ck(:,:,k)))-Ck(:,:,k); Akk=Akk./norm(Akk);
        Akk(Akk>ThValue)=1; 
        Akk(Akk<ThValue)=0;
        Ak(:,:,k)=Akk; clear Akk
                        
        Aii=diag(diag(Ci(:,:,k)))-Ci(:,:,k); Aii=Aii/norm(((Aii)));
        Aii(Aii>ThValue)=1; 
        Aii(Aii<ThValue)=0; 
        Ai(:,:,k)=Aii; clear Aii
        end


end