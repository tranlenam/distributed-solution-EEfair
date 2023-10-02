function [eta, mu,u,W,obj,res]=MaxMin_ADMM_mosek_global_two_c(nUsers,nBS,zeta_n,xi_n,Tau_n,eta_n,c1,c2)


varlength = nUsers*nBS*nBS+1 + (2*(nBS-1)*nUsers+1)*nBS+3;

%%  x = mu-tau => tau = mu - x  x1=tau-mu => x+mu=tau
a1=[];
for b=1:nBS
for iBS=1:nBS
    if b~=iBS
        for i=1:nUsers
            temp=zeros(1,varlength);
            temp(1,[(iBS-1)*nUsers*nBS+(b-1)+(i-1)*nUsers+i])=1;
            a1=[a1;temp];     
        end
    
    end
end
for i=1:nUsers    
        for iBS=1:nBS
            if b~=iBS
            temp=zeros(1,varlength);
            temp(1,[(b-1)*nUsers*nBS+(iBS-1)+(i-1)*nUsers+i])=1;
            a1=[a1;temp];     
        end
    
    end
end
end   
nlength1 = nUsers*nBS*nBS;
for i=1:(size(a1,1))
    a1(i,nlength1+i)=sqrt(2/c2);
end
nlength2 = nlength1+size(a1,1);
%% eta_n - eta = x2 => eta+x2=eta_n;
a2=[];
for i=1:nBS
    temp = zeros(1,varlength);
    temp(1,[nlength2+1,nlength2+1+i])=[1 ,sqrt(2/c1)];
    a2=[a2;temp];
end
nlength3 = nlength2+1+nBS;
%% obj >= -sum(eta_n)+zeta1*x1-zeta2*x2...+ c/2(sum[()^2]
% obj+sum(eta_n) - sqrt(2/c2)*zeta1*x1-sqrt(2/c1)*zeta2*x2... >= (sum[()^2]
%( obj+sum(eta_n) - sqrt(2/c2)*zeta1*x1-zeta2*x2... )>= 2X
zeta = myvec(zeta_n);
xi = myvec(xi_n);
a3=zeros(1,varlength);
a3(1,[nlength3+1,nlength1+1:nlength2,nlength2+2:nlength3,nlength3+2])=[1 -sqrt(2/c2)*zeta', -sqrt(2/c1)*xi',-2];
nlength4=nlength3+2; % for obj and X
prob.cones{1}.type = 'MSK_CT_RQUAD';
prob.cones{1}.sub = [nlength4,nlength4+1,nlength1+1:nlength2,nlength3-nBS+1:nlength3];
            


%% 
prob.c=[zeros(1,nlength3) 1 zeros(1,2) ];
prob.a=sparse([a1;a2;a3]);

prob.blc = [myvec(Tau_n)',... %a1
            myvec(eta_n)',... %a2
           -sum(eta_n),... %a3
           ];   

prob.buc = [myvec(Tau_n)',... %a1
            myvec(eta_n)',... %a2
           -sum(eta_n),... %a3
          ]; 
      
           


prob.blx = [ zeros(1,nUsers*nBS*nBS),...%mu
             -inf*ones(1,(2*(nBS-1)*nUsers)*nBS),...%x1
             0,... %eta
             -inf*ones(1,nBS),...%x2
             -inf,...%obj
             0 1,...%2X and 1
             ];
             
prob.bux = [ inf*ones(1,nUsers*nBS*nBS),...%mu
             inf*ones(1,(2*(nBS-1)*nUsers)*nBS),...%x1
             inf,... %eta
             inf*ones(1,nBS),...%x2
             inf,...%obj
             inf 1,...%2X and 1
             ];
             


% tic
[r,res]       = mosekopt('minimize echo(0) info',prob);     
% toc
W = res.sol.itr.xx;
% time_solving = [time_solving res.info.MSK_DINF_OPTIMIZER_TIME];
% no_flops = no_flops + res.info.MSK_DINF_INTPNT_FACTOR_NUM_FLOPS;
% res.sol.itr.solsta
 mu=reshape(W(1:nUsers*nBS*nBS),nBS,nUsers,nBS);
 eta=W(nlength2+1);
 for b=1:nBS
 u(:,b) = [ myvec(mu(b,:,1:nBS~=b))' myvec(mu(1:nBS~=b,:,b))']';
 end
%  obj1=W(nlength3+1);
 

obj = sum(-eta_n)+sum(myvec(zeta_n.*(Tau_n-u)))+sum(xi_n.*(eta_n-eta))+c1/2*norm(eta_n-eta)^2+c2/2*norm(myvec(Tau_n-u))^2;


    

% if all((alpha-(jj-1)*maxrange)<=0) || all((jj*maxrange-alpha)<=0)
    


