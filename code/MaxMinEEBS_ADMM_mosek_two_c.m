function [eta_b, beamformer_b, z_b, t_b, q_b, tau, tau_pr,Tau,W,obj,time_solving,noflop,res]=MaxMinEEBS_ADMM_mosek_two_c(channel, power,Pdyn,Psta, scalefactor,beamformer_i,q_i,z_i,t_i,eta, u, b, nBS, nUsers,nTx,PAeff, xi_n, zeta_n, c1,c2)
param = []; 
 
% Primal feasibility tolerance for the primal solution 
param.MSK_DPAR_INTPNT_CO_TOL_PFEAS = 1.0e-9; 
 
% Dual feasibility tolerance for the dual solution 
param.MSK_DPAR_INTPNT_CO_TOL_DFEAS = 1.0e-9; 
 
% Relative primal-dual gap tolerance. 
param.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1.0e-9; 
% b=1;
m=10;
% u = u(:,b);
varlength = 2*nUsers*nTx+1 + 2*nUsers*(nBS-1)*nUsers+nUsers*(nBS-1) + 2*nUsers*(nBS-1) + nUsers*(nBS-1)+nUsers+ nUsers +nUsers+nUsers...
            + 2*nUsers*(nUsers-1)+nUsers + 3 + 2*nUsers*nTx + nUsers+1+1+1 + 2*nUsers + 3*nUsers+3*(m+3)*nUsers + 1 + 2*nUsers*(nBS-1)+1 + 3;

nlength3=2*nUsers*nTx+1;% 1 due to adding nBS variables for sqrt(power) in conic constraints

%% norm(w)< power

    prob.cones{1}.type = 'MSK_CT_QUAD';
    prob.cones{1}.sub = [2*nUsers*nTx+1,1:2*nTx*nUsers];

% cindex3 = size(a1,1)+size(a2,1)+size(a3,1);
%% q_k>= |SINR_kb|^2


% tau = nUser*nBS
%%
% tau > sum(| |^2)
a4=[]; % x=Re(wh) x=Im(wh) 
for iBS=1:nBS
   if iBS~=b
      for iUser=1:nUsers
          for i=1:nUsers
                    temp = zeros(2,varlength);
                    temp(1,(i-1)*2*nTx+1:(i)*2*nTx)...
                          =[ myvec([ real(channel(:,iUser,iBS,b))';  imag(channel(:,iUser,iBS,b))'])'];
                    temp(2,(i-1)*2*nTx+1:(i)*2*nTx )...
                          =[ myvec([-imag(channel(:,iUser,iBS,b))';  real(channel(:,iUser,iBS,b))'])'];
                    a4= [a4;temp]; 
          end                    
      end
   end
end
for i=1:size(a4,1)    
    a4(i,nlength3+i)=-1;
end
nlength4=nlength3+2*nUsers*(nBS-1)*nUsers+nUsers*(nBS-1); %(nBS-1)*nUsers for vars sigma=1;

% 2*tau*1/2 >= Re^2 + im^2
for i=1:nUsers*(nBS-1)    
    prob.cones{1+i}.type = 'MSK_CT_RQUAD';
    prob.cones{1+i}.sub =  [nlength4+i,nlength4+nUsers*(nBS-1)+i,nlength3+2*nUsers*(i-1)+1:nlength3+2*nUsers*i];
end
nlength4=nlength4+2*nUsers*(nBS-1);%nUsers*(nBS-1) vars for tau and nUsers*nBS x=1/2;
nlength4temp=nlength4;
%%
% (q_bk-sum(tau_pr_bk)+1)/2 >= 
%   (q_b(iUser)-sum(myvec(tau_pr(1:nBS~=b,iUser)))+1)/2...
%                       >= norm([1 channel(:,iUser,b,b)'*beamformer_b(:,1:nUsers~=iUser) (q_b(iUser)-sum(myvec(tau_pr(1:nBS~=b,iUser)))-1)/2]);  
% sum(myvec(tau_pr(1:nBS~=b,iUser))) <= x
a41=[]; % =
for i=1:nUsers
    temp=zeros(1,varlength);
    temp(1,[nlength4+(i-1)*(nBS-1)+1:nlength4+(i)*(nBS-1),nlength4+nUsers*(nBS-1)+i])=[ones(1,nBS-1),-1];
    a41=[a41;temp];
end
nlength41=nlength4+nUsers*(nBS-1)+nUsers+ nUsers +nUsers+nUsers;  %tau_pr and x and q_bk and x' and 1
%q_bk - x >= norm(1 | |^2)
%q_bk-x>=2x';
a42=[];
for i=1:nUsers
    temp=zeros(1,varlength);
    temp([nlength41-3*nUsers+i,nlength41-4*nUsers+i,nlength41-2*nUsers+i])=[1 -1 -2];
    a42=[a42;temp];
end
% Re = k1  Im=k2
a43=[]; 
      for iUser=1:nUsers
          for i=1:nUsers
              if i~=iUser
                    temp = zeros(2,varlength);
                    temp(1,(i-1)*2*nTx+1:(i)*2*nTx)...
                          =[ myvec([ real(channel(:,iUser,b,b))';  imag(channel(:,iUser,b,b))'])']; % Re
%                     temp(1,nlength3++1)=-1;
                    temp(2,(i-1)*2*nTx+1:(i)*2*nTx )...
                          =[ myvec([-imag(channel(:,iUser,b,b))';  real(channel(:,iUser,b,b))'])']; % Im
%                     temp(2,nlength3+nUsers*(iBS-1)+2)=-1;
                    a43= [a43;temp]; 
              end 
          end
      end
for i=1:size(a43,1)
    a43(i,nlength41+i)=-1;
end
nlength43=nlength41+2*nUsers*(nUsers-1)+nUsers; %2*nUsers*(nUsers-1) for Re Im and nUsers for vars sigma;

% 2x'*1 >= norm(1+|Re|^2+|Im|^2)
for i=1:nUsers    
    prob.cones{1+nUsers*(nBS-1)+i}.type = 'MSK_CT_RQUAD';
    prob.cones{1+nUsers*(nBS-1)+i}.sub =  [nlength41-2*nUsers+i,nlength41-nUsers+i,nlength41+2*(nUsers-1)*(i-1)+1:nlength41+2*(nUsers-1)*(i),nlength43-nUsers+i];
end
nlength4=nlength43;
%% PAeff*(t_b-P0)*scalefactor^2=2x => PAeff*scalefactor^2*tb-2x=P0*PAeff*scalefactor^2

a5=[];
%     for iUser=1:nUsers
        temp=zeros(1,varlength);
        temp([nlength4+1,nlength4+1+1])=[PAeff*scalefactor^2 -2];
        a5=[a5;temp];

nlength5=nlength4+3;% t_b and 2x and 1;
%% 2x*1>=sum(w)
% x6=w
a6=[];
    for ii=1:2*nUsers*nTx
        temp = zeros(1,varlength);
        temp([ii,nlength5+ii])=[1 -1];
        a6=[a6;temp];
    end

% 2x*1>=sum(w)

    prob.cones{1+nUsers*(nBS-1)+nUsers+1}.type = 'MSK_CT_RQUAD';
    prob.cones{1+nUsers*(nBS-1)+nUsers+1}.sub = [nlength5-1,nlength5,nlength5+1:nlength5+2*nUsers*nTx];

nlength6=nlength5+2*nUsers*nTx;%w
% cindex6=size(a1,1)+size(a2,1)+size(a3,1)+size(a4,1)+size(a5,1)+size(a6,1);
%% sum(log(1+g_kb) >= z_b^2 => sum(x7)>=z^2; 1+g_bk >= exp(x_7k)
% sum(x7)=x71+...+x7n=2x7'
a7=zeros(1,varlength);

    a7(1,[nlength6+1:nlength6+nUsers nlength6+nUsers+1])=[ones(1,nUsers) -2];

nlength7=nlength6+nUsers+1+1+1;%nUsers for x7, 1 for 2x7', 1 for 1, 1 for z_b
% cindex7=cindex6+size(a7,1);
%2x7'*1>=z^2

    prob.cones{1+nUsers*(nBS-1)+nUsers+1+1}.type = 'MSK_CT_RQUAD';
    prob.cones{1+nUsers*(nBS-1)+nUsers+1+1}.sub = [nlength6+nUsers+1,nlength6+nUsers+2,nlength6+nUsers+3];


%%
% 1+g_kb>=exp(x7);
% 1+g_kb>=x8 => g_kb-x8>=-1
a8=zeros(1,varlength);
for i=1:nUsers
    a8(i,[nlength7+i nlength7+nUsers+i])=[1 -1];
end
nlength8=nlength7+2*nUsers;%g_b x8
%% x8 >=exp(x7)
% m=10;


conelength=1+nUsers*(nBS-1)+nUsers+1+1;
a9=zeros(1,varlength);
a10=[];
a11=[];
a12=zeros((m+3)*nUsers,varlength);
a12(:,[nlength8+3*nUsers+1:nlength8+3*nUsers+nUsers*(m+3), nlength8+3*nUsers+2*nUsers*(m+3)+1:nlength8+3*nUsers+3*nUsers*(m+3)])= [eye(nUsers*(m+3)) -eye(nUsers*(m+3))];
for i=1:nUsers
    a9(i,[nlength7+nUsers+i nlength8+i])=[1 -1]; %x8-t>=0; t=nlength8+1:nlength8+nUsers*nBS
    for ii=1:m+3
        if ii<5
            if ii==1
                temp=zeros(1,varlength);
                temp([nlength6+i nlength8+nUsers+i])=[1/2^m -1];%x7/2^m-x9_1=-1; x9_1=nlength8+nUsers*nBS+1:nlength8+2*nUsers*nBS
                prob.cones{conelength+(i-1)*(m+3)+ii}.type = 'MSK_CT_RQUAD';
                prob.cones{conelength+(i-1)*(m+3)+ii}.sub = [nlength8+3*nUsers+(i-1)*(m+3)+ii,nlength8+3*nUsers+nUsers*(m+3)+(i-1)*(m+3)+ii,nlength8+nUsers+i];% u1_kb * (2x=1) >=x9_kb^2
                a10=[a10;temp];
            end
            if ii==2
                temp=zeros(1,varlength);
                temp([nlength6+i nlength8+2*nUsers+i])=[1/2^(m+1) -1];%x7/2^(m+1)x9_2=-5/6; x9_2=length8+2*nUsers*nBS+1:length8+3*nUsers*nBS
                prob.cones{conelength+(i-1)*(m+3)+ii}.type = 'MSK_CT_RQUAD';
                prob.cones{conelength+(i-1)*(m+3)+ii}.sub = [nlength8+3*nUsers+(i-1)*(m+3)+ii,nlength8+3*nUsers+nUsers*(m+3)+(i-1)*(m+3)+ii,nlength8+2*nUsers+i];% u2_kb * (2x=1) >=x9_kb^2
                a10=[a10; temp];
            end
            if ii==3
%                 a10((i-1)*(m+3)+ii,[nlength6+i nlength8+nUsers*nBS+(i-1)*nUsers*nBS+ii])=[0 0];% nothing
                prob.cones{conelength+(i-1)*(m+3)+ii}.type = 'MSK_CT_RQUAD';
                prob.cones{conelength+(i-1)*(m+3)+ii}.sub = [nlength8+3*nUsers+(i-1)*(m+3)+ii,nlength8+3*nUsers+nUsers*(m+3)+(i-1)*(m+3)+ii,nlength8+3*nUsers+2*nUsers*(m+3)+(i-1)*(m+3)+1];% u3_kb * (2x=1) >=u1_kb^2
            end
            if ii==4
                temp=zeros(1,varlength);
                temp([nlength8+3*nUsers+(i-1)*(m+3)+2, nlength8+3*nUsers+(i-1)*(m+3)+3, nlength8+3*nUsers+(i-1)*(m+3)+4])=[1 1/24 -1];% u2_kb+1/24*u3_kb-u4 <= -19/72
                a11=[a11;temp];
            end
            % note:
            % u= nlength8+3*nUsers*nBS+1:nlength8+3*nUsers*nBS+(m+3)*nUsers*nBS
            % 2x=nlength8+3*nUsers*nBS+(m+3)*nUsers*nBS+1:nlength8+3*nUsers*nBS+2*(m+3)*nUsers*nBS
            % u'=u=nlength8+3*nUsers*nBS+2*(m+3)*nUsers*nBS+1:nlength8+3*nUsers*nBS+3*(m+3)*nUsers*nBS
        else 
            %a10((i-1)*(m+3)+ii,:)=0;% nothing
            prob.cones{conelength+(i-1)*(m+3)+ii-1}.type = 'MSK_CT_RQUAD';
            prob.cones{conelength+(i-1)*(m+3)+ii-1}.sub = [nlength8+3*nUsers+(i-1)*(m+3)+ii,nlength8+3*nUsers+nUsers*(m+3)+(i-1)*(m+3)+ii-1,nlength8+3*nUsers+2*nUsers*(m+3)+(i-1)*(m+3)+ii-1];% u_i_kb * (2x=1) >=u_i-1_kb^2        
        end
    end
    prob.cones{conelength+(i-1)*(m+3)+ii}.type = 'MSK_CT_RQUAD';
    prob.cones{conelength+(i-1)*(m+3)+ii}.sub = [nlength8+i,nlength8+3*nUsers+nUsers*(m+3)+(i-1)*(m+3)+ii,nlength8+3*nUsers+2*nUsers*(m+3)+(i-1)*(m+3)+ii];% t* (2x=1) >=u_m+3_kb^2
end
nlength10=nlength8+3*nUsers+3*(m+3)*nUsers; %t,x9_1, x9_2 , 2* u, 1
conelength = conelength + (m+3)*nUsers;
bound10=myvec([-ones(1,nUsers);-5/6*ones(1,nUsers)])';

%% phi_b = ()z_b -()t_b - eta_b>= 0
a13=zeros(1,varlength);

    a13(1,[nlength4+1,nlength7,nlength10+1])=[-(z_i(b)/t_i(b))^2, 2*z_i(b)/t_i(b),-1];

nlength12=nlength10+1; %1 for eta_b
% cindex11=cindex9+size(a10,1);
%% psi=2real()-()*q_kb>=g_kb
% 
a14=[];

    for iUser=1:nUsers
        whh = beamformer_i(:,iUser,b)'*channel(:,iUser,b,b)*channel(:,iUser,b,b)';
        temp=zeros(1,varlength);
        temp((iUser-1)*2*nTx+1:(iUser)*2*nTx ) =2/q_i(iUser,b)*[ myvec([ real(whh); -imag(whh)])'];  %Re
        temp(nlength41-3*nUsers+iUser)= -(abs(channel(:,iUser,b,b)'*beamformer_i(:,iUser,b))/q_i(iUser,b))^2;
        temp(nlength7+iUser)=-1;
        a14=[a14;temp];
    end

nlength14=nlength12;
%% objective : obj >= -eta_b + zeta'*(theta_b-mu_b)+c/2*||theta_b-mu_b||^2
% theta_b=[tau tau_pr zeta_b]
% mu_b = [mu zeta];
% 2/c(eta_b-zeta*(theta_b-mu_b)+obj) >= sum[(theta_b(i)-mu_b(i))^2]
%  theta_b  - mu_b =  sqrt(2/c2)*x1
a_obj_theta =[];
for i=1:nUsers*(nBS-1)
    temp=zeros(1,varlength);
    temp(1,[nlength4temp-2*nUsers*(nBS-1)+i,nlength12+i])=[1 -sqrt(2/c2)];
    a_obj_theta=[a_obj_theta;temp];
end
for i=1:nUsers*(nBS-1)
    temp=zeros(1,varlength);
    temp(1,[nlength41-nUsers*(nBS-1)-4*nUsers+i,nlength12+nUsers*(nBS-1)+i])=[1 -sqrt(2/c2)];
    a_obj_theta=[a_obj_theta;temp];
end
% eta_n - sqrt(2/c1)*x1' = eta
    temp=zeros(1,varlength);
    temp(1,[nlength12,nlength12+2*nUsers*(nBS-1)+1])=[1 -sqrt(2/c1)];
    a_obj_theta=[a_obj_theta;temp];

nlength_obj_theta = nlength12+2*nUsers*(nBS-1)+1; % 2*nUsers*(nBS-1)+1 for x1
% 2/c(eta_b-sum[zeta*(theta_b-mu_b)]+obj) >= x1 =>
% (eta_b-sum[zeta*sqrt(2/c2)*x1]-sum[xi*sqrt(2/c1)*x1']+obj) >= 2x2
a_obj_linear = zeros(1,varlength);
a_obj_linear(1,[nlength12, nlength12+1:nlength12+2*nUsers*(nBS-1)+1,nlength_obj_theta+1,nlength_obj_theta+2]) ...
    = [1, -sqrt(2/c2)*zeta_n(:,b)',-sqrt(2/c1)*xi_n(b), 1 ,-2];
nlength_obj_linear=nlength_obj_theta+3; %1 for obj, 1 for 2 x2,  1 for 1
%2x2*1 >= sum (x1^2)

   prob.cones{conelength+1}.type = 'MSK_CT_RQUAD';
   prob.cones{conelength+1}.sub = [nlength_obj_theta+2,nlength_obj_theta+3,nlength12+1:nlength12+2*nUsers*(nBS-1)+1];


%% 
prob.c=[zeros(1,nlength_obj_theta) 1 zeros(1,2) ];
prob.a=sparse([a4;a41;a42;a43;a5;a6;a7;a8;a9;a10;a11;a12;a13;a14;a_obj_theta;a_obj_linear]);

prob.blc = [zeros(1,size(a4,1)),... %a4
           0*ones(1,size(a41,1)),... %a41
           zeros(1,size(a42,1)),... %a42
           zeros(1,size(a43,1)),... %a43
           (Psta + Pdyn*nTx)*PAeff*scalefactor^2*ones(1,size(a5,1)),... %a5
           zeros(1,size(a6,1)),... %a6
           zeros(1,size(a7,1)),... %a7
           -1*ones(1,size(a8,1)),... %a8
           zeros(1,size(a9,1)),... %a9
           bound10,... %a10
           -inf*ones(1,size(a11,1)),... %a11
           zeros(1,size(a12,1)),... %a12
           zeros(1,size(a13,1)),... %a13
           zeros(1,size(a14,1)),... %a14
           u',eta,... %atheta
           zeros(1,size(a_obj_linear,1))]; %aobj
    

prob.buc = [zeros(1,size(a4,1)),... %a4
           zeros(1,size(a41,1)),... %a41
           inf*ones(1,size(a42,1)),... %a42
           zeros(1,size(a43,1)),... %a43
           (Psta + Pdyn*nTx)*PAeff*scalefactor^2*ones(1,size(a5,1)),... %a5
           zeros(1,size(a6,1)),... %a6
           inf*ones(1,size(a7,1)),... %a7
           inf*ones(1,size(a8,1)),... %a8
           inf*ones(1,size(a9,1)),... %a9
           bound10,... %a10
           -19/72*ones(1,size(a11,1)),... %a11
           zeros(1,size(a12,1)),... %a12
           inf*ones(1,size(a13,1)),... %a13
           inf*ones(1,size(a14,1)),... %a14
           u',eta,... %atheta
           inf*ones(1,size(a_obj_linear,1))]; %aobj
           


prob.blx = [ -inf*ones(1,2*nTx*nUsers),... % w 
             sqrt(power),... %sqrt(power)
             -inf*ones(1,2*nUsers*(nBS-1)*nUsers),... % intercell-interference of tau
             ones(1,nUsers*(nBS-1)),... % sigma
             zeros(1,nUsers*(nBS-1)),... % tau
             0.5*ones(1,nUsers*(nBS-1)),... % x=1/2
             zeros(1,nUsers*(nBS-1)),... % tau_pr
             zeros(1,nUsers),... % x=sum(tau_pr)
             zeros(1,nUsers),... %q_bk
             zeros(1,nUsers),... %2x'=q_bk-x
             ones(1,nUsers),... % 1
             -inf*ones(1,2*nUsers*(nUsers-1)),... %Re Im itracell-interference of tau_pr
             ones(1,nUsers),... %sigma
             0 0 1,... % t_b, 2x, 1
             -inf*ones(1,2*nTx*nUsers),... % w
             -inf*ones(1,nUsers) -inf,... % x7, x7' 
             1 0,... % 1, z_b
             -ones(1,nUsers) zeros(1,nUsers),... % g_kb, x8
             -inf*ones(1,nUsers),... % kappa_0
             -inf*ones(1,2*nUsers),... % x9_1, x9_2
             -inf*ones(1,(m+3)*nUsers),... % kappa_i
             0.5*ones(1,(m+3)*nUsers),... % 1/2
             -inf*ones(1,(m+3)*nUsers),... % kappa_i'
             0,... %eta_b
             -inf*ones(1,(2*nUsers*(nBS-1)+1)),... % x1= theta-mu
             -inf,... %obj
             0 1];% x2 and 1
             
prob.bux = [ inf*ones(1,2*nTx*nUsers), ... % w 
             sqrt(power),... %sqrt(power)
             inf*ones(1,2*nUsers*(nBS-1)*nUsers),... % interference of tau
             ones(1,nUsers*(nBS-1)),... % sigma
             inf*ones(1,nUsers*(nBS-1)),... % tau
             0.5*ones(1,nUsers*(nBS-1)),... % x=1/2
             inf*ones(1,nUsers*(nBS-1)),... % tau_pr
             inf*ones(1,nUsers),... % x=sum(tau_pr)
             inf*ones(1,nUsers),... %q_bk
             inf*ones(1,nUsers),... %2x'=q_bk-x
             ones(1,nUsers),... % 1
             inf*ones(1,2*nUsers*(nUsers-1)),... %Re Im itracell-interference of tau_pr
             ones(1,nUsers),... %sigma
             inf inf 1,... % t_b, 2x, 1
             inf*ones(1,2*nTx*nUsers), ... % w 
             inf*ones(1,nUsers) inf,... % x7, x7'
             1 inf,... % 1, z_b
             inf*ones(1,nUsers) inf*ones(1,nUsers),... % g_kb, x8
             inf*ones(1,nUsers),... % kappa_0
             inf*ones(1,2*nUsers),... % x9_1, x9_2
             inf*ones(1,(m+3)*nUsers),... % kappa_i
             0.5*ones(1,(m+3)*nUsers),... % 1/2
             inf*ones(1,(m+3)*nUsers),... % kappa_i'
             inf,... %eta_b
             inf*ones(1,(2*nUsers*(nBS-1)+1)),... % x= theta-mu
             inf,... %obj
             inf 1]; %x2 and 1
             


% tic
[r,res]       = mosekopt('minimize echo(0) info',prob,param);     
% toc
W = res.sol.itr.xx;
% res.sol.itr.solsta
time_solving = [res.info.MSK_DINF_OPTIMIZER_TIME];
noflop=res.info.MSK_DINF_INTPNT_FACTOR_NUM_FLOPS;
% time_solving = [time_solving res.info.MSK_DINF_OPTIMIZER_TIME];
% no_flops = no_flops + res.info.MSK_DINF_INTPNT_FACTOR_NUM_FLOPS;

    for iUser=1:nUsers        
            reBeamformer = W((iUser-1)*nTx*2+1:2:iUser*nTx*2,1);
            imBeamformer = W((iUser-1)*nTx*2+2:2:iUser*nTx*2,1);
            beamformer_b(:,iUser) = reBeamformer+j*imBeamformer;        
    end
eta_b = W(nlength12,1);
tau_pr=zeros(nBS,nUsers);
tau_pr(1:nBS~=b,:) = reshape(W(nlength4temp+1:nlength4temp+nUsers*(nBS-1)),nBS-1,nUsers);
tau = zeros(nUsers,nBS);
tau(:,1:nBS~=b) = reshape(W(nlength4temp-2*nUsers*(nBS-1)+1:nlength4temp-nUsers*(nBS-1)),nUsers,nBS-1);
z_b = W(nlength7);
t_b = W(nlength4+1);
q_b = W(nlength41-3*nUsers+1:nlength41-2*nUsers);
 obj1 = W(nlength_obj_theta+1);
Tau= [W(nlength4temp-2*nUsers*(nBS-1)+1:nlength4temp-nUsers*(nBS-1)); W(nlength4temp+1:nlength4temp+nUsers*(nBS-1))];
 
obj = -eta_b+zeta_n(:,b)'*(Tau-u)+xi_n(b)*(eta_b-eta)+c1/2*(eta_b-eta)^2+c2/2*norm(Tau-u)^2;

% a=-(z_i(b)/t_i(b))^2*t_b+ 2*z_i(b)/t_i(b)*z_b-eta_b
% -1+xi_n(b)+c1*(eta_b-eta)
% -res.sol.itr.y(79)
% eta_b

    

% if all((alpha-(jj-1)*maxrange)<=0) || all((jj*maxrange-alpha)<=0)
    


