function [EEbeamforming,beamformer,iIteration,total_time,res, tau,no_flops, ...
    beamformer_i,mybeta_i,z_i,t_i,eta,sumrate]=MaxMinEEBS_SCA_mosek_approx(nBS, nTx, nUsers,P,PAeff,...
    No, channel, scalefactor, beamformerinit,power)
Psta = P(2);
Pdyn = P(3);
BD = zeros(nUsers*nBS*nBS,nUsers,nBS);
for iBS=1:nBS
    for iUser=1:nUsers
        BD(:,iUser,iBS)= myvec(kron(eye(nBS),ones(nUsers,1)));
        BD((iBS-1)*nUsers*(1+nBS)+iUser,iUser,iBS)= 0;
    end
end

mybeta = zeros(nUsers,nBS);

myalpha = zeros(nUsers,nBS);


for iBS=1:nBS
    for iUser = 1:nUsers
        myalpha(iUser,iBS)= abs(channel(:,iUser,iBS,iBS)'*beamformerinit(:,iUser,iBS))^2;
        mybeta(iUser,iBS) =  real(norm([myvec((reshape(myvec(channel(:,iUser,iBS,:)),nTx,nBS)'*reshape(myvec(beamformerinit),nTx,nBS*nUsers))').*BD(:,iUser,iBS);1])^2);
        mygamma(iUser,iBS) = myalpha(iUser,iBS)/mybeta(iUser,iBS);
        mytheta(iUser,iBS) = log(1+mygamma(iUser,iBS));
    end
    z(iBS,1) = sqrt(sum(mytheta(:,iBS)));
    t(iBS,1) = Psta + Pdyn*nTx + (1/PAeff)*norm(myvec(beamformerinit(:,:,iBS)/scalefactor))^2; % zeta
end
beamformer_i=beamformerinit;



mybeta_i = real(mybeta);
z_i = z;
t_i = t;
nIterations =40;
uptime=[];
seqobj = zeros(nIterations+1,1);
for iBS=1:nBS
    for iUser = 1:nUsers
        SINR(iUser,iBS) = (abs(channel(:,iUser,iBS,iBS)'*beamformerinit(:,iUser,iBS)))^2/...
            real(norm([myvec((reshape(myvec(channel(:,iUser,iBS,:)),nTx,nBS)'*reshape(myvec(beamformerinit),nTx,nBS*nUsers))').*BD(:,iUser,iBS);1])^2); % calculates beta (=IUI plus noise power)
        rate(iUser,iBS) = log2(1+SINR(iUser,iBS));
    end
    EEbeamforming(1,iBS)= sum(rate(:,iBS))/(Psta + Pdyn*nTx + (1/PAeff)*norm(myvec(beamformerinit(:,:,iBS)/scalefactor))^2); % zeta

end
% EEbeamforming(1) = computeEE_real(channel,beamformerinit,Psta + Pdyn*nTx,PAeff,beamformerinit/scalefactor); %% EE of the initial iteration
iIteration = 1;
seqobj(1) = min(t);
m=10;
maxrange=8;
tau(1)=0;
time_solving = [];
no_flops = 0;
Z=[];
T=[];
Beta=[];
%%

varlength = 2*nBS*nUsers*nTx+(nUsers*nBS-1)*nUsers*nBS*2+nBS+nBS*nUsers+2*nUsers*nBS+nBS*3+2*nUsers*nTx*nBS+nUsers*nBS+nBS+nBS+nBS+2*nBS*nUsers+3*nUsers*nBS+3*(m+3)*nUsers*nBS +1;

nlength3=2*nBS*nUsers*nTx+nBS;% nBS due to adding nBS variables for sqrt(power) in conic constraints
%% norm(w)< power
for i=1:nBS
    socp_con{i}.type = 'MSK_CT_QUAD';
    socp_con{i}.sub = [2*nBS*nUsers*nTx+i,(i-1)*2*nTx*nUsers+1:(i)*2*nTx*nUsers];
end
% cindex3 = size(a1,1)+size(a2,1)+size(a3,1);
%% q_k>= |SINR_kb|^2
a4=[];
for iBS=1:nBS
    for iUser=1:nUsers
        for jj=1:nBS
            for i=1:nUsers
                if jj==iBS
                    if i~=iUser
                        temp = zeros(2,varlength);
                        temp(1,(jj-1)*nUsers*2*nTx+(i-1)*2*nTx+1:(jj-1)*nUsers*2*nTx+(i)*2*nTx )...
                            =[ myvec([ real(channel(:,iUser,iBS,jj))';  imag(channel(:,iUser,iBS,jj))'])'];
                        %                     temp(1,2*nTx*nBS*nUsers+nBS*nUsers+(iBS-1)*(nUsers*nBS-1)+(jj-1)*nUsers+i)=-1;
                        temp(2,(jj-1)*nUsers*2*nTx+(i-1)*2*nTx+1:(jj-1)*nUsers*2*nTx+(i)*2*nTx )...
                            =[ myvec([-imag(channel(:,iUser,iBS,jj))';  real(channel(:,iUser,iBS,jj))'])'];
                        a4= [a4;temp];
                    end
                else
                    temp = zeros(2,varlength);
                    temp(1,(jj-1)*nUsers*2*nTx+(i-1)*2*nTx+1:(jj-1)*nUsers*2*nTx+(i)*2*nTx )...
                        =[ myvec([ real(channel(:,iUser,iBS,jj))';  imag(channel(:,iUser,iBS,jj))'])'];
                    temp(2,(jj-1)*nUsers*2*nTx+(i-1)*2*nTx+1:(jj-1)*nUsers*2*nTx+(i)*2*nTx )...
                        =[ myvec([-imag(channel(:,iUser,iBS,jj))';  real(channel(:,iUser,iBS,jj))'])'];
                    a4= [a4;temp];
                end
            end
        end
    end
end
for i=1:size(a4,1)
    a4(i,nlength3+i)=-1;
end

nlength4=nlength3+(nUsers*nBS-1)*nUsers*nBS*2+nBS*nUsers; %nBS*nUsers for vars x=1;
for i=1:nUsers*nBS
    socp_con{nBS+i}.type = 'MSK_CT_RQUAD';
    socp_con{nBS+i}.sub = [nlength4+i,nlength4+nUsers*nBS+i,nlength3+2*(i-1)*(nUsers*nBS-1)+1:nlength3+2*i*(nUsers*nBS-1),nlength4-nBS*nUsers+i];
end
nlength4=nlength4+2*nUsers*nBS;%nUsers*nBS vars for q_kb and nUsers*nBS x=1/2;
%% PAeff*(t_b-P0)*scalefactor^2=2x => PAeff*scalefactor^2*tb-2x=P0*PAeff*scalefactor^2

a5=[];
for iBS=1:nBS
    %     for iUser=1:nUsers
    temp=zeros(1,varlength);
    temp([nlength4+iBS,nlength4+nBS+iBS])=[PAeff*scalefactor^2 -2];
    a5=[a5;temp];
end
nlength5=nlength4+nBS*3;% t_b and 2x and 1;
%% 2x*1>=sum(w)
% x6=w
a6=[];
for i=1:nBS
    for ii=1:2*nUsers*nTx
        temp = zeros(1,varlength);
        temp([(i-1)*2*nUsers*nTx+ii,nlength5+(i-1)*2*nUsers*nTx+ii])=[1 -1];
        a6=[a6;temp];
    end
end
% 2x*1>=sum(w)
for i=1:nBS
    socp_con{nBS+nUsers*nBS+i}.type = 'MSK_CT_RQUAD';
    socp_con{nBS+nUsers*nBS+i}.sub = [nlength5-2*nBS+i,nlength5-nBS+i,nlength5+(i-1)*2*nUsers*nTx+1:nlength5+(i)*2*nUsers*nTx];
end
nlength6=nlength5+2*nUsers*nTx*nBS;%w
% cindex6=size(a1,1)+size(a2,1)+size(a3,1)+size(a4,1)+size(a5,1)+size(a6,1);
%% sum(log(1+g_kb) >= z_b^2 => sum(x7)>=z^2; 1+g_b >= exp(x7)
% sum(x7)=x71+...+x7n=2x7'
a7=zeros(nBS,varlength);
for iBS=1:nBS
    a7(iBS,[nlength6+(iBS-1)*nUsers+1:nlength6+iBS*nUsers nlength6+nUsers*nBS+iBS])=[ones(1,nUsers) -2];
end
nlength7=nlength6+nUsers*nBS+nBS+nBS+nBS;%nUsers*nBS for x7, nBS for 2x7', nBS for 1, nBS for z_b
% cindex7=cindex6+size(a7,1);
%2x7'*1>=z^2
for i=1:nBS
    socp_con{nBS+nUsers*nBS+nBS+i}.type = 'MSK_CT_RQUAD';
    socp_con{nBS+nUsers*nBS+nBS+i}.sub = [nlength6+nUsers*nBS+i,nlength6+nUsers*nBS+nBS+i,nlength6+nUsers*nBS+2*nBS+i];
end

%%
% 1+g_kb>=exp(x7);
% 1+g_kb>=x8 => g_kb-x8>=-1
a8=zeros(1,varlength);
for i=1:nBS*nUsers
    a8(i,[nlength7+i nlength7+nBS*nUsers+i])=[1 -1];
end
nlength8=nlength7+2*nBS*nUsers;%g_b x8
%% x8 >=exp(x7)
% m=10;

tic
conelength=nBS+nUsers*nBS+nBS+nBS;
a9=zeros(1,varlength);
a10=[];
a11=[];
a12=zeros((m+3)*nUsers*nBS,varlength);
a12(:,[nlength8+3*nUsers*nBS+1:nlength8+3*nUsers*nBS+nUsers*nBS*(m+3), nlength8+3*nUsers*nBS+2*nUsers*nBS*(m+3)+1:nlength8+3*nUsers*nBS+3*nUsers*nBS*(m+3)])= [eye(nUsers*nBS*(m+3)) -eye(nUsers*nBS*(m+3))];
for i=1:nBS*nUsers
    a9(i,[nlength7+nBS*nUsers+i nlength8+i])=[1 -1]; %x8-t>=0; t=nlength8+1:nlength8+nUsers*nBS
    for ii=1:m+3
        if ii<5
            if ii==1
                temp=zeros(1,varlength);
                temp([nlength6+i nlength8+nUsers*nBS+i])=[1/2^m -1];%x7/2^m-x9_1=-1; x9_1=nlength8+nUsers*nBS+1:nlength8+2*nUsers*nBS
                socp_con{conelength+(i-1)*(m+3)+ii}.type = 'MSK_CT_RQUAD';
                socp_con{conelength+(i-1)*(m+3)+ii}.sub = [nlength8+3*nUsers*nBS+(i-1)*(m+3)+ii,nlength8+3*nUsers*nBS+nUsers*nBS*(m+3)+(i-1)*(m+3)+ii,nlength8+nUsers*nBS+i];% u1_kb * (2x=1) >=x9_kb^2
                a10=[a10;temp];
            end
            if ii==2
                temp=zeros(1,varlength);
                temp([nlength6+i nlength8+2*nUsers*nBS+i])=[1/2^(m+1) -1];%x7/2^(m+1)x9_2=-5/6; x9_2=length8+2*nUsers*nBS+1:length8+3*nUsers*nBS
                socp_con{conelength+(i-1)*(m+3)+ii}.type = 'MSK_CT_RQUAD';
                socp_con{conelength+(i-1)*(m+3)+ii}.sub = [nlength8+3*nUsers*nBS+(i-1)*(m+3)+ii,nlength8+3*nUsers*nBS+nUsers*nBS*(m+3)+(i-1)*(m+3)+ii,nlength8+2*nUsers*nBS+i];% u2_kb * (2x=1) >=x9_kb^2
                a10=[a10; temp];
            end
            if ii==3
                %                 a10((i-1)*(m+3)+ii,[nlength6+i nlength8+nUsers*nBS+(i-1)*nUsers*nBS+ii])=[0 0];% nothing
                socp_con{conelength+(i-1)*(m+3)+ii}.type = 'MSK_CT_RQUAD';
                socp_con{conelength+(i-1)*(m+3)+ii}.sub = [nlength8+3*nUsers*nBS+(i-1)*(m+3)+ii,nlength8+3*nUsers*nBS+nUsers*nBS*(m+3)+(i-1)*(m+3)+ii,nlength8+3*nUsers*nBS+2*nUsers*nBS*(m+3)+(i-1)*(m+3)+1];% u3_kb * (2x=1) >=u1_kb^2
            end
            if ii==4
                temp=zeros(1,varlength);
                temp([nlength8+3*nUsers*nBS+(i-1)*(m+3)+2, nlength8+3*nUsers*nBS+(i-1)*(m+3)+3, nlength8+3*nUsers*nBS+(i-1)*(m+3)+4])=[1 1/24 -1];% u2_kb+1/24*u3_kb-u4 <= -19/72
                a11=[a11;temp];
            end
            % note:
            % u= nlength8+3*nUsers*nBS+1:nlength8+3*nUsers*nBS+(m+3)*nUsers*nBS
            % 2x=nlength8+3*nUsers*nBS+(m+3)*nUsers*nBS+1:nlength8+3*nUsers*nBS+2*(m+3)*nUsers*nBS
            % u'=u=nlength8+3*nUsers*nBS+2*(m+3)*nUsers*nBS+1:nlength8+3*nUsers*nBS+3*(m+3)*nUsers*nBS
        else
            %a10((i-1)*(m+3)+ii,:)=0;% nothing
            socp_con{conelength+(i-1)*(m+3)+ii-1}.type = 'MSK_CT_RQUAD';
            socp_con{conelength+(i-1)*(m+3)+ii-1}.sub = [nlength8+3*nUsers*nBS+(i-1)*(m+3)+ii,nlength8+3*nUsers*nBS+nUsers*nBS*(m+3)+(i-1)*(m+3)+ii-1,nlength8+3*nUsers*nBS+2*nUsers*nBS*(m+3)+(i-1)*(m+3)+ii-1];% u_i_kb * (2x=1) >=u_i-1_kb^2
        end
    end
    socp_con{conelength+(i-1)*(m+3)+ii}.type = 'MSK_CT_RQUAD';
    socp_con{conelength+(i-1)*(m+3)+ii}.sub = [nlength8+i,nlength8+3*nUsers*nBS+nUsers*nBS*(m+3)+(i-1)*(m+3)+ii,nlength8+3*nUsers*nBS+2*nUsers*nBS*(m+3)+(i-1)*(m+3)+ii];% t* (2x=1) >=u_m+3_kb^2
end
nlength10=nlength8+3*nUsers*nBS+3*(m+3)*nUsers*nBS; %t,x9_1, x9_2 , 2* u, 1
conelength = conelength + (m+3)*nUsers*nBS;
bound10=myvec([-ones(1,nUsers*nBS);-5/6*ones(1,nUsers*nBS)])';
%

while(1)
    %% phi_b = ()z_b -()t_b - tau>= 0
    a13=zeros(nBS,varlength);
    for iBS=1:nBS
        a13(iBS,[nlength4+iBS,nlength7-nBS+iBS,varlength])=[-(z_i(iBS)/t_i(iBS))^2, 2*z_i(iBS)/t_i(iBS),-1];
    end
    nlength12=nlength10;
    % cindex11=cindex9+size(a10,1);
    %% psi=2real()-()*q_kb>=g_kb
    %
    a14=[];
    for iBS=1:nBS
        for iUser=1:nUsers
            whh = beamformer_i(:,iUser,iBS)'*channel(:,iUser,iBS,iBS)*channel(:,iUser,iBS,iBS)';
            temp=zeros(1,varlength);
            temp((iBS-1)*nUsers*2*nTx+(iUser-1)*2*nTx+1:(iBS-1)*nUsers*2*nTx+(iUser)*2*nTx ) =2/mybeta_i(iUser,iBS)*[ myvec([ real(whh); -imag(whh)])'];
            temp(nlength4-2*nUsers*nBS+(iBS-1)*nUsers+iUser)= -(abs(channel(:,iUser,iBS,iBS)'*beamformer_i(:,iUser,iBS))/mybeta_i(iUser,iBS))^2;
            temp(nlength7+(iBS-1)*nUsers+iUser)=-1;
            a14=[a14;temp];
        end
    end
    nlength14=nlength12;
    %%
    clear prob;
    prob.c=[zeros(1,varlength-1) 1];
    prob.a=sparse([a4;a5;a6;a7;a8;a9;a10;a11;a12;a13;a14]);
    prob.blc= [zeros(1,size(a4,1))    (Psta + Pdyn*nTx)*PAeff*scalefactor^2*ones(1,size(a5,1)) 0*ones(1,size(a6,1)) zeros(1,size(a7,1)) -1*ones(1,size(a8,1)) zeros(1,size(a9,1)) bound10 -inf*ones(1,size(a11,1)) zeros(1,size(a12,1)) zeros(1,size(a13,1)) zeros(1,size(a14,1)) ];
    prob.buc= [zeros(1,size(a4,1)) (Psta + Pdyn*nTx)*PAeff*scalefactor^2*ones(1,size(a5,1)) zeros(1,size(a6,1)) inf*ones(1,size(a7,1)) inf*ones(1,size(a8,1)) inf*ones(1,size(a9,1)) bound10 -19/72*ones(1,size(a11,1)) 0*ones(1,size(a12,1))  inf*ones(1,size(a13,1)) inf*ones(1,size(a14,1)) ];
    prob.blx= [-inf*ones(1,2*nTx*nUsers*nBS),... % w
        sqrt(power)*ones(1,nBS),... % sqrt(power)
        -inf*ones(1,(nUsers*nBS-1)*nUsers*nBS*2),... % interference
        ones(1,nBS*nUsers),... % 1
        zeros(1,nUsers*nBS),... % q_kb
        0.5*ones(1,nUsers*nBS),... % 1/2
        zeros(1,nBS) zeros(1,nBS) ones(1,nBS),... % t_b, 2x, 1
        -inf*ones(1,2*nUsers*nTx*nBS),... % x=w
        -inf*ones(1,nUsers*nBS) -inf*ones(1,nBS),... % x7, x7'
        ones(1,nBS) zeros(1,nBS),... % 1, z_b
        -ones(1,nUsers*nBS) zeros(1,nUsers*nBS),... % g_kb, x8
        -inf*ones(1,nUsers*nBS),... % t
        -inf*ones(1,2*nUsers*nBS),... % x9_1, x9_2
        -inf*ones(1,(m+3)*nUsers*nBS),... % u
        0.5*ones(1,(m+3)*nUsers*nBS),... % 1/2
        -inf*ones(1,(m+3)*nUsers*nBS),... % u'
        0      ]; %tau
    prob.bux= [inf*ones(1,2*nTx*nUsers*nBS),... % w
        sqrt(power)*ones(1,nBS),... % sqrt(power)
        inf*ones(1,(nUsers*nBS-1)*nUsers*nBS*2),... % interference
        ones(1,nBS*nUsers),... % 1
        inf*ones(1,nUsers*nBS),... % q_kb
        0.5*ones(1,nUsers*nBS),... % 1/2
        inf*ones(1,nBS) inf*ones(1,nBS) ones(1,nBS),... % t_b, 2x, 1
        inf*ones(1,2*nUsers*nTx*nBS),... % x=w
        inf*ones(1,nUsers*nBS) inf*ones(1,nBS),... % x7, x7'
        ones(1,nBS) inf*ones(1,nBS),... % 1, z_b
        inf*ones(1,nUsers*nBS) inf*ones(1,nUsers*nBS),... % g_kb, x8
        inf*ones(1,nUsers*nBS),... % t
        inf*ones(1,2*nUsers*nBS),... % x9_1, x9_2
        inf*ones(1,(m+3)*nUsers*nBS),... % u
        0.5*ones(1,(m+3)*nUsers*nBS),... % 1/2
        inf*ones(1,(m+3)*nUsers*nBS),... % u'
        inf]; % tau

    prob.cones=socp_con;
    %        2*nBS*nUsers*nTx+nBS*nUsers+(nUsers*nBS-1)*nUsers*nBS*2+nBS*nUsers+nBS*nUsers+nUsers*nBS+nUsers*nBS+nBS+2*nUsers*nBS+nUsers*nBS+nBS+1;
    % param=[];
    uptime_temp=toc;
    uptime=[uptime uptime_temp];
    % tic
    [r,res]       = mosekopt('maximize  echo(0) info',prob);
    % toc
    W = res.sol.itr.xx;
    % [res.info.MSK_DINF_OPTIMIZER_TIME res.info.MSK_DINF_SIM_TIME]
    time_solving = [time_solving res.info.MSK_DINF_OPTIMIZER_TIME];
    no_flops1 =  res.info.MSK_DINF_INTPNT_FACTOR_NUM_FLOPS;
    for iBS=1:nBS
        for iUser=1:nUsers

            reBeamformer = W((iBS-1)*nUsers*2*nTx+(iUser-1)*2*nTx+[1:2:2*nTx-1],1);
            imBeamformer = W((iBS-1)*nUsers*2*nTx+(iUser-1)*2*nTx+[2:2:2*nTx],1);
            beamformer(:,iUser,iBS) = reBeamformer+j*imBeamformer;

        end
    end
    tau(iIteration+1) = res.sol.itr.pobjval;
    eta=res.sol.itr.pobjval;
    alpha = W(nlength6+1:nlength6+nBS*nUsers);
    z_i = W(nlength7-nBS+1:nlength7);
    t_i = W(nlength4+1:nlength4+nBS);
    mybeta_i = reshape(W(nlength4-2*nUsers*nBS+1:nlength4-nUsers*nBS),nUsers,nBS);
    beamformer_i = beamformer;

    for iBS=1:nBS
        for iUser = 1:nUsers
            SINR(iUser,iBS) = (abs(channel(:,iUser,iBS,iBS)'*beamformer(:,iUser,iBS)))^2/...
                real(norm([myvec((reshape(myvec(channel(:,iUser,iBS,:)),nTx,nBS)'*reshape(myvec(beamformer),nTx,nBS*nUsers))').*BD(:,iUser,iBS);1])^2); % calculates beta (=IUI plus noise power)
            rate(iUser,iBS) = log2(1+SINR(iUser,iBS));
        end
        EEbeamforming(iIteration+1,iBS)= sum(rate(:,iBS))/(Psta + Pdyn*nTx + (1/PAeff)*norm(myvec(beamformer(:,:,iBS)/scalefactor))^2); % zeta
        sumrate(iBS)=sum(rate(:,iBS));
    end
    %         EEbeamforming

    if iIteration >1
        if abs(tau(iIteration+1) - tau(iIteration)) <= 1e-6
            break
        end
    end
    %
    %         end
    iIteration = iIteration + 1;
    if iIteration == nIterations
        break
    end

    temp1=beamformer_i;
    temp2=mybeta_i;
    temp3=z_i;
    temp4=t_i;
    temp5=res.sol.itr.pobjval;


end
tau=tau';
total_time=(time_solving');
% if all((alpha-(jj-1)*maxrange)<=0) || all((jj*maxrange-alpha)<=0)



