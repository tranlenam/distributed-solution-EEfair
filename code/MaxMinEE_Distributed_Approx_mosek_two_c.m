clear;
clc;
close all
%% Initial
nBS = 3 ;%number of BS
nUsers = 2 ;%number of MS
nTx =4 ;%number of BS antennas
%
P(1) = 10^(35/10)*10^-3; %[W] BS power constraint
P(2) = 10^(33/10)*10^-3; %[W] circuit power consumption
PAeff = 0.35; % power amplifier efficiency
No = -164;
P(3) = 10^(40/10)*10^-3;
snrthres=10^(0.0);

Pcon = P(1);
Psta = P(2);
Pdyn = P(3);
for iBS=1:nBS
    for iUser=1:nUsers
        BD(:,iUser,iBS)= myvec(kron(eye(nBS),ones(nUsers,1)));
        BD((iBS-1)*nUsers*(1+nBS)+iUser,iUser,iBS)= 0;
    end
end



load('../data/ChannelRealization_20.mat');

channel= channel*scalefactor;
scalefactor = sqrt(10^((-No)/10))/scalefactor;
power = P(1)*scalefactor^2;
[beamformerinit] = Initial_Beamformer(nBS,nUsers,nTx,channel,scalefactor,power,snrthres);

%% Centralized solutions
[ EEbeamforming_SCA, beamformer_SCA,iIteration_SCA,time_solving,res2,eta_centralized,...
    no_flops,beamformer,mybeta_i,z,t,eta,sumrate] = MaxMinEEBS_SCA_mosek_approx(nBS, nTx, ...
    nUsers,P,PAeff,No, channel, scalefactor, beamformerinit,power);


%% Distributed Solutions

%% Initialization

BD = zeros(nUsers*nBS*nBS,nUsers,nBS);
matrx =  myvec(kron(eye(nBS),ones(nUsers,1)));

for iBS=1:nBS
    for iUser=1:nUsers
        BD(:,iUser,iBS)= myvec(kron(eye(nBS),ones(nUsers,1)));
        BD((iBS-1)*nUsers*(1+nBS)+iUser,iUser,iBS)= 0;
    end
end



tau_i=zeros(nBS,nUsers,nBS);
tau_pr_i=zeros(nBS,nUsers,nBS);
q = zeros(nUsers,nBS);
for iBS=1:nBS
    for iUser = 1:nUsers
        myalpha(iUser,iBS)= norm(channel(:,iUser,iBS,iBS)'*beamformerinit(:,iUser,iBS))^2;
        q_i(iUser,iBS) =  real(norm([myvec((reshape(myvec(channel(:,iUser,iBS,:)),nTx,nBS)'*reshape(myvec(beamformerinit),nTx,nBS*nUsers))').*BD(:,iUser,iBS);1])^2);
        for jBS=1:nBS
            if iBS ~= jBS
                tau_i(jBS,iUser,iBS) = norm(channel(:,iUser,iBS,jBS)'*beamformerinit(:,1:nUsers,jBS))^2;
                tau_pr_i(iBS,iUser,jBS) = tau_i(jBS,iUser,iBS);
            end
        end
        q(iUser,iBS) = norm([1 channel(:,iUser,iBS,iBS)'*beamformerinit(:,1:nUsers~=iUser,iBS)])^2;
        g(iUser,iBS) = myalpha(iUser,iBS)/q(iUser,iBS);
        mytheta(iUser,iBS) = log(1+g(iUser,iBS));
    end
    z_i(iBS,1) = sqrt(sum(mytheta(:,iBS)));
    t_i(iBS,1) = Psta + Pdyn*nTx + (1/PAeff)*norm(myvec(beamformerinit(:,:,iBS)/scalefactor))^2; % zeta
end
for b=1:nBS
    Theta_n(:,b) = [myvec(tau_i(b,:,1:nBS~=b)) ; zeros(nUsers*(nBS-1),1)]; % no ICI information from the other nodes
end
% mybeamformer = beamformerinit;
q_i=q;
beamformer_i=beamformerinit;
% initial local variable
eta_n = (z_i.^2./t_i);
Nu=zeros((nBS-1)*nUsers*2,nBS);
% penalty parameter
c1=10;
c2=0.2;
% Lagrangian multiplier
xi_n=0*ones(nBS,1);
zeta_n=0*ones((nBS-1)*nUsers*2,nBS);

% storage vector
eta_sca=0;
obval=[];
obval1=[];
SINR = zeros(nUsers,nBS);
rate = zeros(nUsers,nBS);
EEbeamforming = zeros(1,nBS);
for iBS=1:nBS
    for iUser = 1:nUsers
        SINR(iUser,iBS) = (abs(channel(:,iUser,iBS,iBS)'*beamformerinit(:,iUser,iBS)))^2/...
            real(norm([myvec((reshape(myvec(channel(:,iUser,iBS,:)),nTx,nBS)'...
            *reshape(myvec(beamformerinit),nTx,nBS*nUsers))').*BD(:,iUser,iBS);1])^2); % calculates beta (=IUI plus noise power)
        rate(iUser,iBS) = log2(1+SINR(iUser,iBS));
    end
    EEbeamforming(iBS)= sum(rate(:,iBS))/(Psta + Pdyn*nTx + (1/PAeff)*norm(myvec(beamformerinit(:,:,iBS)/scalefactor))^2); % zeta

end
%  EEbeamforming
iIteration = 1;
EE=[];
% EE=[EE;EEbeamforming];
EE1=[];
% EE1=[EE1;EEbeamforming];


% nIteration =250;
Max_iteration=3000; % I_ADMM , I_ADMM is effective if set small

ee=[];
eta_i=[];
delta=[];
while(1) %for iIteration =1:nIterations

    % solve distributed system

    iIteration;
    % if iIteration >1
    %% Global update (using optimization method)
    [eta, mu1,Nu,W,optval]=MaxMin_ADMM_mosek_global_two_c(nUsers,nBS,zeta_n,xi_n,Theta_n,eta_n,c1,c2);
    % remark: eta, mu1, Nu can be calculated using close-form expression
    % provided in the paper

    ee=[ee; eta];
    %% local update
    valtemp = [];
    for b=1:nBS % for each BS
        [eta_b, beamformer_b, z_b, t_b, q_b, tau, tau_pr,Theta,W,obj,solTime_b,noflop_b,res]=MaxMinEEBS_ADMM_mosek_two_c(channel, power,Pdyn,Psta, scalefactor,beamformer_i,q_i,z_i,t_i,eta, Nu(:,b), b, nBS, nUsers,nTx,PAeff, xi_n, zeta_n, c1,c2);
        eta_n(b,1)=eta_b; % vector  eta_b = [eta_1, eta_2,...,eta_B]
        z(b,1)=z_b;       % vector  z = [z_1, z_2,...,z_B]
        t(b,1)=t_b;       % vector  t = [t_1, t_2,...,t_B]
        tau_pr_n(:,:,b)=tau_pr; % tau_prime
        tau_n(b,:,:) = tau';    % tau
        q(:,b)=q_b;       % vector  q = [q_1, q_2,...,q_B]
        beamformer(:,:,b)= beamformer_b; % vector w_b
        Theta_n(:,b)=Theta; % arranged vector \theta_b
    end
    eta_i=[eta_i;eta_n' eta];

    %% update Lagrangian multiplier
    % compute
    zeta_n = zeta_n + c2*(Theta_n - Nu);
    xi_n = xi_n + c1*(eta_n-eta);


    %% calculate the objective
    delta=[delta ; sum(-eta_n)+sum(myvec(zeta_n).*myvec((Theta_n-Nu)))+sum(xi_n.*(eta_n-eta))+c1/2*norm(eta_n-eta)^2+c2/2*norm(myvec(Theta_n-Nu))^2];

    % eta
    %% recalculate actual SINR, rate and EE
    for iBS=1:nBS
        for iUser = 1:nUsers
            SINR(iUser,iBS) = (norm(channel(:,iUser,iBS,iBS)'*beamformer(:,iUser,iBS)))^2/...
                real(norm([myvec((reshape(myvec(channel(:,iUser,iBS,:)),nTx,nBS)'*reshape(myvec(beamformer),nTx,nBS*nUsers))').*BD(:,iUser,iBS);1])^2); % calculates beta (=IUI plus noise power)
            rate(iUser,iBS) = log2(1+SINR(iUser,iBS));

        end
        EEbeamforming(iBS)= sum(rate(:,iBS))/(Psta + Pdyn*nTx + (1/PAeff)*norm(myvec(beamformer(:,:,iBS)/scalefactor))^2); % zeta
        sumrate(iBS) = sum(rate(:,iBS));
    end
    EEbeamforming;
    EE=[EE; EEbeamforming];

    %%
    if iIteration >1
        if (iIteration == Max_iteration) || norm(eta_i(iIteration,1:nBS)-eta_i(iIteration-1,1:nBS))<=1e-5

            Max_iteration = 3000 + iIteration;

            % udate SCA
            z_i = z;
            t_i = t;
            beamformer_i = beamformer;
            q_i=q;

            EE1=[EE1;EEbeamforming iIteration];
            eta_sca=[eta_sca; eta];
            %                 c2=c2+0.0005;
            %                 c1=c1/2;c2=c2/2;
            if abs(mean(eta_sca(end))-eta_sca(end-1))<1e-5
                break;
            end
        end
    end

    %%
    if  (iIteration==20000)% | ((eta >= 0.9*eta_centralized(end)) & (delta(iIteration) <= - 0.7264) )
        break;
    end

    iIteration = iIteration + 1;

end
disp(strcat('Centralized EE solution:   ', num2str(eta_centralized(end))));
disp(strcat('Distributed EE solution:   ', num2str(eta)))
disp(strcat('No of iterations:   ', num2str(iIteration)))
figure
% plot sequence of the objective
semilogx(1:length(delta),delta)
xlabel('Iteration index');
ylabel('Objective value');
title('Convergence of the objective');
saveas(gcf,'../results/convergence.png')