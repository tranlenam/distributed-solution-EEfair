function [beamformer] = Initial_Beamformer(nBS,nUsers,nTx,channel,scalefactor,power,snrthres)

[r, res] = mosekopt('symbcon'); 

prob = [];
% varlength = 2*nBS*nUsers*nTx+nBS*nUsers+(nUsers*nBS-1)*nUsers*nBS*2+nBS+nUsers*nBS+(nUsers*nBS-1)*nUsers*nBS*2+2*nUsers*nBS+nUsers*nBS+nBS+nBS+1;
varlength = 2*nBS*nUsers*nTx+nBS*nUsers+(nUsers*nBS-1)*nUsers*nBS*2+nBS+nUsers*nBS;
a1=zeros(nUsers*nBS,varlength);
for iBS=1:nBS
    for iUser=1:nUsers
        a1((iBS-1)*nUsers+iUser,(iBS-1)*nUsers*2*nTx+(iUser-1)*2*nTx+1:(iBS-1)*nUsers*2*nTx+(iUser)*2*nTx ) =[ myvec([-imag(channel(:,iUser,iBS,iBS))'; real(channel(:,iUser,iBS,iBS))'])'];
    end
end
nlength1=2*nBS*nUsers*nTx;
% Real() >= a2
a2=zeros(nUsers*nBS,varlength);
for iBS=1:nBS
    for iUser=1:nUsers
        a2((iBS-1)*nUsers+iUser,(iBS-1)*nUsers*2*nTx+(iUser-1)*2*nTx+1:(iBS-1)*nUsers*2*nTx+(iUser)*2*nTx ) =[ myvec([ real(channel(:,iUser,iBS,iBS))'; imag(channel(:,iUser,iBS,iBS))'])'];
        a2((iBS-1)*nUsers+iUser,2*nTx*nBS*nUsers+(iBS-1)*nUsers+iUser)=-1; %
    end
end
%interference
nlength2=nlength1+nBS*nUsers;
a3=[];
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
                    a3= [a3;temp];  
                end
                else 
                    temp = zeros(2,varlength);
                    temp(1,(jj-1)*nUsers*2*nTx+(i-1)*2*nTx+1:(jj-1)*nUsers*2*nTx+(i)*2*nTx )...
                          =[ myvec([ real(channel(:,iUser,iBS,jj))';  imag(channel(:,iUser,iBS,jj))'])'];
                    temp(2,(jj-1)*nUsers*2*nTx+(i-1)*2*nTx+1:(jj-1)*nUsers*2*nTx+(i)*2*nTx )...
                          =[ myvec([-imag(channel(:,iUser,iBS,jj))';  real(channel(:,iUser,iBS,jj))'])'];
                    a3= [a3;temp];
                end
            end
        end
    end
end
for i=1:size(a3,1)
    a3(i,nlength2+i)=-1/sqrt(snrthres);
end

nlength3=2*nBS*nUsers*nTx+nBS*nUsers+(nUsers*nBS-1)*nUsers*nBS*2+nBS;% nBS due to adding nBS variables for sqrt(power) in conic constraints
nlength4=nlength3+nBS*nUsers; %nBS*nUsers x=sigma 


prob.c = [];
prob.a = sparse([a1;a2;a3]);

prob.blc = [zeros(1,size(a1,1)) zeros(1,size(a2,1)) zeros(1,size(a3,1)) ]; 
prob.buc = [zeros(1,size(a1,1)) inf*ones(1,size(a2,1)) zeros(1,size(a3,1))];
prob.blx = [-inf*ones(1,2*nBS*nUsers*nTx) 0*ones(1,nBS*nUsers), -inf*ones(1,(nUsers*nBS-1)*nUsers*nBS*2) sqrt(power)*ones(1,nBS) ones(1,nUsers*nBS) ]; 
prob.bux = [inf*ones(1,2*nBS*nUsers*nTx+nBS*nUsers+(nUsers*nBS-1)*nUsers*nBS*2) sqrt(power)*ones(1,nBS) ones(1,nUsers*nBS) ]; 

prob.cones = cell(nUsers*nBS+nBS,1);

%real() >= snrthress*sqrt(interference)
for i=1:nUsers*nBS
    prob.cones{i}.type = res.symbcon.MSK_CT_QUAD;
    prob.cones{i}.sub = [2*nBS*nUsers*nTx+i,2*nTx*nBS*nUsers+nUsers*nBS+2*(i-1)*(nUsers*nBS-1)+1:2*nTx*nBS*nUsers+nUsers*nBS+2*i*(nUsers*nBS-1),nlength3+i];
end
%norm(w)< power
for i=1:nBS
    prob.cones{nUsers*nBS+i}.type = res.symbcon.MSK_CT_QUAD;
    prob.cones{nUsers*nBS+i}.sub = [nlength3-nBS+i,(i-1)*2*nTx*nUsers+1:(i)*2*nTx*nUsers];
end


[r,res]       = mosekopt('maximize echo(0)',prob); 

W = res.sol.itr.xx;
for iBS=1:nBS
    for iUser=1:nUsers
        
            reBeamformer = W((iBS-1)*nUsers*2*nTx+(iUser-1)*2*nTx+[1:2:2*nTx-1],1);
            imBeamformer = W((iBS-1)*nUsers*2*nTx+(iUser-1)*2*nTx+[2:2:2*nTx],1);
            beamformer(:,iUser,iBS) = reBeamformer+j*imBeamformer;
        
    end
end
end
%%
% beamformerinit=mybeamformer;