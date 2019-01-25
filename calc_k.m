function [ k,dKdp,Revs,dRevsdp ] = calc_k( model,x )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

nk = model.p.nk;
k = zeros(nk,1);
x1 = x;
ne = cell2mat(model.ensemble.ne);
nei = cell2mat(model.ensemble.nei);
nvr = cell2mat(model.ensemble.nvr);
nf = ne-1;
np = nf+nei+nvr;
nr = length(ne);


nu = model.ensemble.nu;
%u = x1(1:nu);
%x1(1:nu) = [];
dKdp = zeros(nk,sum(np));
Revs = zeros(nr,1);
dRevsdp = zeros(nr,sum(np));

% main loop to compute k

for i = 1:nr
    nk = sqrt(length(model.ensemble.dELdp{i}(:,1)));
    xr = x1(1:np(i));
    x1(1:np(i)) = [];
    pos = 0;
    V = model.ensemble.V{i};%*u;
    %npr = np(i);
    if model.p.exch(i)
        %npr = npr+1;
    %if model.p.exch(i)&&~model.p.sub(i)
    if V(1)>0
        xr = [xr;0];
    else
        xr = [xr;abs(V(1))];
    end
    end
    p = model.ensemble.N{i}*xr;
    p(1) = p(1)+1;
    npr = length(p);
    %IL = model.ensemble.IL{i};
    %IR = model.ensemble.IR{i};
    %dRLdp = model.ensemble.dRLdp{i};
    dELdp = model.ensemble.dELdp{i};
    %dELdp = [sparse(nk*nk,nu),dELdp];
    dELdp = dELdp;
    dVRdp = model.ensemble.dVRdp{i};
    %dVRdp = [sparse(nk,nu),dVRdp];
    dVRdp = dVRdp;
    %dVnetdp = [model.ensemble.V{i},sparse(nk,npr)];
    dVnetdp = sparse(nk,npr);
    %p = [u;p];
    %RL = IL+(reshape(dRLdp*p,nk,nk));
    EL = reshape(dELdp*p,nk,nk);
    VR = (reshape(dVRdp*p,nk,1));
    K = EL\(VR+V);
    DK = reshape(dELdp,nk,nk*(npr));
    DK = mat2cell(DK,nk,nk*ones(1,npr));
    DK = cell2mat(DK');
    DK = DK*K;
    DK = reshape(DK,nk,(npr));
    dkdp = EL\(dVnetdp + dVRdp - DK);
    %dkdu = EL\model.ensemble.V{i};
    %dKdu(model.p.kblocks(i)+1:model.p.kblocks(i+1),:) = dkdu;
    if nei(i)>0
        %dEildp = [sparse(nei(i)*nei(i),nu),model.ensemble.dEildp{i}];
        dEildp = model.ensemble.dEildp{i};
        Eil = reshape(dEildp*p,nei(i),nei(i));
        %dEirdp = [sparse(nei(i),nu),model.ensemble.dEirdp{i}];
        dEirdp = model.ensemble.dEirdp{i};
        Eir = dEirdp*p;
        Ki = Eil\Eir;
        DK = reshape(dEildp,nei(i),nei(i)*(npr));
        DK = mat2cell(DK,nei(i),nei(i)*ones(1,npr));
        DK = cell2mat(DK');
        DK = DK*Ki;
        DK = reshape(DK,nei(i),npr);
        dKidp = Eil\(dEirdp-DK);
        K = [K;Ki];
        dkdp = [dkdp;dKidp];
    end
    %dkdu = dkdp(:,1:nu);
    %dkdp(:,1:nu) = [];
    dkdp = dkdp*model.ensemble.N{i};
    %dkdp = [dkdu,dkdp];
    %account for irreversible last step in exchange reactions
    if model.p.exch(i)
        dkdp(:,end) = [];
    end
    k(model.p.kblocks(i)+1:model.p.kblocks(i+1)) = K;
    %dKdp(model.p.kblocks(i)+1:model.p.kblocks(i+1),1:nu) = dkdu;
    dKdp(model.p.kblocks(i)+1:model.p.kblocks(i+1),model.ensemble.ppos{i}) = dkdp;
    
    %computing reversibilities
    vrefind = find(any(model.ensemble.V{i},2));
    if isempty(vrefind)
        vnet = 0;
    else
        vnet = V(vrefind(1));
    end
    if ~model.p.exch(i)
        vx = 1:ne(i);
        vxr = 2*vx;
        vxf = vxr-1;
        dVrdp = dVRdp(vxr,:);
        dVfdp = dVRdp(vxf,:);% + dVnetdp(vxf,:);
        Vr = dVrdp*p;
        Vf = dVfdp*p;
        Vf = Vf + vnet;
        D1 = prod(Vf);
        D2 = prod(Vr);
        D1p = D1./Vf;
        D2p = D2./Vr;
        sVF = diag(D1p)*dVfdp;
        sVr = diag(D2p)*dVrdp;
        dPVfdp = sum(sVF,1);
        dPVrdp = sum(sVr,1);

        if vnet<0
            R = D1/D2;
            dRdp = (dPVfdp-(R*dPVrdp))/D2;
        else
            R = D2/D1;
            dRdp = (dPVrdp-(R*dPVfdp))/D1;
        end

        Revs(i) = R;
        %drdu = dRdp(1:nu);
        %dRdp(1:nu) = [];
        dRdp = dRdp*model.ensemble.N{i};
        %dRevsdp(i,1:nu) = drdu;
        dRevsdp(i,model.ensemble.ppos{i}) = dRdp;
    else
        Revs(i) = 0;
        dRevsdp(i,model.ensemble.ppos{i}) = 1;
    end
        
end
    



end

