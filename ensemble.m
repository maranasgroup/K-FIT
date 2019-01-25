function [model] = ensemble(model)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here



nr = length(model.rid);
Nflx = model.p.N;
nu = length(Nflx(1,:));
model.ensemble.Nflx = Nflx;
model.ensemble.nu = nu;
model.ensemble.ne = cell(nr,1);
model.ensemble.nei = cell(nr,1);
%model.ensemble.nrev = cell(nr,1);
model.ensemble.nvr = cell(nr,1);
%model.ensemble.IL = cell(nr,1);
%model.ensemble.dRLdp = cell(nr,1);
model.ensemble.dELdp = cell(nr,1);
%model.ensemble.dRRdp = cell(nr,1);
model.ensemble.dVRdp = cell(nr,1);
%model.ensemble.IR = cell(nr,1);
model.ensemble.V = cell(nr,1);
model.ensemble.dEildp = cell(nr,1);
model.ensemble.dEirdp = cell(nr,1);
model.ensemble.dEirdp = cell(nr,1);
model.ensemble.N = cell(nr,1);
model.ensemble.ppos = cell(nr,1);

csum = 0;
vref = model.d.flx{1};


for i = 1:nr
    ne = sqrt(length(model.p.dEdk{i}(:,1)));
    nk = length(model.kinetic(i).S(1,:));
    nei = 0;
    if ~isempty([model.kinetic(i).c_in,model.kinetic(i).uc_in,model.kinetic(i).nc_in])
        nei = length(model.kinetic(i).I(1,:));
    end
    %if vref(i)>=0
        rflag = true;
    %else
    %    rflag = false;
    %end
    %if vref(i) == 0
    %    vref(i) = 0.01;
    %end
    nvr = ne;
    np = ne+nei+nvr-1;
    nef = np-nvr;
    N = [-1*ones(1,nef);eye(nef)];
    
    %{
    % constructing RL
    dRLdr = zeros(nk*nk,nrev);
    IL = eye(nk);
    for j = 1:nrev
        dRLdr(nk*(2*j-2)+(2*j-1),j) = -1;
        dRLdr(nk*(2*j-1)+(2*j),j) = -1;
    end
    dRLdp = [zeros(nk*nk,(ne+nei)),dRLdr];
    %}
    %constructing EL
    dELdne = sparse(nk*nk,ne);
    for j = 1:ne
        dELdne(nk*(2*j-2)+(2*j-1),j) = 1;
        if j < ne
            dELdne(nk*(2*j-1)+(2*j),j+1) = 1;
        else
            dELdne(nk*(2*j-1)+(2*j),1) = 1;
        end
    end
    dELdp = [dELdne,sparse(nk*nk,nei+nvr)];
    
    %constructing VR
    dVRdr = sparse(nk,nvr);
    %IR = zeros(nk,nk);
    for j = 1:nvr
        %if ~rflag
            dVRdr((2*j-1),j) = 1;
            dVRdr((2*j),j) = 1;
            %IR(2*j,2*j) = 1;
        %{
        %else
            dRRdr(nk*(2*j-1)+(2*j),j) = 1;
            IR(2*j-1,2*j-1) = 1;
        end
            %}
    end
    dVRdp = [sparse(nk,(ne+nei)),dVRdr];
    
    %computing V
    V = vref(i)*ones(nk,1);
    
    %V = vref(logical(model.d.rmap{1}(:,i)))*ones(nk,1);
    %V = repmat(Nflx(i,:),nk,1);
    ind = 2:2:nk;
    %V(ind) = 0;
    V(ind,:) = 0;
    
    %matrices for computing Ki
    dEirdei = speye(nei);
    dEilde = sparse(nei*nei,ne);
    if nei>0
        for j = 1:nei
            dEilde(nei*(j-1)+j,:) = model.kinetic(i).I(:,j)';
        end
    end
    dEirdp = [sparse(nei,ne),dEirdei,sparse(nei,nvr)];
    dEildp = [dEilde,sparse(nei*nei,nei),sparse(nei*nei,nvr)];
    
    %collect and store required matrices
    
    model.ensemble.ne{i} = ne;
    model.ensemble.nei{i} = nei;
    %model.ensemble.IL{i} = IL;
    %model.ensemble.dRLdp{i} = dRLdp;
    model.ensemble.dELdp{i} = dELdp;
    model.ensemble.dVRdp{i} = dVRdp;
    %model.ensemble.IR{i} = IR;
    %model.ensemble.V{i} = abs(vref(i))*ones(nk,1);
    model.ensemble.V{i} = V;
    model.ensemble.dEildp{i} = dEildp;
    model.ensemble.dEirdp{i} = dEirdp;
    model.ensemble.N{i} = blkdiag(N,eye(nvr));
    if model.p.exch(i)
    %if model.p.exch(i)&&~model.p.sub(i)
        nvr = nvr-1;
    end
    model.ensemble.nvr{i} = nvr;
    np = ne+nei+nvr-1;
    model.ensemble.ppos{i} = csum+1:csum+np;
    csum = csum+np;
end
    
        

end

