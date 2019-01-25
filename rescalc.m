function [r,W,J,vop,cs] = rescalc(x,model)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
ncond = length(model.d.vpert(1,:));
cind = 69;
nr = length(model.p.dEdk);
v = zeros(nr,ncond);
sthresh = 1e-4;
[k,dkdp,revs,drevsdp] = calc_k(model,x);
nrevs = length(model.d.ridx);
if nrevs > 0
    rpred = revs(model.d.ridx);
    rmeas = model.d.revs;
    std = model.d.rerr;
    res = rpred-rmeas;
    drdp = drevsdp(model.d.ridx,:);
else
    res = zeros(nrevs,1);
    drdp = zeros(nrevs,length(drevsdp(1,:)));
    std = res;
end

model.p.k = k;
%tspan = [0,10000];
nm = length(model.metprop);
cs = zeros(nm,ncond);
dvdk = cell(1,ncond);
dcsdp = dvdk;
desdp = dcsdp;
nenz = model.ssubs.eblocks(end);
es = zeros(nenz,ncond);
%{
for i = 1:ncond
    done = false;
    c0 = ones(1,nm);
    ctr = 0;
    while ~done
        complete = false;
        css = c0;
        iter = 0;
        while ~complete
            iter = iter+1;
            [css,~,complete] = svsucsubs(css,model,i);
            if ~complete
                [css,complete] = SIFOInteg(css,model,i);
            %else
            %    css = ci;
            end
            if iter>1 && ~complete
                opt = odeset('Jacobian',@(t,x) jacfn(t,x,model,i));
                [~,C] = ode15s(@(t,x) svinteg(t,x,model,i),[0,1e4],css,opt);
                css = C(end,:)';
                [dc] = svinteg(1,css,model,i);
                if max(abs(dc))<1e-5
                    complete = true;
                end
            end
        end
        [dx,vx] = svinteg(1,css,model,i);
        if max(abs(dx))<1e-5
            done = true;
            [vx,dvxdk] = flxcalc(model,css);
            vx = vx.*model.d.vpert(:,i);
            dvxdk = diag(full(model.d.vpert(:,i)))*dvxdk;
            %vsc = 100*vx/abs(vx(12));
            %dvscdk = (100/abs(vx(12)))*(dvxdk-vsc*dvxdk(12,:));
            %Css(:,i) = c(end,:)';
        else
            ctr = ctr+1;
            c0 = css;
            sc = 100/abs(vx(12));
            if ctr>20 && max(abs(dx))*sc<0.1
                done = true;
                [vx,dvxdk] = flxcalc(model,css);
                vx = vx.*model.d.vpert(:,i);
                dvxdk = diag(full(model.d.vpert(:,i)))*dvxdk;
            end
                
        end
    end
    %vx = vx.*model.d.vpert(:,i);
    v(:,i) = vx;
    dvdk{i} = dvxdk;
end
%}
parfor i = 1:ncond
    done = false;
    c0 = ones(1,nm);
    ctr = 0;
    while ~done
        complete = false;
        %css = c0;
        iter = 0;
        while ~complete
            iter = iter+1;
            [cs(:,i),~,complete] = svsucsubs(c0,model,i,sthresh);
            if ~complete
                %[css,complete] = LAInteg(css,model,i);
                %if complete
                %    proceed = true;
                %else
                    proceed = false;
                %end
                while ~proceed
                    [cs(:,i),complete] = LAInteg(cs(:,i),model,i,sthresh,true,cind);
                    jfail = false;
                    if ~complete
                        [cs(:,i),complete,jfail] = SIFOInteg(cs(:,i),model,i,sthresh);
                        if jfail
                            iter = iter-1;
                            [cs(:,i),complete] = LAInteg(cs(:,i),model,i,sthresh,false,cind);
                        end
                    end
                    if complete || ~jfail
                        proceed = true;
                    end
                end
            %else
            %    css = ci;
            end
            if iter>3
                complete = true;
            end
            %if iter>4 && ~complete
            %    opt = odeset('Jacobian',@(t,x) jacfn(t,x,model,i));
            %    [~,C] = ode15s(@(t,x) svinteg(t,x,model,i),[0,1e7],css,opt);
            %    css = C(end,:)';
            %    [dc] = svinteg(1,css,model,i);
            %    if max(abs(dc))<1e-6
            %        complete = true;
            %    end
            %end
        end
        [dx,vx] = svinteg(1,cs(:,i),model,i);
        if max(abs(dx))<sthresh || (abs(vx(cind))<0.1)% && max(abs(dx))<1)
            done = true;
            es(:,i) = efrac(cs(:,i),model,i);
            [dcsdp{i},desdp{i}] = spsens(model,k,cs(:,i),es(:,i));
            %cs(:,i) = css(:);
            %es(:,i) = ess(:);
            %dcsdp{i} = Sc;
            %desdp{i} = Se;
            [v(:,i),dvdk{i}] = flxcalc(model,cs(:,i),es(:,i),k,dcsdp{i},desdp{i});
            v(:,i) = v(:,i).*model.d.vpert(:,i);
            dvdk{i} = diag(full(model.d.vpert(:,i)))*dvdk{i};
            %vsc = 100*vx/abs(vx(12));
            %dvscdk = (100/abs(vx(12)))*(dvxdk-vsc*dvxdk(12,:));
            %Css(:,i) = c(end,:)';
        %else
        %    ctr = ctr+1;
        %    c0 = css;
        %    sc = 100/abs(vx(63));
        %    if ctr>2 || max(abs(dx))*sc<0.1
        %        done = true;
        %        ess = efrac(css,model,i);
        %        [Sc,Se] = spsens(model,k,css,ess);
        %        cs(:,i) = css(:);
        %        es(:,i) = ess(:);
        %        dcsdp{i} = Sc;
        %        desdp{i} = Se;
        %        [vx,dvxdk] = flxcalc(model,css,ess,k,Sc,Se);
        %        vx = vx.*model.d.vpert(:,i);
        %        dvxdk = diag(full(model.d.vpert(:,i)))*dvxdk;
        %    end
                
        end
    end
    %vx = vx.*model.d.vpert(:,i);
    %v(:,i) = vx;
    %dvdk{i} = dvxdk;
end
vop = v;
%{
v = v(:,2:end);
dvdk = dvdk(:,2:end);
v1 = v(:);
dvdk = dvdk(:);
dvdk = cell2mat(dvdk);
dvdp = dvdk*dkdp;
%{
dcsdp = dcsdp(:);
dcsdp = cell2mat(dcsdp);
dcsdp = dcsdp*dkdp;
cs = cs(:);
desdp = desdp(:);
desdp = cell2mat(desdp);
desdp = desdp*dkdp;
es = es(:);
%}
idx = model.d.idx;
idz = idx>length(model.rid);
idx = idx(idz);
idx = idx-length(model.rid);
vpred = v1(idx);
vmeas = model.d.flx(idz);
drdk = dvdk(idx,:);
J = [drdp;drdk*dkdp];
std = [std;model.d.err(idz)];
std = std.^2;
W = diag(1./std);
%save('drdp.mat','drdp');
%save('W.mat','W');
%std = 1./std;
r = [res;vpred-vmeas];
%ssres = res'*W*res;
%g = res'*W*drdp;
%g = 2*g';
%}
for i = 1:ncond
    vj = v(:,i);
    vj = vj.*model.d.vpert(:,i);
    dvjdk = dvdk{i};
    dvjdk = diag(model.d.vpert(:,i))*dvjdk*dkdp;
    rmp1 = model.d.rmap{i};
    r1 = rmp1*vj;
    r1 = r1-(model.d.flx{i});
    dr1dp = rmp1*dvjdk;
    res = [res;r1];
    drdp = [drdp;dr1dp];
    std = [std;model.d.err{i}];
end
r = res;
J = drdp;
std = std.^2;
W = diag(1./std);




end


function [efc] = efrac(x,model,ipert)
%x = max(x,0);
%x(x<1e-9) = 0;
c = [1;x(:)];
k = model.p.k;
nr = length(model.p.dEdk);
ncond = 1;
eblocks = model.ssubs.eblocks;
efc = zeros(eblocks(end),1);
%dedc = zeros(eblocks(end),length(x));

%dvdk = zeros(nr,model.p.nk);
for i = 1:nr
    kr = k(model.p.kblocks(i)+1:model.p.kblocks(i+1));
    %dE1dk = model.p.dEdk{i};
    %le = length(dE1dk(:,1));
    %dE1dvup = zeros(le,nupv*ncond);
    %dE1dc = cell(1,ncond-1);
    %dE1dc(1:ncond-1) = zeros(le,nc);
    for j = 1:ncond
        kmult = [model.p.dkdc{i}]'*c(:,j);
        kp = kr.*kmult;
        le = sqrt(length(model.p.dEdk{i}(:,1)));
        E = model.p.E{i} + reshape([model.p.dEdk{i}]*kp,le,le);
        e = E\model.p.B{i};
        e = e*model.d.vpert(i,ipert);
        efc(eblocks(i)+1:eblocks(i+1)) = e;
        %{
        % computing jacobian
        t1 = model.p.dEdk{i};
        t2 = model.p.dkdc{i}(2:end,:);
        t1 = t1*t2';
        %calculating cmult
        t2 = [full(~any(t2,2)),t2];
        k1 = [1;kr];
        cmult = t2*k1;
        dEdc = t1*diag(cmult);
        lc = length(dEdc(1,:));
        dEdcT = reshape(dEdc,le,le*lc);
        E1 = E\dEdcT;
        E1 = mat2cell(E1,le,le*ones(1,lc));
        E1 = E1';
        E1 = cell2mat(E1);
        de = -E1*e;
        de = reshape(de,le,lc);
        dedc(eblocks(i)+1:eblocks(i+1),:) = de;
        %}
    end
end
end


