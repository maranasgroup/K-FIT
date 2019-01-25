function [ css,vss,complete ] = svsucsubs( c0,model,ipert,sthresh )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
c0 = c0(:);
dKLdk = model.ssubs.dKLdk;
dKRdk = model.ssubs.dKRdk;
dkdE = model.ssubs.dkdE;
k = model.p.k;
%sthresh = 1e-6;
[dc,v] = svinteg(1,c0,model,ipert);
if max(abs(dc))<sthresh
    done = true;
    complete = true;
    css = c0;
    vss = v;
else
    done = false;
end
ck = c0;
iter = 0;
ctr = 0;
q = max(abs(dc));
w = ones(size(ck));
while ~done
    
    q0 = q;
    iter = iter+1;
    if iter>1
        x1k = xk;
        g1k = gk;
    end
    e = efrac(ck,model,ipert);
    kmult = dkdE*e;
    kp = k.*kmult;
    KL = dKLdk*kp;
    KR = dKRdk*kp;
    cki = diag(KL)\KR;%./KL;
    xk = ck;
    gk = cki;
    if iter>1
        dx = (xk-x1k);
        dx(dx==0)=1e-7;
        s = (gk-g1k)./dx;
    else
        s = 100*ones(size(xk));
    %    if max(abs(s)) < 0.8
    %        w = 1./(1-s);
    %    else
    %        w = ones(size(ck));
    %    end
    end
    %protecting against negative concentrations
    %a =(gk-xk)<0;
    %if any(a)
    %    w(a) = 0.99*(1+(gk(a)./(xk(a)-gk(a))));
    %end
    s = s(s>0);
    ck1 = ((1-w).*ck+w.*cki);
    %check for termination:
    [dc,v] = svinteg(1,ck1,model,ipert);
    q = max(abs(dc));
    sx = s>0.995 & s<1.005;
    if (abs(q-q0)<0.01 || abs(q-q0)/q0 < 1e-4) || iter>1000
    %if abs(q-q0)/q0 < 1e-3 && q < 5
    %if 10000*abs(q-q0)<1
        ctr = ctr+1;
    else
        ctr = 0;
    end
    if max(abs(dc))<sthresh
    %if max(abs(ck1-ck))<1e-10 || max(abs(dc))<1e-8
        done = true;
        complete = true;
        css = ck1;
        vss = v;
    %elseif iter>2500 && max(abs(dc))<1
    elseif ctr >5
        done = true;
        complete = false;
        css = ck1;
        vss = v;
    else
        ck = ck1;
    end
end
end

function [efrc] = efrac(x,model,ipert)
%x = max(x,0);
%x(x<1e-9) = 0;
c = [1;x(:)];
k = model.p.k;
nr = length(model.p.dEdk);
ncond = 1;
eblocks = model.ssubs.eblocks;
%efrc = zeros(eblocks(end),1);
efrc = cell(length(eblocks)-1,1);
%dedc = zeros(eblocks(end),length(x));

%dvdk = zeros(nr,model.p.nk);
for i = 1:nr
    kr = k(model.p.kblocks(i)+1:model.p.kblocks(i+1));
    %dE1dk = model.p.dEdk{i};
    %le = length(dE1dk(:,1));
    %dE1dvup = zeros(le,nupv*ncond);
    %dE1dc = cell(1,ncond-1);
    %dE1dc(1:ncond-1) = zeros(le,nc);
    %for j = 1:ncond
        kmult = [model.p.dkdc{i}]'*c(:,1);
        kp = kr.*kmult;
        le = sqrt(length(model.p.dEdk{i}(:,1)));
        E = model.p.E{i} + reshape([model.p.dEdk{i}]*kp,le,le);
        e = E\model.p.B{i};
        e = e*model.d.vpert(i,ipert);
        %efrc(eblocks(i)+1:eblocks(i+1)) = e;
        efrc{i} = e;
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
    %end
end
efrc = cell2mat(efrc); 
end
