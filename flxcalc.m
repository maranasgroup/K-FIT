function [ v,dvdk ] = flxcalc( model,css,ess,k,Sc,Se )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

nk = length(k);
nr = length(model.p.dEdk);
v = zeros(nr,1);
dvdk = zeros(nr,nk);
%ncond = 1;
kblocks = model.p.kblocks;
eblocks = model.ssubs.eblocks;
kemult = model.grads.dkde*ess;
dkdc = [double(~any(model.grads.dkdc,2)),model.grads.dkdc];
kcmult = dkdc*([1;css]);
kc = k.*kcmult;
ke = k.*kemult;
ce = kemult.*kcmult;
%main loop
dk = eye(nk);
for i = 1:nr
    %kr = k(kblocks(i)+1:kblocks(i+1));
    vdef = model.p.vdef{i};
    kcr = kc(kblocks(i)+1:kblocks(i+1));
    ker = ke(kblocks(i)+1:kblocks(i+1));
    cer = ce(kblocks(i)+1:kblocks(i+1));
    dkr = dk(kblocks(i)+1:kblocks(i+1),:);
    dc = model.p.dkdc{i};
    dc = dc(2:end,:);
    dc = dc';
    dc = dc*Sc;
    er = ess(eblocks(i)+1:eblocks(i+1));
    de = Se(eblocks(i)+1:eblocks(i+1),:);
    j = vdef(1,1);
    v(i) = kcr(vdef(1,1))*er(vdef(1,2));
    dvdk(i,:) = cer(j)*dkr(j,:) + ker(j)*dc(j,:) + kcr(j)*de(vdef(1,2),:);
    if length(vdef(:,1))>1
        v(i) = v(i)-(kcr(vdef(2,1))*er(vdef(2,2)));
        j = vdef(2,1);
        dvdk(i,:) = dvdk(i,:)- (cer(j)*dkr(j,:) + ker(j)*dc(j,:) + kcr(j)*de(vdef(2,2),:));
    end
end


end

