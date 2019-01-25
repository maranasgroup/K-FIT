function [ model ] = ccalcsetup( model )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


nc = model.p.nc;
nk = model.p.nk;
nr = length(model.rid);
nex = cell2mat(model.ensemble.ne);
ne = sum(nex);
eblocks = zeros(1,nr+1);

for i = 1:nr
    eblocks(i+1) = eblocks(i)+nex(i);
end

%setting up a platform to map kinetic parameters to enzyme fractions

S = zeros(ne,nk);
for i = 1:nr
    S1 = model.kinetic(i).S;
    nkr = length(S1(1,:));
    S(eblocks(i)+1:eblocks(i+1),model.p.kblocks(i)+1:model.p.kblocks(i)+nkr) = S1;
end

dkdE = -min(S,0);
dkdE = dkdE';
S = S';

dcdk = cell2mat([model.p.dkdc]');
dcdk(1,:) = [];
l2 = ~any(S,2);
dcdk(:,l2) = 0;
dKRdk = sparse(nc,nk);
for i = 1:nc
    f1 = find(dcdk(i,:));
    f1 = [f1-1;f1;f1+1];
    f1 = f1(:);
    f1(f1<1) = [];
    f1(f1>nk) = [];
    S1 = S;
    l1 = false(nk,1);
    l1(f1) = true;
    S1(~l1,:) = 0;
    kx = ismember(S1,-S1(logical(dcdk(i,:)),:),'rows');
    dKRdk(i,:) = double(kx');
end

model.ssubs.dKLdk = dcdk;
model.ssubs.dKRdk = dKRdk;
model.ssubs.eblocks = eblocks;
model.ssubs.dkdE = sparse(dkdE);



end

