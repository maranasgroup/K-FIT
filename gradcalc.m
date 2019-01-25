function [model] = gradcalc(model)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

nr = length(model.rid);
nc = model.p.nc;
nenz = cell2mat(model.ensemble.ne);

% mapping kinetic parameters to concentration fold changes:
dkdc = cell2mat([model.p.dkdc]');
dkdc = dkdc(2:end,:);
dkdc = dkdc';

% mapping kinetic parameters to enzyme fractions:
eblocks = model.ssubs.eblocks;
kblocks = model.p.kblocks;
nk = kblocks(end);
ne = eblocks(end);
dkde = sparse(nk,ne);
for i = 1:nr
    S = [-model.kinetic(i).S,model.kinetic(i).I]';
    S = max(S,0);
    dkde(kblocks(i)+1:kblocks(i+1),eblocks(i)+1:eblocks(i+1)) = S;
end

%setting up ensemble equations:
B = sparse((nc+ne),nk);
%enzyme fraction balances:
ctr = 0;
for i = 1:nr
    inds = 1:(nenz(i)^2);
    inds = reshape(inds,nenz(i),nenz(i));
    inds = inds';
    inds = inds(:);
    dEdk = model.p.dEdk{i}(inds,:);
    for j = 1:nenz(i)
        B(ctr+j,kblocks(i)+1:kblocks(i+1)) = sum(dEdk((j-1)*nenz(i)+1:(j-1)*nenz(i)+nenz(i),:),1);
    end
    ctr = ctr+nenz(i);
end

%substrate balances:
B(ctr+1:end,:) = model.ssubs.dKRdk-model.ssubs.dKLdk;

% dcdk part of the derivative equation

dCLdk = sparse(nc*(nc+ne),nk);
for i = 1:(nc+ne)
    dCLdk((i-1)*nc+1:i*nc,:) = (diag(B(i,:))*dkdc)';
end
ind = 1:(nc*(ne+nc));
ind = reshape(ind,nc,nc+ne);
ind = ind';
ind = ind(:);
dCLdk = dCLdk(ind,:);

%dedk part of the derivative equation

dELdk = sparse(ne*(nc+ne),nk);
for i = 1:(nc+ne)
    dELdk((i-1)*ne+1:i*ne,:) = (diag(B(i,:))*dkde)';
end
ind = 1:(ne*(ne+nc));
ind = reshape(ind,ne,nc+ne);
ind = ind';
ind = ind(:);
dELdk = dELdk(ind,:);

%accounting for sum(enzyme fractions) = 1
A = sparse(ne+nc,ne);
ctr = 0;
for i = 1:nr
    A(ctr+1,eblocks(i)+1:eblocks(i+1))=1;
    ctr = ctr+nenz(i);
end

%RHS of derivatice equation
%note that k here is a product of c and e that uniquely maps to the corresponding k
%in dkdc and dkde
dRdk = sparse(nk*(ne+nc),nk);
ctr = 0;
for i = 1:(ne+nc)
    dRdk(ctr+1:ctr+nk,:) = -diag(B(i,:));
    ctr = ctr+nk;
end
ind = 1:(nk*(ne+nc));
ind = reshape(ind,nk,nc+ne);
ind = ind';
ind = ind(:);
dRdk = dRdk(ind,:);

%collect required matrices
model.grads.A = A;
model.grads.B = B;
model.grads.dCLdk = dCLdk;
model.grads.dELdk = dELdk;
model.grads.dRdk = dRdk;
model.grads.dkde = dkde;
model.grads.dkdc = dkdc;

end

