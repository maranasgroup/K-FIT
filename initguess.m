function [ x  ] = initguess( model,A,b,xlb,xub )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% computing total number of variables
ne = cell2mat(model.ensemble.ne);
nei = cell2mat(model.ensemble.nei);
nvr = cell2mat(model.ensemble.nvr);
nf = ne-1;
np = nf+nei+nvr;
nvar = sum(np);

%{
%constructing constraint matrices

A = zeros(length(ne),nvar);
b1 = ones(length(ne),1);
b = b1*1e-7;
pos = 0;
inds = zeros(1,nvar);
vrind = inds;
for i = 1:length(ne)
    A(i,pos+1:pos+nf(i)+nei(i)) = 1;
    inds(1,pos+1:pos+nf(i)+nei(i)) = 1;
    vrind(1,pos+nf(i)+nei(i)+1:pos+np(i)) = model.d.flx(i);
    pos = pos+np(i);
end
inds = logical(inds);
A1 = eye(nvar);
A1(inds,:) = [];
vrind(inds) = [];
vrind = vrind';
A = [A;-A;-A1];
b = [b1;-b;vrind];
%}
pos = 0;
inds = zeros(1,nvar);
vrind = inds;
for i = 1:length(ne)
    inds(1,pos+1:pos+nf(i)+nei(i)) = 1;
    vrind(1,pos+nf(i)+nei(i)+1:pos+np(i)) = model.d.flx{1}(i);
    pos = pos+np(i);
end
x = rand(nvar,1);
x(~inds) = 100*x(~inds);
%{
xlb = 1e-7*ones(size(x));
xub = ones(size(xlb));
xub(~inds) = 100000000;
%}
%f = @(x) dgerr(x,model);
%x = fmincon(f,x,A,b,[],[],1e-7*ones(size(x)),ones(size(x)),[],optimset('Display','iter','MaxFunEvals',10000000));
x = fmincon(@(x) 1,x,A,b,[],[],xlb,xub,[],optimset('Display','iter','MaxFunEvals',10000000));
end

function r = dgerr(x,model)

ne = cell2mat(model.ensemble.ne);
nei = cell2mat(model.ensemble.nei);
nrev = cell2mat(model.ensemble.nrev);
nf = ne-1;
np = nf+nei+nrev;

r = zeros(size(ne));
nx = 0;
for i = 1:length(r)
    nx = nx+np(i);
    r(i) = prod(x(nx-nrev(i)+1:nx));
end
r = r-model.dG;
r = r(~model.p.exch)./model.dgsdv(~model.p.exch);
r = r'*r;
end

