function [ Ain,bin,xlb,xub ] = constraints( model )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% computing total number of variables
ne = cell2mat(model.ensemble.ne);
nei = cell2mat(model.ensemble.nei);
nvr = cell2mat(model.ensemble.nvr);
nf = ne-1;
np = nf+nei+nvr;
nvar = sum(np);

A = zeros(length(ne),nvar);
b1 = ones(length(ne),1);
b = b1*1e-3;
b1 = b1-1e-3;
pos = 0;
inds = zeros(1,nvar);
indsnei = zeros(1,nvar);
vrind = inds;
for i = 1:length(ne)
    A(i,pos+1:pos+nf(i)+nei(i)) = 1;
    inds(1,pos+1:pos+nf(i)+nei(i)) = 1;
    indsnei(1,pos+1:pos+nf(i)+1:pos+1:pos+nf(i)+nei(i)) = 1;
    if model.d.flx{1}(i) == 0
        f = 1;
    else
        f = (model.d.flx{1}(i)/abs(model.d.flx{1}(i)));
    end
    vrind(1,pos+nf(i)+nei(i)+1:pos+np(i)) = model.d.flx{1}(i)+(f*1e-7);
    pos = pos+np(i);
end
inds = logical(inds);
indsnei = logical(indsnei);
A1 = eye(nvar);
A1(inds,:) = [];
vrind(inds) = [];
vrind = vrind';
A = [A;-A;-A1];
b = [b1;-b;vrind];

x = rand(nvar,1);
x(~inds) = 100*x(~inds);
xlb = 1e-7*ones(size(x));
xlb(inds) = 1e-3;
xlb(indsnei) = 1e-7;
xub = ones(size(xlb));
xub(~inds) = 1e4;

Ain = A;
bin = b;
end

