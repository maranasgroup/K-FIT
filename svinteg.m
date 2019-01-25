function [ dx,v ] = svinteg( t,x,model,ipert )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%x = max(x,0);
%x(x<1e-12) = 0;
c = [1;x(:)];
k = model.p.k;
nr = length(model.p.dEdk);
ncond = 1;
vkp = zeros(nr,ncond);
warning('off','all');
%dvdk = zeros(nr,model.p.nk);
for i = 1:nr
    kr = k(model.p.kblocks(i)+1:model.p.kblocks(i+1));
    vdef = model.p.vdef{i};
    %dE1dk = model.p.dEdk{i};
    %le = length(dE1dk(:,1));
    %dE1dvup = zeros(le,nupv*ncond);
    %dE1dc = cell(1,ncond-1);
    %dE1dc(1:ncond-1) = zeros(le,nc);
    for j = 1:ncond
        kmult = [model.p.dkdc{i}]'*c(:,j);
        kp = kr.*kmult;
        le = sqrt(length(model.p.dEdk{i}(:,1)));
        E = model.p.E{i} + [reshape([model.p.dEdk{i}]*kp,le,le)];
        [Q1,R1] = qr(E);
        e = R1\(Q1'*model.p.B{i});
        %e = E\model.p.B{i};
        vkp(i,j) = kp(vdef(1,1))*e(vdef(1,2));
        if length(vdef(:,1))>1
            vkp(i,j) = vkp(i,j) - (kp(vdef(2,1))*e(vdef(2,2)));
        end
        %{
        % computing gradients
        t1 = dEdk{i}*diag(kmult);
        lk = length(kp);
        dEdkT = reshape(t1,le,le*lk);
        E1 = E\dEdkT;
        E1 = mat2cell(E1,le,le*ones(1,lk));
        E1 = E1';
        E1 = cell2mat(E1);
        dedk = E1*e;
        dedk = reshape(dedk,le,lk);
        z = zeros(1,length(kp));
        z1 = z;
        z1(vdef(1,1)) = 1;
        dvkp = z1*e(vdef(1,2)) + (dedk(vdef(1,2),:)*kp(vdef(1,1)));
        if length(vdef(:,1))>1
            z1 = z;
            z1(vdef(2,1)) = 1;
            dvkp = dvkp-(z1*e(vdef(2,2)) + (dedk(vdef(2,2),:)*kp(vdef(2,1))));
        end
        dvdk(i,model.p.kblocks(i)+1:model.p.kblocks(i+1)) = dvkp;
        %}
    end
end
%J = model.p.S*dvdk;
v = vkp;
% scaling integration to improve performance

vin = max(model.p.S,0)*max(v,0) + min(model.p.S,0)*min(v,0);
vout = min(model.p.S,0)*max(v,0) + max(model.p.S,0)*min(v,0);
vout = -vout;
xsc = max(vin,vout);
xsc = max(xsc,1);

%W = eye(length(v));
W = diag(model.d.vpert(:,ipert));

dx = model.p.S*W*v;
v = W*v;
t1 = max(abs(dx));
t2 = min(x);
%sc = 1000/abs(v(12));
%if max(abs(dx))<(1e-3)/sc
%    dx = 0*dx;
%end

d1 = diag(xsc);
%if max(abs(dx))>1
%    dx = d1\dx;
%end
end

