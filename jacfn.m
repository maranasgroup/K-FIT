function [ J,dx,dvdk,v ] = jacfn( t,x,model,ipert )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%x(x<1e-8) = 0;
warning('off','all');
c = [1;x(:)];
k = model.p.k;
nr = length(model.p.dEdk);
ncond = 1;
vkp = zeros(nr,ncond);
dvdk = zeros(nr,length(x));
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
        [Qz,Rz] = qr(E);
        e = Rz\(Qz'*model.p.B{i});
        %e = E\model.p.B{i};
        vkp(i,j) = kp(vdef(1,1))*e(vdef(1,2));
        if length(vdef(:,1))>1
            vkp(i,j) = vkp(i,j) - (kp(vdef(2,1))*e(vdef(2,2)));
        end
        %
        % computing jacobian
        t1 = model.p.dEdk{i};
        t1 = t1*diag(kr);
        t2 = model.p.dkdc{i}(2:end,:);
        t1 = t1*t2';
        %calculating cmult
        t2 = [full(~any(t2,2)),t2];
        %k1 = [1;kr];
        %cmult = t2*k1;
        dEdc = t1;
        lc = length(dEdc(1,:));
        dEdcT = reshape(dEdc,le,le*lc);
        E1 = Rz\(Qz'*dEdcT);
        %E1 = E\dEdcT;
        E1 = mat2cell(E1,le,le*ones(1,lc));
        E1 = E1';
        E1 = cell2mat(E1);
        dedc = -E1*e;
        dedc = reshape(dedc,le,lc);
        t2 = t2(:,2:end)*diag(kr);
        dvkp = [t2(:,vdef(1,1))]'*e(vdef(1,2)) + (dedc(vdef(1,2),:)*kp(vdef(1,1)));
        if length(vdef(:,1))>1
            dvkp = dvkp - ([t2(:,vdef(2,1))]'*e(vdef(2,2)) + (dedc(vdef(2,2),:)*kp(vdef(2,1))));
        end
        dvdk(i,:) = dvkp;        
               
            
        
        
        %t1 = model.p.dEdk{i}*diag(kmult);
        %lk = length(kp);
        %dEdkT = reshape(t1,le,le*lk);
        %E1 = E\dEdkT;
        %E1 = mat2cell(E1,le,le*ones(1,lk));
        %E1 = E1';
        %E1 = cell2mat(E1);
        %dedk = E1*e;
        %dedk = reshape(dedk,le,lk);
        %z = zeros(1,length(kp));
        %z1 = z;
        %z1(vdef(1,1)) = 1;
        %dvkp = z1*e(vdef(1,2)) + (dedk(vdef(1,2),:)*kp(vdef(1,1)));
        %if length(vdef(:,1))>1
        %    z1 = z;
        %    z1(vdef(2,1)) = 1;
        %    dvkp = dvkp-(z1*e(vdef(2,2)) + (dedk(vdef(2,2),:)*kp(vdef(2,1))));
        %end
        %dvdk(i,model.p.kblocks(i)+1:model.p.kblocks(i+1)) = dvkp;
        %}
    end
end
%W = eye(length(v));
v = vkp;

vin = max(model.p.S,0)*max(v,0) + min(model.p.S,0)*min(v,0);
vout = min(model.p.S,0)*max(v,0) + max(model.p.S,0)*min(v,0);
vout = -vout;
xsc = max(vin,vout);
xsc = max(xsc,1);
d1 = diag(xsc);

W = diag(model.d.vpert(:,ipert));
dx = model.p.S*W*v;

J = model.p.S*W*dvdk;
%J = d1\J;

end

