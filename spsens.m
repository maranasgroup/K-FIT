function [dcdk,dedk] = spsens(model,k,c,e)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%nr = length(model.rid);
ne = length(e);
nc = length(c);
nk = length(k);
kemult = model.grads.dkde*e;
dkdc = [double(~any(model.grads.dkdc,2)),model.grads.dkdc];
kcmult = dkdc*([1;c]);
kc = k.*kcmult;
ke = k.*kemult;
ce = kemult.*kcmult;
CL = model.grads.dCLdk*ke;
CL = reshape(CL,ne+nc,nc);
EL = model.grads.dELdk*kc;
EL = reshape(EL,ne+nc,ne);
EL = EL+model.grads.A;
R = model.grads.dRdk*ce;
R = reshape(R,ne+nc,nk);
L = [EL,CL];

% try QR decomposition
[Q,Rx] = qr(L);
S = Rx\(Q'*R);

if any(any(isnan(S)))
    %try LDL solver
    [Lw,Dg] = ldl((L+L')/2);
    Y = Lw\R;
    Z = Dg\Y;
    S = (Lw')\Z;
end

dedk = S(1:ne,:);
dcdk = S(ne+1:ne+nc,:);

end

