function [ yf,complete ] = LAInteg( y0,model,ipert,sthresh,first,cind )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%if first
%    h = 1;
%else
    h = 2*(10^-6);
%end
warning('off','all');
etol = 1;
yf = y0;
ymin = 1e-20;
%yf = max(yf,1e-8);

yf(yf<ymin) = 0;
%sthresh = 1e-6;
%yk =y0;
%ykf = y0;
t0 = 0;
%tk = 0;
%tkf = 0;
tf = 0;
%jfail = false;
I0 = eye(length(y0));
[J,dy0] = jacfn(t0,y0,model,ipert);
q = max(abs(dy0)); 
if q <sthresh
    done = true;
    complete = true;
else
    done = false;
    complete = false;
end

   % t = 1;
    %w = 0;
%
dy1 = (I0-h*J)\(h*dy0);
%alf = 1;
%a = (yf+dy1)<=1e-16;
%if any(a)
%       alf = min(max((yf(a)-ymin),0)./abs(dy1(a)));
%       alf = min(alf,1);
%       h = 0.95*alf*h;
%end
    
N0 = I0;
yact = true(size(y0));

%    iter = 0;
    %qbest = q;
    %ybest = y0;
%main loop
%step = 1;
%ctx = 0;
while ~done
    yk = yf;
    tk = tf;
    Nx = N0(:,yact);
    %I = eye(sum(yact));
    I = I0;
    %single step calc
    [J,dyk] = jacfn(tk,yk,model,ipert);
    %dy1 = (w*I-t*J)\(t*dyk);
    %J = Nx'*J*Nx;   %projected gradient
    %dy1 = (I-h*J)\(h*Nx'*dyk);
    Jx = (I-h*J);
    [Qj,Rj] = qr(Jx);
    dy1 = Rj\Qj'*(h*dyk);
    %yf1 = yk+Nx*dy1;
    yf1 = yk+dy1;
    %yf1(yf1<1e-8)=0;
    %two step calc
    %dy1 = (I-(h/2)*J)\((h/2)*Nx'*dyk);
    Jx = (I-(h/2)*J);
    [Qj,Rj] = qr(Jx);
    dy1 = Rj\Qj'*((h/2)*dyk);
    %yfa = yk+Nx*dy1;
    yfa = yk+dy1;
    yfa=max(yfa,0);
    [J,dyk] = jacfn(tk+(h/2),yfa,model,ipert);
    %J = Nx'*J*Nx;
    %dy1 = (I-(h/2)*J)\((h/2)*Nx'*dyk);
    Jx = (I-(h/2)*J);
    [Qj,Rj] = qr(Jx);
    dy1 = Rj\Qj'*((h/2)*dyk);
    %yf2 = yfa+Nx*dy1;
    yf2 = yfa+dy1;
    
    %check error
    err = abs(yf2-yf1)./(max(yf,1)*etol);
    err = rms(err);
    %if q<1
    %    err = 0.9;
    %end
    
    if err > 1  %reject step
        h = h/max(sqrt(err),2);
    elseif isnan(err)   %restart
        h = 1e-6;
        yf = ones(size(yf));
    else
        yf = yf2;
        yf=max(yf,0);
        tf = tk+h;
        [J,dyk,~,flx] = jacfn(tf,yf,model,ipert);
        q = max(abs(dyk));
        qchk = abs(flx(cind));
    %    h1 = 2*h;
        h1 = h/sqrt(err);
        hx = h;
        hmax = 3*h;
        
        if h1<=1.5*h
            h1 = h; %do not rescale step for very small steps
        else
            h1 = min(hmax,h1);
        end
        %if q<1
        %    t = 1;
        %    w = 0;
        %else
        %    t = min(h1,1e9);
        %    w = 1;
        %end
        %dy1 = (w*I-t*J)\(t*dyk);
        dy1 = (I0-h1*J)\(h1*dyk);
        %yf(yf<1e-8) = 0;
        yact = yf<=ymin & dy1<(1e-4*ymin);
        yact = ~yact;
        % check for negative concentrations
        %a = (yf+dy1)<0; %blocking constraints
        %alf = 1;
        %if any(a(yact))
        %       alf = min((yf(a&yact)-ymin)./abs(dyk(a&yact)));
        %       alf = min(alf,1);
               %h1 = alf;
        %end
        h = h1;%*alf;
        if h<=0
            h = 1e-6;
        end
        %if q<1
            %t = 1;
            %w = 0;
        %else
        %    t = min(h,1e9);
        %    w = 1;
        %end
    end
    
    %check for completion
    if first && err < 1 && h > 1e12 && h1/hx < 1.2
        done = true;
        complete = false;
    elseif q<sthresh || (qchk<0.1)% && q<1)
        done = true;
        complete = true;
        %jfail = false;
        %tf = tk;
        %yf = yk;
    elseif err<1 && h>1e15 && h1/hx < 1.2
        done = true;
        %jfail = false;
        complete = false;
    end
    
end
    



end

