function [ ybest,complete,jfail,tf ] = SIFOInteg( y0,model,ipert,sthresh )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

h = 2*(10^-6);
warning('off','all');
%sthresh = 1e-6;
ymin = 1e-20;
%etol = 10000;
yf = max(y0,ymin);
ymax = 1e50;
%yf(yf<1e-8) = 0;
%yk =y0;
%ykf = y0;
t0 = 0;
%tk = 0;
%tkf = 0;
tf = 0;
I0 = eye(length(y0));
[~,dy0] = jacfn(t0,y0,model,ipert);
q = max(abs(dy0)); 
if q <sthresh
    done = true;
    complete = true;
else
    done = false;
end
jfail = false;
   % t = 1;
    %w = 0;
%{
dy1 = (-J)\(dy0);
alf = 1;
a = (yf+dy1)<0;
if any(a)
       alf = min(yf(a)./abs(dy1(a)));
       alf = min(alf,1);
       %h = alf*h;
end
  %}  
N0 = I0;
yact = true(size(y0));

    iter = 0;
    qbest = q;
    ybest = y0;
%main loop
step = 1;
ctx = 0;
while ~done
    ctx = ctx+1;
    yk = yf;
    tk = tf;
    Nx = N0(:,yact);
    %I = Nx'*I0*Nx;
    %single step calc
    [J,dyk] = jacfn(tk,yk,model,ipert);
    %dy1 = (w*I-t*J)\(t*dyk);
    J = Nx'*J*Nx;   %projected gradient
    dy1 = (-J'*J)\(J'*Nx'*dyk);
    dy1 = Nx*dy1;
    a1 = (yf+dy1)<0; %blocking constraints
    a2 = (yf+dy1)>1e50;
    sldef = 1;
    sudef = 1;
    if any(a1)
        sldef = min((ymin-yk(a1))./dyk(a1));
    end
    if any(a2)
        sudef = min((ymax-yk(a2))./dyk(a2));
    end
    alf = 0.7*min([sldef,sudef,1]);
    %{
    a = (a1&yact)|(a2&yact);
    alf = 1;
    if any(a)
        sl = yf-ymin;
        %sl = yf;
        su = 1e50-yf;
        dl = [dy1;dy1];
        ax = [a1&yact;a2&yact];
        s = [sl;su];
        alf = min(s(ax)./abs(dl(ax)));
           alf = min(alf,1);
           %h1 = alf;
    end
    ax = yf>ymin;
    if max(abs(dy1(ax))./(abs(yf(ax))))<1e-20 
        done = true;
        complete = false;
        jfail = false;
        ybest = yf;
    end
        %}
    %stepf = false;
    %while ~stepf
        yf1 = yk+step*alf*dy1;
        %yf1(yf1<1e-8) = 0;
    %    [~,dyk] = jacfn(tf,yf1,model,ipert);
    %    qt = max(abs(dyk));
    %    if qt<=q+1
    %        stepf = true;
    %    else
    %        step = step*0.5;
    %    end
    %end
    %yf1(yf1<1e-8)=0;
    %{
    %two step calc
    dy1 = (I-(h/2)*J)\((h/2)*Nx'*dyk);
    yfa = yk+Nx*dy1;
    [J,dyk] = jacfn(tk+(h/2),yfa,model,ipert);
    J = Nx'*J*Nx;
    dy1 = (I-(h/2)*J)\((h/2)*Nx'*dyk);
    yf2 = yfa+Nx*dy1;
    %}
    %check error
    %err = abs(yf2-yf1)/etol;
    %err = rms(err);
    %if q<1
    %    err = 0.9;
    %end
    
    %if err > 1  %reject step
    %    h = h/sqrt(err);
    %else
        yf = yf1;
        %yf(yf<1e-8) = 0;
        tf = tk+h;
        [J,dyk] = jacfn(tf,yf,model,ipert);
        q = max(abs(dyk));
        if q<qbest
            qbest = q;
            ybest = yf;
            iter = 0;
        else
            iter = iter+1;
        end
        if det(J'*J)==0
        %if rank(J) < length(yf)
            jfail = true;
            done = true;
            complete = false;
        end
        if iter>=25
            jfail = true;
            complete = false;
            done = true;
        end
    %    h1 = 2*h;
        %h1 = h/sqrt(err);
        %hmax = 2*h;
        %h1 = min(hmax,h1);
        %if h1<=1.5*h
        %    h1 = h; %do not rescale step for very small steps
        %end
        %if q<1
        %    t = 1;
        %    w = 0;
        %else
        %    t = min(h1,1e9);
        %    w = 1;
        %end
        %dy1 = (w*I-t*J)\(t*dyk);
        %{
        if ~jfail
            dy1 = (-J'*J)\(J'*dyk);
            yact = (yf<=ymin & dy1<1e-4*ymin) | (yf==1e50 & dy1>=0);
            actf = false;
            while ~actf
                Nz = N0(:,~yact);
                Jx = Nz'*J*Nz;
                dy2 = (-Jx'*Jx)\(Jx'*Nz'*dyk);
                dy2 = Nz*dy2;
                yact2 = (yf<=ymin & dy2<1e-4*ymin) | (yf==1e50 & dy2>=0);
                if isequal(yact,yact2)
                    actf = true;
                else
                    ac1 = find(yact);
                    ac2 = find(yact2);
                    yact = yact2;
                end
            end

            yact = ~yact;
        end
        %}
        % check for negative concentrations
        %{
        a1 = (yf+dy1)<0; %blocking constraints
        a2 = (yf+dy1)>1e12;
        alf = 1;
        if any(a1(yact)|a2(yact))
               alf = min(yf((a1&yact)|(a2&yact))./abs(dy1((a1&yact)|(a2&yact))));
               alf = min(alf,1);
               %h1 = alf;
        end
        %}
        %h = h1*alf;
        %if q<1
            %t = 1;
            %w = 0;
        %else
        %    t = min(h,1e9);
        %    w = 1;
        %end
    %end
    
    %check for completion
    if q<sthresh
        done = true;
        complete = true;
        jfail = false;
        %tf = tk;
        %yf = yk;
    elseif ctx>200
        done = true;
        jfail = false;
        complete = false;
    end
    
end
    



end

