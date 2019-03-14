function [r_name,km_name,km,dkm] = eval_mm_params(model,res)

%initialize
r_name = cell(0,1);
km_name = cell(0,1);
km = zeros(0,1);
dkm = zeros(0,1);

%create list of metabolites
mets = [model.metprop.metid]';

%loop over all reactions
for i = 1:length(model.rid)
    %assemble kinetic parameters
    k = zeros(0,1);
    dk = k;
    %assembling catalytic k
    for j = 1:length(res.kinetic_params(i).kf)
        k = [k;res.kinetic_params(i).kf(j);res.kinetic_params(i).kr(j)];
        dk = [dk;res.kinetic_params(i).kf_sd(j);res.kinetic_params(i).kr_sd(j)];
    end
    %assembling regulatory k
    k = [k;res.kinetic_params(i).Ki];
    dk = [dk;res.kinetic_params(i).Ki_sd];
    
    % Defining symbolic variables
    nk = length(k);
    ksym = sym('k',[nk,1]);
    
    %determine symbolic expression for the denominator of rate law-
    %contains all Km terms
    nr = sqrt(length(model.p.dEdk{i}(:,1)));
    A = model.p.E{i} + (reshape(model.p.dEdk{i}*ksym,nr,nr));
    B = model.p.B{i};
    e = A\B;
    vdef = model.p.vdef{i};
    V = (ksym(vdef(1,1))*e(vdef(1,2)));
    if size(vdef,1)>1
        for q = 2:size(vdef,1)
            V = V - (ksym(vdef(q,1))*e(vdef(q,2)));
        end
    end
    %e = e(end);
    [num,den] = numden(V);
    num = simplify(num);
    den = simplify(den);
    
    % identifying kinetic parameters linked to metabolite concentraitons
    vars = ~logical(full(model.p.dkdc{i}(1,:)))';
    vars = ksym(vars);
    
    % Separating denominator into component terms
    [cx,tx] = coeffs(den,vars);
    cterm = cx(ismember(tx,1));
    cx = cx(~ismember(tx,1));
    tx = tx(~ismember(tx,1));
    %J = jacobian(tx,vars);
    
    % map ks to metabolite names
    dkdcx = model.p.dkdc{i}(2:end,:);
    mlist = cell(length(vars),1);
    for j = 1:length(vars)
        pos = find(ismember(ksym,vars(j)));
        pos = find(dkdcx(:,pos));
        mlist(j) = mets(pos);
    end
    
    %assembling Kms
    for j = 1:length(tx)
        Km = (cx(j)*tx(j))/cterm;
        %kmctr = kmctr+1;
        % identifying met combinations in Kms
        done = false;
        kname = 'Km';
        Jac = tx(j);
        while (~done)
            Jac = jacobian(Jac,vars);
            pos = find(~ismember(Jac,0));
            if ~isempty(pos)
                pos = pos(1);
                kname = [kname,'_',mlist{pos}];
                Jac = Jac(pos);
            else
                done = true;
            end
        end
        km_name = [km_name;{kname}];
        Kmval = double(subs(Km,ksym,k));
        dKmval = double(subs(Km,ksym,(k+dk)));
        dKm = abs(dKmval-Kmval);
        km = [km;Kmval];
        dkm = [dkm;dKm];
        r_name = [r_name;model.rid(i)];
    end
    
    %assembling Vmax
    
    [cx,tx] = coeffs(num,vars);
    for j = 1:length(tx)
        vmax = (cx(j)*tx(j))/cterm;
        vval = double(subs(vmax,ksym,k));
        if vval<0
            vmax = -vmax;
            vnm = 'Vmax_rev';
        elseif vval>0
            vnm = 'Vmax_fwd';
        else
            vnm = 'Vmax';
        end
        Jac = tx(j);
        done = false;
        while (~done)
            Jac = jacobian(Jac,vars);
            pos = find(~ismember(Jac,0));
            if ~isempty(pos)
                pos = pos(1);
                vnm = [vnm,'_',mlist{pos}];
                Jac = Jac(pos);
            else
                done = true;
            end
        end
        km_name = [km_name;{vnm}];
        Vmval = double(subs(vmax,ksym,k));
        km = [km;Vmval];
        dvmval = double(subs(vmax,ksym,(k+dk)));
        dvm = abs(dvmval-Vmval);
        dkm = [dkm;dvm];
        r_name = [r_name;model.rid(i)];
    end
end





end