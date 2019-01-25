function res = compileresult(xopt,model)


% fmin
[r,W,~,v,css] = rescalc(xopt,model);
res.fmin = r'*W*r;

% kinetic parameters and reversibilities
[k,~,rev] = calc_k(model,xopt);

kin = struct('id','','kf',[],'kr',[],'Ki',[],'kact',[]);
nr = length(model.rid);

[res.kinetic_params(1:nr)] = deal(kin);
%x1 = xopt;
%nu = model.ensemble.nu;
%u = x1(1:nu);
%x1(1:nu) = [];
ne = cell2mat(model.ensemble.ne);
%nei = cell2mat(model.ensemble.nei);
%nvr = cell2mat(model.ensemble.nvr);
%nf = ne-1;
%np = nf+nei+nvr;
for i = 1:nr
    res.kinetic_params(i).id = model.rid(i);
    eind = 1:ne(i);
    krevind = 2*eind;
    kfwdind = krevind-1;
    krxn = k(model.p.kblocks(i)+1:model.p.kblocks(i+1));
    kcatal = krxn(1:2*ne(i));
    res.kinetic_params(i).kr = kcatal(krevind);
    res.kinetic_params(i).kf = kcatal(kfwdind);
    krxn(1:2*ne(i)) = [];
    res.kinetic_params(i).Ki = krxn;
end

% Predicted fluxes and concentrations
flx = struct('id','','val',[]);
met = struct('id','','val',[]);
mut = struct('gene_KO','','fluxes',[],'concentrations',[]);
ncond = length(v(1,:));
nmet = length(model.metprop);
[m1(1:ncond)] = deal(mut);
for i = 1:ncond
    %m1 = mut;
    pert = model.d.vpert(:,i);
    if ~any(~pert)
        m1(i).gene_KO = {'WT'};
    else
        m1(i).gene_KO = model.rid(~pert);
    end
    [f(1:nr)] = deal(flx);
    for j = 1:nr
        f(j).id = model.rid(j);
        f(j).val = v(j,i);
    end
    [m(1:nmet)] = deal(met);
    for j = 1:nmet
        m(j).id = model.metprop(j).metid;
        m(j).val = css(j,i);
    end
    m1(i).fluxes = f;
    m1(i).concentrations = m;
end
res.predictions = m1;

% Lack-of-Fit

rever = struct('flxid','','val',[],'data',[],'WRES',[],'SRES',[]);
%flxft = struct('flxid','','val',[],'data',[],'WRES',[],'SRES',[]);

%reversibilities
nrevs = length(model.d.ridx);
[rx(1:nrevs)] = deal(rever);
for i = 1:nrevs
    rx(i).flxid = model.rid(model.d.ridx(i));
    rx(i).val = rev(model.d.ridx(i));
    rx(i).data = model.d.revs(i);
    rx(i).WRES = (rx(i).val-rx(i).data)/model.d.rerr(i);
    rx(i).SRES = (rx(i).WRES).^2;
end
res.residuals.reversibility = rx;

%fluxes
flx = struct('measid','','combination','','mut','','data',[],'val',[],'WRES',[],'SRES',[]);
%nflx = length(model.d.id);
nflx = 0;
for i = 1:ncond
    nflx = nflx+length(model.d.id{i});
end
[flxs(1:nflx)] = deal(flx);
ctr = 0;
%{
idx = model.d.idx;
fdat = model.d.flx;
erdat = model.d.err;
for i = 1:ncond
    relid = idx(idx<=length(model.rid));
    fd1 = fdat(relid);
    err1 = erdat(relid);
    pert = model.d.vpert(:,i);
    for j = 1:length(relid)
        ctr = ctr+1;
        if ~any(~pert)
            flxs(ctr).mut = {'WT'};
        else
            flxs(ctr).mut = model.rid(~pert);
        end
        flxs(ctr).measid = model.rid(relid(j));
        flxs(ctr).combination = model.rid(relid(j));
        flxs(ctr).val = v(relid(j),i);
        flxs(ctr).data = fd1(j);
        flxs(ctr).WRES = (flxs(ctr).val - flxs(ctr).data)/err1(j);
        flxs(ctr).SRES = (flxs(ctr).WRES).^2;
    end
    fdat(idx<=length(model.rid)) = [];
    erdat(idx<=length(model.rid)) = [];
    idx(idx<=length(model.rid)) = [];
    idx = idx-length(model.rid);
end
%}
for i = 1:ncond
    pert = model.d.vpert(:,i);
    for j = 1:length(model.d.id{i})
        ctr = ctr+1;
        if ~any(~pert)
            flxs(ctr).mut = {'WT'};
        else
            flxs(ctr).mut = model.rid(~pert);
        end
        flxs(ctr).measid = model.d.id{i}(j);
        flxs(ctr).combination = model.rid(find(model.d.rmap{i}(j,:)));
        flxs(ctr).val = model.d.rmap{i}(j,:)*v(:,i);
        flxs(ctr).data = model.d.flx{i}(j);
        flxs(ctr).WRES = (flxs(ctr).val - flxs(ctr).data)/model.d.err{i}(j);
        flxs(ctr).SRES = (flxs(ctr).WRES).^2;
    end
end


res.residuals.fluxes = flxs;        


%reinitialization
res.reinit_data = xopt;





end