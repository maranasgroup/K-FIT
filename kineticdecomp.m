function [ model ] = kineticdecomp( model,kinetic_file )
% KINETICDECOMP parses the kinetic properties of reactions detailed in
% KINETIC_FILE and performs an elementary step decomposition for the
% corresponding reactions in MODEL. 
% MODEL is updated with the kinetic field for storing kinetic parameters
% and elementary step decompositions
%KINETIC_FILE contains the following fields:
%       Reaction ID corresponding to IDs in the model
%       mechanism: 'seq', 'ppg', 'act', 'pass', 'diff'
%       Substrate binding order: Substrate names separated by ";"
%       Product release order: Product names separated by ";"
%       Competitive inhibitors
%       Uncompetitive inhibitors
%       Non-competitive inhibitors
%       Activators
%       Flux value
%       Flux standard deviation


% Initializing structure fields
rx = struct('id','','mech','','subs',{''},'pdts',{''},'c_in',{''},'uc_in',[],'nc_in',[],'act',[],'cs',[],'cp',[],'vdef',[],'dEdk',[],'I',[],'S',[],'dkdc',[],'nk',[],'exch',[],'sb',[]);
%mx = struct('metid','','x',[],'dkcdx',[],'dkpdx',[],'vmeas',[],'signv',[],'drcdx',[],'drpdx',[],'kmap',[]);

% read kinetic file
[~,txt,raw] = xlsread(kinetic_file);
entries = raw(1,:);
ent = [{'ID'};{'mechanism'};{'SBO'};{'PRO'};{'CI'};{'UCI'};{'NCI'};{'act'};{'exch'};{'sub'}];
if length(entries) < 6 || sum(ismember(ent,entries)) < 10
    error('Insufficient data to perform kinetic decomposition. Please complete mechanism file')
end
txt = txt(2:end,:);
raw = raw(2:end,:);
rs = txt(:,ismember(entries,{'ID'}));
ms = txt(:,ismember(entries,{'mechanism'}));
sb = txt(:,ismember(entries,{'SBO'}));
pr = txt(:,ismember(entries,{'PRO'}));
ci = txt(:,ismember(entries,{'CI'}));
u_ci = txt(:,ismember(entries,{'UCI'}));
n_ci = txt(:,ismember(entries,{'NCI'}));
act = txt(:,ismember(entries,{'act'}));
model.param = false(length(model.rid),1);
model.param(ismember(model.rid,rs)) = true;
nr = length(rs);
nm = length([model.metprop.metid]');
exc = cell2mat(raw(:,ismember(entries,{'exch'})));
sub = cell2mat(raw(:,ismember(entries,{'sub'})));
mets = [model.metprop.metid]';
cx = zeros(nm,1);   % cx will be used to adjust kinetic parameters in 
                    % mutant cases to account for metabolite concentration 
                    % fold changes

% Initilize kinetic structures
[model.kinetic(1:nr,1)] = deal(rx);
%[model.kinetic.mets(1:nm,1)] = deal(mx);

% Null-space matrix to handle free fluxes
S = model.S;
[R,jb] = rref(S);
vfind = 1:length(S(1,:));
vfind(jb) = [];
Nfree = R(1:length(jb),vfind);
N = zeros(length(S(1,:)),length(vfind));
N(jb,:) = -Nfree;
N(vfind,:) = eye(length(vfind));


% Elementary step decomposition
for i = 1:nr
    r = rx;
    r.id = rs{i};
    r.mech = ms{i};
    subs = sb{i};
    pdt = pr{i};
    if ~isempty(subs)
    subs = regexp(subs,';','Split');
    else
        subs = '';
    end
    if ~isempty(pdt)
    pdt = regexp(pdt,';','Split');
    else
        pdt = '';
    end
    if ~isempty(ci{i})
    in_c = regexp(ci{i},';','Split');
    else
        in_c = '';
    end
    if ~isempty(u_ci{i})
    in_uc = regexp(u_ci{i},';','Split');
    else
        in_uc = '';
    end
    if ~isempty(n_ci{i})
    in_nc = regexp(n_ci{i},';','Split');
    else
        in_nc = '';
    end
    if ~isempty(act{i})
    actv = regexp(act{i},';','Split');
    else
        actv = '';
    end
    r.subs = subs;
    r.pdts = pdt;
    r.c_in = in_c;
    r.uc_in = in_uc;
    r.nc_in = in_nc;
    r.act = actv;
    ns = length(subs);
    np = length(pdt);
    
     = length(actv);
    nci = length(in_c);
    nuci = length(in_uc);
    nnci = length(in_nc);
    r.cs = zeros(ns,1);
    r.cp = zeros(np,1);
    rind = ismember(model.rid,{r.id});
    rev = model.rxnprop(rind).rev;
    switch r.mech
        case 'seq'
            %construction of elementary S-matrix
            nk = 2*(ns+np+1);
            ne = ns+np+1;
            Se = sparse(ne+1,nk);
            for j = 1:1:ne
                Se(j,2*j-1) = -1;
                Se(j+1,2*j-1) = 1;
                Se(:,2*j) = -Se(:,2*j-1);
            end
            Se(1,:) = Se(1,:) + Se(end,:);
            Se(end,:) = [];
            dkdc = sparse(nm,nk);
            for j = 1:1:ns
                c = cx;
                c(ismember(mets,subs{j})) = 1;
                dkdc(:,(2*j-1)) = c;
            end
            for j = 1:1:np
                c = cx;
                c(ismember(mets,pdt{j})) = 1;
                dkdc(:,(2*ns+2)+(2*j)) = c;
            end
            if ~rev
                j = 2*ns+2:2:nk;
                Se(:,j) = [];
                dkdc(:,j) = [];
                r.vdef = [2*ns+np+1,ne];    % definition of steady-state flux in terms of kinetics
            else
                r.vdef = [nk-1,ne;nk,1];
            end
            
            % Inhibitors
            cI = sparse(ne,nci+nnci);
            dkcidc = sparse(nm,nci+nnci);
            ncI = sparse(ne,nnci);
            dkncidc = sparse(nm,nnci);
            ucI = sparse(ne,nuci+nnci);
            dkucidc = sparse(nm,nuci+nnci);
            if nci > 0
                cI(1,1:nci) = 1;
                for j = 1:nci
                    c = cx;
                    c(ismember(mets,in_c{j})) = 1;
                    dkcidc(:,j) = c;
                end
            end
            if nnci > 0
                ncI(1:(ns+1),:) = 1;
                cI(1,nci+1:nci+nnci) = 1;
                ucI(2,nuci+1:nuci+nnci) = 1;
                for j = 1:nnci
                    c = cx;
                    c(ismember(mets,in_nc{j})) = 1;
                    dkncidc(:,j) = c;
                    dkcidc(:,nci+j) = c;
                    dkucidc(:,nuci+j) = c;
                end
            end
            if nuci > 0
                %ucI(ns+1,:) = 1;
                ucI(2,1:nuci) = 1;
                for j = 1:nuci
                    c = cx;
                    c(ismember(mets,in_uc{j})) = 1;
                    dkucidc(:,j) = c;
                end
            end
            %Ireg = [cI,ncI,ucI];
            Ireg = [cI,ucI];
            %dkidc = [dkcidc,dkncidc,dkucidc];
            dkidc = [dkcidc,dkucidc];
            ni = size(Ireg,2);
            
            % Activators
            if nact > 0
                Se(1,end+1) = 1;
                Se(end+1,end) = -1;
                Se(:,end+1) = Se(:,end); % should this not be negative of Se(:,end)?
                ne = ne+1;
                nk = nk+2*nact;
                dkadc = sparse(nm,2*nact);
                for j = 1:1:nact
                    dkadc(ismember(mets,actv{j}),2*j-1) = 1;
                end
                dkdc = [dkdc,dkadc];
            end
            
            %setting up dEdk
            % first row:
            dEdk = [sparse(ne,nk),Ireg];
            
            %including inhibitor Ks in the other terms
            e2 = sparse(ne*(ne-1),ni);
            
            %other elements
            e1 = sparse(ne*(ne-1),nk);
            S1 = min(Se,0);
            S2 = max(Se,0);
            for j = 1:ne-1
                e1((j-1)*ne+1:j*ne,:) = abs(Se*diag(S2(j+1,:)));
                e1((j-1)*ne+j+1,:) = S1(j+1,:);
            end
            dEdk = [dEdk;e1,e2];
            r.dEdk = dEdk;
            r.S = Se;
            r.I = Ireg;
            r.dkdc = [dkdc,dkidc];
            r.nk = nk+ni;

            %{
            nstep = ns+np+1;
            % Definition of Enzyme-Metabolite complexes
            r.E = cell(nstep,1);
            r.E{1} = 'E';
            for j = 2:ns+1
                r.E{j} = [r.E{j-1},'-S',num2str(j-1)];
            end
            r.E{end} = ['E-P',num2str(np)];
            for j = nstep-1:-1:2+ns
                r.E{j} = [r.E{j+1},'-P',num2str(nstep-j+1)];
            end
            r.ce = zeros(size(r.E));
            if rev
                nk = (2*ns)+np+1;
                r.k = zeros(nk,1);
                r.ks = ones(nk,1);
                r.ksb = [r.k,1000*r.ks];
                r.ksd = [r.k,r.k];
                for j = 1:ns
                    r.ksd(2*j-1,:) = [2*j-1,j];
                    r.ksd(2*j,:) = [2*j,j+1];
                end
                j = j+1;
                for j1 = 2*ns+1:nk
                    r.ksd(j1,:) = [j1,j];
                    j = j+1;
                end
                %}

        %case 'ppg'
            
        %case 'act'
            
        %case 'pass'
            
        %case 'diff'
            
        otherwise
            error(['Unknown mechanism for reaction ',r.id])
    end
    r.exch = logical(exc(i));
    r.sb = logical(sub(i));
    model.kinetic(i) = r;
end

% Construction of kinetic model parametrization structure

dEdk = cell(nr,1);
%dEdkT = dEdk;
A = cell(nr,1);
B = cell(nr,1);
dkdc = cell(nr,1);
vmap = zeros(nr,1);
vdef = cell(nr,1);
vup = ~model.param;
kblocks = zeros(1,nr+1);
exch = [model.kinetic.exch]';
sub = [model.kinetic.sb]';

for i = 1:nr
    dEdk{i} = model.kinetic(i).dEdk;
    le = length(dEdk{i}(:,1));
    ishf = 1:le;
    ishf = ishf(:);
    ishf = reshape(ishf,sqrt(le),sqrt(le));
    ishf = ishf';
    ishf = ishf(:);
    dEdk{i} = dEdk{i}(ishf,:);
    %t = dEdk{i};
    %l1 = length(t(:));
    %le = sqrt(le);
    %l1 = l1/le;
    %dEdkT{i} = reshape(t,le,l1);
    A{i} = zeros(sqrt(length(model.kinetic(i).dEdk(:,1))));
    A{i}(1,:) = 1;
    B{i} = A{i}(:,1);
    dkdc{i} = [full(~any(model.kinetic(i).dkdc,1));model.kinetic(i).dkdc];
    kblocks(i+1) = kblocks(i)+model.kinetic(i).nk;
    vdef{i} = model.kinetic(i).vdef;
end
nc = nm;
%kblocks(1) = 1;
S = model.S;
model.p = struct();
model.p.dEdk = dEdk;
%model.p.dEdkT = dEdkT;
model.p.dkdc = dkdc;
model.p.vmap = vmap;
model.p.vup = vup;
model.p.kblocks = kblocks;
model.p.nc = nc;
model.p.S = S;
model.p.vdef = vdef;
model.p.kid = rs;
model.p.E = A;
model.p.B = B;
model.p.nk = kblocks(end);
model.p.exch = exch;
model.p.sub = sub;
model.p.N = N;
model.p.vfind = vfind;
%model.p = struct('dEdk',dEdk,'dkdc',dkdc,'vmap',vmap,'vup',vup,'kblocks',kblocks,'nc',nc,'S',S);
end



