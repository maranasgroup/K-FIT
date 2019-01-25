function [ model ] = loaddata( model,data_filename )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[~,~,rawx] = xlsread(data_filename);
rawx = rawx(2:end,:);
type = rawx(:,5);
raw = rawx(ismember(type,'flx'),:);
muts = raw(:,2);
%d.err = [];
conds = unique(muts);
ncond = length(conds);
rxns = model.rid;
d.id = cell(ncond,1);
d.flx = cell(ncond,1);
d.err = cell(ncond,1);
d.gko = cell(ncond,1);
d.rmap = cell(ncond,1);
d.vpert = zeros(length(rxns),ncond);
measmaps = unique(raw(:,end));
rmps = zeros(length(measmaps),length(rxns));
for i = 1:length(measmaps)
    rz = regexp(measmaps{i},';','Split');
    rmps(i,:) = double(ismember(rxns,rz)');
end



% Wild-type fluxes
% WT index is always 1

mutmask = ismember(muts,{'WT'});
rs = raw(mutmask,1);
mmp = raw(mutmask,6);
vko = ones(length(rxns),1);
d.id{1} = rs;
d.gko{1} = vko;
d.vpert(:,1) = vko;
d.flx{1} = cell2mat(raw(mutmask,3));
d.err{1} = cell2mat(raw(mutmask,4));
d.rmap{1} = zeros(length(d.flx{1}),length(rxns));
for i = 1:length(rs)
    k = ismember(measmaps,mmp{i});
    d.rmap{1}(i,:) = rmps(k,:);
end
%d.id = rs;


%csum = 0;
%[~,flxidx] = ismember(rs,rxns);
%d.idx = flxidx;

nr = length(rxns);

raw(mutmask,:) = [];
%csum = csum+nr;
muts = raw(:,2);
mutlist = unique(muts);
nmut = length(mutlist);
%d.vpert = zeros(nr,nmut+1);
%d.vpert(:,1) = vko;
%d.rmap = cell(

for i = 1:nmut
    mutmask = ismember(muts,mutlist(i));
    rs = raw(mutmask,1);
    mmp = raw(mutmask,6);
    d.id{i+1} = rs;
    %v = vko;
    %v(ismember(rxns,mutlist(i))) = 0;
    mz = regexp(mutlist{i},';','Split');
    gko = full(~ismember(rxns,mz));
    d.gko{i+1} = gko;
    d.vpert(:,i+1) = gko;
    
    %[~,flxidx] = ismember(rs,rxns);
    %flxidx = flxidx + csum;
    %d.idx = [d.idx;flxidx];
    d.flx{i+1} = cell2mat(raw(mutmask,3));
    d.err{i+1} = cell2mat(raw(mutmask,4));
    d.rmap{i+1} = zeros(length(rs),nr);
    for j = 1:length(rs)
        k = ismember(measmaps,mmp{j});
        d.rmap{i+1}(j,:) = rmps(k,:);
    end
    
end

%reversibilities

raw = rawx(ismember(type,'rev'),:);
inds = 1:length(model.rid);
inds = inds';
nrevs = length(raw(:,1));
d.revs = cell2mat(raw(:,3));
d.rerr = cell2mat(raw(:,4));
d.ridx = zeros(nrevs,1);
for i = 1:nrevs
    in = ismember(model.rid,raw(i,1));
    d.ridx(i) = double(in')*inds;
end






model.d = d;

end

