function [model] = modcompile(model_file,mech_file,data_file)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

model = xls2MFAmodel(model_file);
model = kineticdecomp(model,mech_file);
model = loaddata(model,data_file);

%vref = model.d.flx(1:length(model.rid));
%model = ensemble(model,vref);
model = ensemble(model);
model = ccalcsetup(model);
model = gradcalc(model);
model.options.reinit = true;
model.options.multistarts = 1;
end

