# K-FIT for kinetic model parameter estimation
- You want to run “modcompile.m” first to create the model file form your excel inputs. Next run “kineticestimate.m”. It reports a “res” file containing the lack-of-fit, predictions and values of the optimization variables at the best solution. 
- The folder also contains three examples: a toy model, medium model, and a core model. Before you run each model, please change the variable “cind” in “rescalc.m” to 14 for the toy model, 1 for the medium model, and 69 for the core model. this is the index pf the carbon substrate uptake flux. 