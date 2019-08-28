INPUTS = [];

INPUTS.da_tot = da_tot;
INPUTS.n_tips = n_tips;
INPUTS.n_mesh = n_mesh;

INPUTS.type_Gs = type_Gs;
INPUTS.type_Hs = type_Hs;

INPUTS.Gs = Gs;
INPUTS.Hs = Hs;

save([NAME,'_inputs'],'NAME','INPUTS')