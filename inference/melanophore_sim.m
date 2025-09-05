function melanophore_sim(alpha,beta,mutant,randomseed,savefolder)

sim_params = zebrafish_default_params();
sim_params.MbornAddPar = alpha;
sim_params.MbornPropXtoMfarlb = beta;

switch mutant
    case "wt"
        % do nothing
    case "pfeffer"
        sim_params.cellsBornRandOrig = 0;
    case "nacre"
        sim_params.cellsBornM = 0;
        sim_params.cellsBornRandOrigM = 0;
    case "shady"
        sim_params.birthOption = 3;
    otherwise
        error("No such mutant type");
end

T=45;
tda_params=tda_default_params();
ts=1:45;
simname=sprintf('zebrafish_melanophores_%s_alpha=%.2f,beta=%.2f,rng=%d',mutant,alpha,beta,randomseed);
[~,~,~] = params_to_tda(sim_params,randomseed,T,tda_params,ts,savefolder,simname,0,1);
end
