function [parLM] = LorentzModel(XmT, spec, parLM_initial)
    if ~exist('parLM_initial','var')
        parLM_initial = LorentzModel_initial  (XmT,spec);
    end
    parLM_mle = LorentzModel_ls      (XmT,spec,parLM_initial);
    parLM = LorentzModel_mleerror (XmT,parLM_mle);
