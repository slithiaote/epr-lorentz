function [parLMD] = LorentzModelDeri(XmT,spec,parLMD_initial)
    if ~exist('parLMD_initial','var')
        parLMD_initial = LorentzModelDeri_initial  (XmT,spec);
    end
    parLMD_mle = LorentzModelDeri_mle      (XmT,spec,parLMD_initial);
    parLMD = LorentzModelDeri_mleerror (XmT,parLMD_mle);
