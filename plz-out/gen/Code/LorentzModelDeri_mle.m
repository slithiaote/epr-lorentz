function [parLMD] = LorentzModelDeri_mle(B,spec,initial)
lnL = @(t)(length(spec)/2)*log(2*pi*t(1)^2) + norm(spec - (t(2)*t(3)*(t(4)-B))./(pi*((B-t(4)).^2+(t(3)^2)/4).^2) )^2 / (2*t(1)^2);   
[t]=fminsearch(lnL,[initial.v,initial.C,initial.FWHM,initial.Br]);
parLMD.v=t(1); parLMD.C=t(2);parLMD.FWHM=t(3);parLMD.Br=t(4);
