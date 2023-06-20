function [parLMD] = LorentzModelDeri_initial(B,spec)
[y_max, index_max]=max(spec);[y_min, index_min]=min(spec);
parLMD.Br = (B(index_max)+B(index_min))/2;
parLMD.FWHM = (B(index_min)-B(index_max))*sqrt(3);
parLMD.C = 1; parLMD.v = 0;
spec_temp = LorentzModelDeri_simulate(B, parLMD);
parLMD.C = spec_temp(:) \ spec(:); % Linear regression, requires col vectors
parLMD.v = std(spec - spec_temp * parLMD.C);
end
