function [parLM] = LorentzModel_initial(B,spec)
    [y_max, index_max]=max(spec);[y_min, index_min]=min(spec);
    parLM.Br = (B(index_max)+B(index_min))/2;
    parLM.FWHM = (B(index_min)-B(index_max))*sqrt(3);
    parLM.C = 1; parLM.v = 0;
    parLM.MA = 0.5; % arbitrary initial guess
    spec_temp = LorentzModel_simulate(B, parLM);
    parLM.C = spec_temp(:) \ spec(:); % Linear regression
    parLM.v = std(spec - spec_temp * parLM.C);
end
