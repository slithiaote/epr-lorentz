function [spec] = LorentzModel_simulate(B,par)
  % LorentzModel with area C, location Br, full width half max FWHM
  % Finite difference / modulation amplitude MA
  % Additive white Gaussian noise of standard deviation v
  % FWHM, Br and MA are in the same units as B
  FWHM_2 = par.FWHM / 2;
  spec = par.C / FWHM_2 * (LorentzShape((B-par.Br+par.MA/2) ./ FWHM_2) - LorentzShape((B-par.Br-par.MA/2)./FWHM_2))+randn(size(B))*par.v;
end
