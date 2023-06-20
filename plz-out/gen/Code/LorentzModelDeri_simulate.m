function [spec] = LorentzModelDeri_simulate(B,parLMD)
BmBr = B - parLMD.Br; % used only to simplify expression
spec = - parLMD.C * parLMD.FWHM / pi * BmBr ./ (1/4*parLMD.FWHM^2 + BmBr .^2 ).^2 + randn(size(B))*parLMD.v;
end
