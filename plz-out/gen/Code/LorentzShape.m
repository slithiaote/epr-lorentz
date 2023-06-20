function [spec] = LorentzShape(B)
  % Normalized Lorentz shape, area=1, width=1
  spec = (B.^2+1).^(-1) ./ pi;
end
