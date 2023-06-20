function [e] = LorentzModelDeri_mleerror(B,par)
a = par.C;    % a: area under ab curve = par.C
s = par.FWHM; % s: sigma or lw   (s = sqrt(3)*p2p) = par.FWHM
m = par.Br;   % m: mu or Br = par.Br
v = par.v;    % v: noise variance = par.v
% x = B;        % x: B magnetic field = B
% y = spec;     % y: spectrum = spec
% v=sqrt(v);
% DerL = - a*s*(x-m)/(pi*((x-m)^2+s^2/4)^2)
% lnL = -(y+a*s*(x-m)/(pi*((x-m)^2+s^2/4)^2))^2/(2*v^2)-n*ln(2*pi*v^2)/2
% hessian(lnL, [a, m, s, v])

xm = B-par.Br;
xm2 = (B-par.Br).^2;
xm2s24 = xm2+s^2/4;

FIM1(1,1) = sum( s^2/pi^2/v^2 * xm2 ./ xm2s24.^4   ); % -e_daa
FIM1(1,2) = sum( s * (4*a*s*xm2-a*s*xm2s24) ./ (pi*xm2s24.^3) .*(B-par.Br)./(pi*v^2*xm2s24.^2)  ); % -e_dami
FIM1(1,3) = sum(s*a/pi* ( xm2 ./ xm2s24.^2 - xm2*s^2 ./ xm2s24.^3 ) ./ (pi*v^2*xm2s24.^2) ); % -e_dasi
FIM1(1,4) = 0;
FIM1(2,1) = FIM1(1,2);
FIM1(2,2) = sum( ((4*a*s*xm2 -a*s*xm2s24) ./(pi*xm2s24.^3)).^2 /v^2 );
FIM1(2,3) = sum( (a*xm./(pi*xm2s24.^2)-a*xm*s^2./(pi*xm2s24.^3)) .*((4*a*xm2*s-a*s*xm2s24)./(pi*xm2s24.^3)) /v^2); % -e_dmsi
FIM1(2,4) = 0;
FIM1(3,1) = FIM1(1,3);
FIM1(3,2) = FIM1(2,3);
FIM1(3,3) = sum( (a*xm./(pi*xm2s24.^2)-a*xm*s^2./(pi*xm2s24.^3)).^2 /v^2 ); % -e_dssi
FIM1(3,4) = 0;
FIM1(4,1) = 0;
FIM1(4,2) = 0;
FIM1(4,3) = 0;
FIM1(4,4) = 2*length(B)/v^2;

if rcond(FIM1) > 1e-20
    cov_mle = inv(FIM1);
    e = par;
    e.eBr = sqrt(cov_mle(2,2));
    e.eFWHM = sqrt(cov_mle(3,3));
    e.eC = sqrt(cov_mle(1,1));
    e.ev = sqrt(cov_mle(4,4));
else
    e = par; e.eBr = -1; e.eFWHM = -1; e.eC = -1; e.ev = -1;
end

