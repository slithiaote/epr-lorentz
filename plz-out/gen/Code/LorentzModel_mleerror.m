function prec = LorentzModel_mleerror(B,par)
  % return the precision estimates as a struct
  n=length(B); gamma = par.FWHM/2;
  Abs   = @(B) LorentzShape(B);
  AbsB  = @(B) -2/pi * B ./ (B.^2 + 1).^2;
  AbsBB = @(B) ( 8/pi * B.^2 - 2/pi * (B.^2 + 1) ) ./ (B.^2 + 1).^3;
  BAbsB = @(B) B .* AbsB(B);
  BAbsBB = @(B) B .* AbsBB(B);
  BBAbsBB = @(B) B.^2 .* AbsB(B);
  ds = @(f) f( (B-par.Br+par.MA/2)/gamma)/gamma - f( (B-par.Br-par.MA/2)/gamma)/gamma;
  dms = @(f) f( (B-par.Br+par.MA/2)/gamma)/gamma/2 + f( (B-par.Br-par.MA/2)/gamma)/gamma/2;
  par0 = par; par0.v = 0;
  EPR=LorentzModel_simulate(B,par0);
dvv = -2*n/par.v^2;
dvC = 0; dvB = 0; dvM = 0; dvg = 0; 
dCv = dvC;
dCC = -1 / par.v^2 * sum(ds(Abs).^2);
dCB = -1 / par.v^2 / gamma * sum( -EPR .* ds(AbsB));
dCM =  1 / par.v^2 / gamma * sum( -EPR .* dms(AbsB));
dCg = -1 / par.v^2 / gamma * sum( -EPR .* (ds(Abs)+ds(BAbsB)));
dBv = dvB;
dBC = dCB;
dBB = -par.C / par.v^2 / gamma^2 * sum( par.C * ds(AbsB).^2 );
dBM = -par.C / par.v^2 / gamma^2 * sum( -par.C * dms(AbsB) .* ds(AbsB) );
dBg = -par.C / par.v^2 / gamma^2 * sum( (EPR + par.C*ds(BAbsB)) .* ds(AbsB)  );
dMv = dvM;
dMC = dCM;
dMB = dBM;
dMM = par.C / par.v^2 / gamma^2 * sum(-par.C * dms(AbsB).^2 );
dMg = par.C / par.v^2 / gamma^2 * sum( par.C * (ds(Abs)+ds(BAbsB)) .* dms(AbsB) );
dgv = dvg;
dgC = dCg;
dgB = dBg;
dgM = dMg;
dgg = par.C / par.v^2 / gamma^2 * sum( -par.C*(ds(Abs)+ds(BAbsB)).^2 );
hessian = [dvv, dvC, dvB, dvM, dvg; dCv, dCC, dCB, dCM, dCg; dBv, dBC, dBB, dBM, dBg; dMv, dMC, dMB, dMM, dMg; dgv, dgC, dgB, dgM, dgg];
prec = par;
if rcond(hessian) > 1e-20
    cov = inv(-hessian);
    prec.ev = sqrt(cov(1,1)); prec.eC = sqrt(cov(2,2)); prec.eBr = sqrt(cov(3,3)); prec.eMA = sqrt(cov(4,4)); prec.eFWHM = 2*sqrt(cov(5,5)); 
else
    prec.ev = -1; prec.eC = -1; prec.eBr = -1; prec.eMA = -1; prec.eFWHM = -1;
end
