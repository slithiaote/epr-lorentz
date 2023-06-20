function [parLS] = LorentzModel_ls(XmT,spec,init)
    % Least-squares estimate
    f = @(t)(norm(spec - LorentzModel_simulate(XmT,vect_to_par(t))));
    parLS = vect_to_par(fminsearch(f,par_to_vect(init)));
    parLS.v = std(spec - LorentzModel_simulate(XmT,parLS));
end

% Local functions

function [par] = vect_to_par(t)
% Convert from vector to parameter struct
par.v = 0; par.C = t(2); par.Br = t(3); par.FWHM = t(4); par.MA = t(5);
end

function [t] = par_to_vect(par)
% Convert from parameter struct to vector
t = [par.v; par.C; par.Br; par.FWHM; par.MA];
end
