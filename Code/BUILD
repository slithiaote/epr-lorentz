lepton("LorentzModelDeri", # Lorentz 1D derivative-limit model, with MLE estimate
      ["//Assets:Test-data"], # required for testing on real data
      ["LorentzModelDeri_simulate.m","LMD_simulate.jpg", # Simulation+test
       "LorentzModelDeri_initial.m",    # Initial estimates
       "LorentzModelDeri_mle.m",        # MLE estimates
       "LorentzModelDeri_mleerror.m",   # MLE error estimates
       "LorentzModelDeri.m", # Simplified function
       "LMD_TAM21.jpg"])

lepton("LorentzModel_hessian",[],[]) # Compute Hessian formulae with sympy
lepton("LorentzModel", # Lorentz 1D full model, with MLE estimate
      ["//Assets:Test-data",":LorentzModel_hessian"],
      ["LorentzModel_simulate.m","LorentzShape.m","LM_simulate.jpg","LM_simulate2.jpg",
       "LorentzModel_initial.m", "LorentzModel_ls.m","LorentzModel_mleerror.m",
       "LorentzModel.m"]) # Simulation+test
