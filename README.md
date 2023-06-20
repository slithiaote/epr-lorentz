# Signal processing and data analysis for EPR experiments

This project analyses experimental results obtained during Electron Paramagnetic Resonance experiments.
We provide matlab source code to perform parametric estimation for Lorentz models of EPR spectra.


# Usage


## Browsing the results on Github

Look inside the [`plz-out/gen/`](plz-out/gen) directory . In particular:
  - Parameter estimation for the Lorentz model in derivative limit [here](plz-out/gen/Code/LorentzModelDeri.md).
  - Parameter estimation for the general Lorentz model [here](plz-out/gen/Code/LorentzModel.md).


## Using the functions in Matlab

- Detailed tutorial [here](plz-out/gen/tutorial.md), including application to real data.

or 

- Instructions for advanced users:
  - You should have a working installation of Matlab. The [Easyspin](https://www.easyspin.org/) toolbox is used to read data from Bruker instruments. No additional toolbox is necessary.
  - Download the code [here](archive/refs/heads/main.zip), and unzip.
  - Add the `plz-out/gen/Code/` directory to your matlab path.
  - The `LorentzModelDeri(B,spec,par_init)` function estimates the parameters of the Lorentz model in derivative limit (see [doc](plz-out/gen/Code/LorentzModelDeri.md) for details).
  - The `LorentzModel(B,spec,par_init)` function estimates the parameters of the general Lorentz model (see [doc](plz-out/gen/Code/LorentzModel.md) for details).

Notes : 
- The Lorentz models are intended for single line, baseline corrected, 1D EPR spectra.
- The initial parameters `par_init` are optional, and coarse estimates are used when they are omitted. 
- Precision estimates can be unreliable due to numerical computation errors. In that case, the functions will set precision to `-1`.


# For developpers


## Reproducing the results with Docker

Docker is a framework for deploying containers, which are similar to lightweight virtual machines.

The file [`matlab-lepton-latex`](matlab-lepton-latex) contains the commands needed to build such a container, and thus completely
describes the set of software used in this project (Matlab, Please, Lepton, Python, etc.). 
To use it, you need to install Docker, and then build 
the a Docker image (similar to a VM template) with the following command.

`docker image build -f matlab-lepton-latex -t matlab-lepton .`

To actually run software, you tell Docker to start a container (the virtual machine) from a given image (the template)
and run a command inside this container. For instance, you can run the `echo` command to display some text, or the `bash` command 
to obtain a prompt running inside of the container.

`docker run matlab-lepton echo "Hello world."`

`docker run -it matlab-lepton bash`

In this project, we use the Please build system to manage the different scripts. Running `pleasew` inside the Docker
container will look for `BUILD` files throughout the repository, and execute the commands provided therein to 
generate the contents of the `plz-out/gen/` directory.

`docker run -it --rm--shm-size=512M -v $(pwd):/home/matlab/project matlab-lepton ./pleasew build --show_all_output`

Note that you need to setup your credentials to run Matlab inside the Docker container. One possibility is to set 
the `mlm-license-file` environment variable in Please by editing the `.plz-config.local` file at repository root:
````
[BuildEnv]
mlm-license-file=truc
````

## Notes on the Please build system

This project uses the [Please build system](https://please.build/) to manage the execution of scripts, store intermediate computations in a cache 
and build all the results. The Please system infers what to do from files with the name `BUILD` located in each directory
of the project. Consequently, these `BUILD` files provide an exhaustive list of all the steps 
(called targets, with naming scheme `//package/sub-package:target-name`) 
involved in this project. We refer to this list in the documentation below.

Related files in this project are located in the same directory, and the directory name is used as the `package` name in the
Please system. We use the following packages : 
- `Assets` : read-only data and metadata, typically the raw files obtained from the EPR instrument,
- `Code` : matlab source files necessary to perform signal processing and model fitting,


## Assets = read-only data and metadata
- `//Assets/Test-data` : sample data for testing methods and algorithms <br>
  `TAM21`: TAM sample in atmospheric conditions, single-line low-noise spectrum


## Code = functions that do not produce output files
- //Code:LorentzModelDeri: Lorentz Model, derivative limit case (see [doc](plz-out/gen/Code/LorentzModelDeri.md) for details).
  - `LorentzModelDeri_simulate(B,par)` : simulate EPR spectrum from given parameters
  - `LorentzModelDeri_initial(B,spec)` : compute initial parameter estimates
  - `LorentzModelDeri_mle(B,spec,par_init)` : estimate parameters starting from initial values
  - `LorentzModelDeri_mleerror(B,par)` : compute estimates of the parameter precision
  - `LorentzModelDeri(B,spec,par_init)` : convenience function to estimate all parameters, initial parameters are optional.


- //Code:LorentzModel: Lorentz Model with overmodulation for 1D EPR spectra (see [doc](plz-out/gen/Code/LorentzModel.md) for details).
  - `LorentzModel_simulate(B,par)` : simulate EPR spectrum from given parameters
  - `LorentzModel_initial(B,spec)` : compute initial parameter estimates
  - `LorentzModel_ls(B,spec,par_init)` : estimate parameters starting from initial values
  - `LorentzModel_mleerror(B,par)` : compute estimates of the parameter precision
  - `LorentzModel(B,spec,par_init)` : convenience function to estimate all parameters, initial parameters are optional.



