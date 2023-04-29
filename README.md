# selfdiff [![go.dev reference](https://img.shields.io/badge/go.dev-reference-007d9c?logo=go&logoColor=white&style=flat-square)](https://pkg.go.dev/github.com/kpotier/selfdiff)
Tools to calculate the self diffusion coefficient. It is able to calculate the mean squared displacement and the velocity autocorrelation function.

### Supported formats

1. Lammps Trajectory (.lammpstrj)

### Usage

1. Install ```Go 1.13```.

2. Go to the ```cmd/msd``` or ```cmd/vac``` directory.

3. Execute ```go build``` or ```go install```.

### Examples
1. ```msd config.yaml```
   This command will calculate the mean squared displacement. Examples of the config.yaml file can be found in the ```test``` directory.

2. ```vac config.yaml```
   This command will calculate the velocity autocorrelation function. Examples of the config.yaml file can be found in the ```test``` directory.
