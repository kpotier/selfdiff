package msd

// Conv is a structure that will be used by the modules. It contains information
// like the number of atoms, the number of molecules.
type Conv struct {
	Traj string
	Out  string

	At    int
	Mol   int
	AtTot int
	Dist  [3]float64
}
