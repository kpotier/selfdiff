package cfg

import (
	"bufio"
	"fmt"
	"os"
	"path/filepath"
	"strings"

	"github.com/kpotier/selfdiff/pkg/msd"
	lammpstrjMSD "github.com/kpotier/selfdiff/pkg/msd/lammpstrj"
	"github.com/kpotier/selfdiff/pkg/vac"
	lammpstrjVAC "github.com/kpotier/selfdiff/pkg/vac/lammpstrj"

	"gopkg.in/yaml.v3"
)

// Method is the method of calculation
type Method string

// Here are the following accepted methods. MSD means mean squared
// displacement. VAC means velocity auto correlation.
var (
	MMSD Method = "msd"
	MVAC Method = "vac"
)

// Type is the type of the trajectory
type Type string

// Here are the accepted types. Lammpstrj is a Lammps Trajectory file.
var (
	TLammpstrj Type = "lammpstrj"
)

// Cfg is a structure containing the parameters specified in the configuration
// file. It can be instanced through the NewCfg method or by "hand". If it is
// instanced by hand, please use the Check method to check if the Cfg meets the
// requirements.
type Cfg struct {
	// Traj is the file containing the configurations
	Traj string `yaml:"traj"`

	// Type is the type of trajectory (e.g: lammpstrj)
	Type Type `yaml:"type"`

	// Method is the method of calculation
	Method Method `yaml:"method"`

	// PBC specifies if the periodic boundary conditions are used in the above file
	PBC bool `yaml:"pbc"`

	// Start is the first configuration that will be read. It must start be
	// greater or equal to 0
	Start int `yaml:"start"`

	// End is the last configuration that will be read. It means that if
	// End = 1000, the 1000th configuration will be read
	End int `yaml:"end"`

	// Mem is the number of configurations that will be put in memory. If it
	// is set to 3, the last 3 configurations will be put in memory (the most used)
	Mem int `yaml:"mem"`

	// Mol is the number of molecules in one configuration
	Mol int `yaml:"mol"`

	// At is the number of atoms in one molecule
	At int `yaml:"at"`

	// Dist is the largest distance between two atoms in one molecule
	Dist [3]float64 `yaml:"msdDist"`

	// Masses are the masses of each atoms in one molecule
	Masses []float64 `yaml:"masses"`

	// Dt is the timestep in whatever unit you want
	Dt float64 `yaml:"dt"`
}

// New opens and decodes the specified configuration file. The file must be
// a YAML file. This method automatically calls the Check method to check the
// integrity of Cfg.
func New(path string) (*Cfg, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	var c Cfg
	r := bufio.NewReader(f)
	dec := yaml.NewDecoder(r)
	err = dec.Decode(&c)
	if err != nil {
		return nil, err
	}

	err = c.Check()
	if err != nil {
		return nil, fmt.Errorf("Check: %w", err)
	}

	return &c, nil
}

// Check checks if Cfg is correct. It returns an error if a field doesn't meet
// the requirements.
func (c *Cfg) Check() error {
	if c.Start < 0 {
		return fmt.Errorf("Start must be greater or equal to 0")
	}

	if c.End <= c.Start {
		return fmt.Errorf("End cannot be lower or equal to Start")
	}

	if (c.End - c.Start) == 1 {
		return fmt.Errorf("End-Start must not be equal to 1")
	}

	if c.Mem > (c.End-c.Start) || c.Mem < 0 {
		return fmt.Errorf("Mem cannot be lower than 0 or greater than End-Start")
	}

	if c.Mol <= 0 || c.At <= 0 {
		return fmt.Errorf("Mol or Att cannot be lower or equal to 0")
	}

	if len(c.Masses) != c.At {
		return fmt.Errorf("the length of the masses slice is not equal to At")
	}

	if c.Dt == 0 {
		return fmt.Errorf("Dt cannot be lower or equal to 0")
	}

	return nil
}

// Conv returns the conversion method. It is usefull to get the non PBC
// trajectory from a PBC trajectory.
func (c *Cfg) Conv() error {
	if !c.PBC {
		return fmt.Errorf("pbc set to false")
	}

	if c.Method != MMSD {
		return fmt.Errorf("msd method is required")
	}

	ext := filepath.Ext(c.Traj)
	filename := strings.TrimSuffix(c.Traj, ext)
	newTraj := fmt.Sprint(filename, "_nopbc", ext)

	conv := &msd.Conv{Traj: c.Traj, Out: newTraj, At: c.At, Mol: c.Mol, AtTot: (c.At * c.Mol), Dist: c.Dist}

	var err error
	switch c.Type {
	case TLammpstrj:
		err = lammpstrjMSD.Conv(conv)
	default:
		err = fmt.Errorf("unsupported type")
	}

	if err != nil {
		return err
	}

	c.PBC = false
	c.Traj = newTraj

	return nil
}

// MSD calculates the mean squared displacement.
func (c *Cfg) MSD() (err error) {
	if c.PBC {
		return fmt.Errorf("pbc set to true")
	}

	if c.Method != MMSD {
		return fmt.Errorf("msd method is required")
	}

	out := fmt.Sprint(c.Traj, "_msd.out")
	msd := &msd.MSD{Traj: c.Traj, Out: out, Start: c.Start, End: c.End, Mem: c.Mem, At: c.At, Mol: c.Mol, Masses: c.Masses, Dt: c.Dt}

	switch c.Type {
	case TLammpstrj:
		msd.Method = lammpstrjMSD.New(msd)
	default:
		err = fmt.Errorf("unsupported type")
		return
	}

	err = msd.Perform()
	if err != nil {
		return
	}

	err = msd.Write()
	return
}

// VAC calculates the velocity autocorrelation function.
func (c *Cfg) VAC() (err error) {
	if c.Method != MVAC {
		return fmt.Errorf("vac method is required")
	}

	out := fmt.Sprint(c.Traj, "_vac.out")
	vac := &vac.VAC{Traj: c.Traj, Out: out, Start: c.Start, End: c.End, Mem: c.Mem, At: c.At, Mol: c.Mol, Masses: c.Masses, Dt: c.Dt}

	switch c.Type {
	case TLammpstrj:
		vac.Method = lammpstrjVAC.New(vac)
	default:
		err = fmt.Errorf("unsupported type")
		return
	}

	err = vac.Perform()
	if err != nil {
		return
	}

	err = vac.Write()
	return
}
