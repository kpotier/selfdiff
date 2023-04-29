package vac

import (
	"fmt"
	"os"
)

// Method is an interface that will be used by the modules.
type Method interface {
	Read() error
	GetCfg(int) ([][3]float64, error)
	End() error
}

// VAC structure is a structure containing information that will be used by the
// modules. It contains the position of the first configuration, the position of
// the last configuration, etc.
type VAC struct {
	Method Method

	Traj string
	Out  string

	Start int
	End   int
	Mem   int

	Tot    int
	MemPos int // Position of the configurations that are in the memory
	AtTot  int

	At     int
	Mol    int
	Masses []float64
	Dt     float64

	Res    []float64
	ResDiv float64
	Int    float64
}

// Perform performs the mean squared displacement.
func (m *VAC) Perform() (err error) {
	m.Tot = m.End - m.Start
	m.Res = make([]float64, m.Tot-1)
	m.AtTot = m.At * m.Mol
	m.MemPos = m.Tot - m.Mem

	err = m.Method.Read()
	if err != nil {
		return
	}

	for i := 0; i < m.Tot-1; i++ {
		fmt.Print("\r> Step ", i+1, "/", m.Tot-1)

		var icfg [][3]float64
		icfg, err = m.Method.GetCfg(i)
		if err != nil {
			return
		}

		for mol := 0; mol < m.Mol; mol++ {
			for k := 0; k < 3; k++ {
				m.ResDiv += icfg[mol][k] * icfg[mol][k]
			}
		}

		for j := i + 1; j < m.Tot; j++ {
			var tcfg [][3]float64

			tcfg, err = m.Method.GetCfg(j)
			if err != nil {
				return
			}

			for mol := 0; mol < m.Mol; mol++ {
				for k := 0; k < 3; k++ {
					m.Res[j-i-1] += icfg[mol][k] * tcfg[mol][k]
				}
			}
		}
	}

	fmt.Print("\033[2K\033[1G")

	m.ResDiv /= float64((m.Tot-1)*m.Mol*3) / 2.
	for i := 0; i < m.Tot-1; i++ {
		m.Res[i] /= float64((m.Tot-1-i)*m.Mol*3) / 2.
		m.Int += m.Res[i]
	}

	return
}

// Write writes the results into Out.
func (m *VAC) Write() error {
	f, err := os.Create(m.Out)
	if err != nil {
		return err
	}

	fmt.Fprintln(f, "Integral", m.Int)
	for i := 0; i < m.Tot-1; i++ {
		fmt.Fprintln(f, float64(i+1)*m.Dt, m.Res[i], m.ResDiv)
	}

	return nil
}
