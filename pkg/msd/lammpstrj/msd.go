package lammpstrj

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
	"strings"

	"github.com/kpotier/selfdiff/pkg/msd"
)

// MSD is a structure specific to a Lammps Trajectory file. It contains the xyz
// coordinates for the configurations that are in the memory and the position of
// the xyz coordinates for the configurations that are not in the memory. It
// also contains the number of columns and the position of the fields xu yu and zu.
type MSD struct {
	*msd.MSD

	f *os.File

	cols    [3]int
	colsTot int

	xyz  [][][3]float64
	xyzB []uint32
}

// New returns an instance of the MSD structure for a Lammps Trajectory file.
func New(c *msd.MSD) *MSD {
	return &MSD{c, nil, [3]int{0, 0, 0}, 0, nil, nil}
}

// Read is part of the MSD interface in the msd package. It reads the
// configurations and put the last ones into memory. It also reads the position
// of the lines of the configurations that won't be put in memory.
func (m *MSD) Read() error {
	var (
		err      error
		bTot     int
		notFirst bool
	)

	m.f, err = os.Open(m.Traj)
	if err != nil {
		return err
	}

	r := bufio.NewReader(m.f)

	// Read and discard the useless configurations
	lines := m.Start * (m.AtTot + 9)
	for l := 0; l < lines; l++ {
		_, bTot = readLine(r, bTot)
	}

	// For the configurations that won't be put in memory
	for c := 0; c < (m.Tot - m.Mem); c++ {
		for l := 0; l < 9; l++ {
			_, bTot = readLine(r, bTot)
		}

		m.xyzB = append(m.xyzB, uint32(bTot))

		for l := 0; l < m.AtTot; l++ {
			_, bTot = readLine(r, bTot)
		}
	}

	// For the configurations that will be put into memory
	for c := 0; c < m.Mem; c++ {
		var xyz [][3]float64

		for l := 0; l < 8; l++ {
			_, bTot = readLine(r, bTot)
		}

		if !notFirst {
			var (
				s     string
				found int
			)

			s, bTot = readLine(r, bTot)
			fields := strings.Fields(s)

			if len(fields) <= 2 {
				return fmt.Errorf("not enough columns")
			}

			fields = fields[2:] // Omission of ITEM: ATOMS
			m.colsTot = len(fields)

			for k, v := range fields {
				switch v {
				case "xu":
					m.cols[0] = k
				case "yu":
					m.cols[1] = k
				case "zu":
					m.cols[2] = k
				default:
					continue
				}
				found++
			}

			if found < 3 {
				return fmt.Errorf("cannot find the columns xu yu, and zu")
			}

			notFirst = true
		} else {
			_, bTot = readLine(r, bTot)
		}

		// We read the position of each atom for each molecule and we determine
		// its center of mass
		for mol := 0; mol < m.Mol; mol++ {
			var tmpXYZ [3]float64
			var mTot float64

			for a := 0; a < m.At; a++ { // Each atom of the molecule m
				var l string
				l, bTot = readLine(r, bTot)

				fields := strings.Fields(l)
				if len(fields) != m.colsTot {
					return fmt.Errorf("number of columns don't match")
				}

				for k := 0; k < 3; k++ {
					pos, _ := strconv.ParseFloat(fields[m.cols[k]], 64)
					tmpXYZ[k] += pos * m.Masses[a]
				}

				mTot += m.Masses[a]
			}

			// Center of mass. Then we add the coordinates into a slice
			for k := 0; k < 3; k++ {
				tmpXYZ[k] /= mTot
			}
			xyz = append(xyz, tmpXYZ)
		}
		m.xyz = append(m.xyz, xyz)
	}

	return nil
}

// GetCfg returns the positions xu yu and zu for a specified configuration.
func (m *MSD) GetCfg(c int) ([][3]float64, error) {
	if c >= m.MemPos {
		return m.xyz[c-m.MemPos], nil
	}

	m.f.Seek(int64(m.xyzB[c]), 0)
	r := bufio.NewReader(m.f)

	// We read the position of each atom for each molecule and we determine
	// its center of mass
	var xyz [][3]float64
	for mol := 0; mol < m.Mol; mol++ {
		var tmpXYZ [3]float64
		var mTot float64

		for a := 0; a < m.At; a++ { // Each atom of the molecule m
			b, _ := r.ReadSlice('\n')

			fields := strings.Fields(string(b)) // Omission of atom type
			if len(fields) != m.colsTot {
				return nil, fmt.Errorf("number of columns don't match")
			}

			for k := 0; k < 3; k++ {
				pos, _ := strconv.ParseFloat(fields[m.cols[k]], 64)
				tmpXYZ[k] += pos * m.Masses[a]
			}

			mTot += m.Masses[a]
		}

		// Center of mass. Then we add the coordinates into a slice
		for k := 0; k < 3; k++ {
			tmpXYZ[k] /= mTot
		}
		xyz = append(xyz, tmpXYZ)
	}

	return xyz, nil
}

// End closes the file opened.
func (m *MSD) End() error {
	return m.f.Close()
}

// readLine reads ONE ligne and returns it with b+(number of bytes in this line).
func readLine(r *bufio.Reader, b int) (string, int) {
	l, _ := r.ReadSlice('\n') // WARNING: ReadSlice doesn't copy l. l will be replaced if another call of ReadSlice is performed
	b += len(l)
	return string(l), b
}
