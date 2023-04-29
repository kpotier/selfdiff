package lammpstrj

import (
	"bufio"
	"bytes"
	"errors"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	"github.com/kpotier/selfdiff/pkg/msd"
)

// Conv converts x y z into xu yu zu from a Lammps Trajectory. See Lammps
// Documentation for the meaning of xu yu and zu.
func Conv(c *msd.Conv) error {
	f, err := os.Open(c.Traj)
	if err != nil {
		return err
	}
	defer f.Close()
	r := bufio.NewReader(f)

	out, err := os.Create(c.Out)
	if err != nil {
		return err
	}
	defer out.Close()

	var (
		colsB   []byte
		cols    [3]int
		colsTot int

		box  [3]float64 // Box of the current configuration
		box2 [3]float64 // Box of the current configuration divided by 2
	)

	corr := make([][3]float64, c.AtTot)    // Correction (incrementation)
	lastXYZ := make([][3]float64, c.AtTot) // Last configuration

	// For the first configuration
	box, box2, err = header(c, r, out)
	if err != nil {
		return fmt.Errorf("header: %w", err)
	}

	colsB, cols, colsTot, lastXYZ, err = firstCfg(c, r, out, box)
	if err != nil {
		return fmt.Errorf("firstCfg: %w", err)
	}

	_, err = r.ReadByte()
	if err != nil {
		if errors.Is(err, io.EOF) {
			return nil
		}
		return err
	}
	r.UnreadByte()

	// For each configuration
	//cfgf := 0
	for {
		//cfgf++
		box, box2, err = header(c, r, out)
		if err != nil {
			return fmt.Errorf("header: %w", err)
		}

		r.ReadSlice('\n')
		out.Write(colsB)

		// For each atom
		for a := 0; a < c.AtTot; a++ {
			l, _ := r.ReadSlice('\n')

			fields := strings.Fields(string(l))
			if len(fields) != colsTot {
				return fmt.Errorf("number of columns don't match")
			}

			//var oldk [3]float64
			for k := 0; k < 3; k++ {
				xyz, _ := strconv.ParseFloat(fields[cols[k]], 64)
				xyz += corr[a][k]
				//oldk[k] = xyz

				dist := lastXYZ[a][k] - xyz
				if dist > box2[k] {
					corr[a][k] += box[k]
					xyz += box[k]
				} else if dist < -box2[k] {
					corr[a][k] -= box[k]
					xyz -= box[k]
				}
				lastXYZ[a][k] = xyz
			}

			//out2, _ := os.OpenFile(fmt.Sprint("./test/", a), os.O_APPEND|os.O_CREATE|os.O_RDWR, 0666)
			//out2.Write([]byte(fmt.Sprintln(cfgf, a, "0", lastXYZ[a][0], lastXYZ[a][1], lastXYZ[a][2], oldk[0], oldk[1], oldk[2])))
			//out2.Close()

			var bytes []byte
			for k, v := range fields {
				switch k {
				case cols[0]:
					bytes = strconv.AppendFloat(bytes, lastXYZ[a][0], 'g', -1, 64)
				case cols[1]:
					bytes = strconv.AppendFloat(bytes, lastXYZ[a][1], 'g', -1, 64)
				case cols[2]:
					bytes = strconv.AppendFloat(bytes, lastXYZ[a][2], 'g', -1, 64)
				default:
					bytes = append(bytes, []byte(v)...)
				}
				bytes = append(bytes, ' ')
			}
			bytes = append(bytes, '\n')
			out.Write(bytes)
		}

		_, err := r.ReadByte()
		if err != nil {
			if errors.Is(err, io.EOF) {
				break
			}
			return err
		}
		r.UnreadByte()
	}

	return nil
}

// readSlice reads until \n and writes it into a file. It also returns the line
// that have been read.
func readSlice(r *bufio.Reader, w io.Writer) []byte {
	b, _ := r.ReadSlice('\n')
	w.Write(b)
	return b
}

// header corresponds to the lines specific to a Lammps trajectory file. It
// contains the size of the box.
func header(c *msd.Conv, r *bufio.Reader, w io.Writer) (box [3]float64, box2 [3]float64, err error) {
	for l := 0; l < 5; l++ {
		readSlice(r, w)
	}

	// Size of the box
	for k := 0; k < 3; k++ {
		b := readSlice(r, w)

		fields := strings.Fields(string(b))
		if len(fields) != 2 {
			err = fmt.Errorf("unable to get the size of the box")
			return
		}

		lmin, _ := strconv.ParseFloat(fields[0], 64)
		lmax, _ := strconv.ParseFloat(fields[1], 64)

		box[k] = lmax - lmin
		box2[k] = box[k] / 2.
	}

	return
}

// firstCfg is the routine that will be performed for the first configuration.
// This configuration will be used as a reference for the other configurations.
func firstCfg(c *msd.Conv, r *bufio.Reader, w io.Writer, box [3]float64) (colsB []byte, cols [3]int, colsTot int, lastXYZ [][3]float64, err error) {
	var (
		found int
		buf   bytes.Buffer
	)

	b, _ := r.ReadSlice('\n')
	fields := strings.Fields(string(b))

	if len(fields) <= 2 {
		err = fmt.Errorf("not enough columns")
		return
	}

	buf.WriteString(fields[0])
	buf.WriteByte(' ')
	buf.WriteString(fields[1])
	buf.WriteByte(' ')

	fields = fields[2:] // Omission of ITEM: ATOMS
	colsTot = len(fields)

	for k, v := range fields {
		switch v {
		case "x":
			cols[0] = k
			buf.WriteString("xu") // unwrapped (see Lammps doc)
		case "y":
			cols[1] = k
			buf.WriteString("yu")
		case "z":
			cols[2] = k
			buf.WriteString("zu")
		default:
			buf.WriteString(v)
			buf.WriteByte(' ')
			continue
		}
		buf.WriteByte(' ')
		found++
	}
	buf.WriteByte('\n')
	w.Write(buf.Bytes())
	colsB = make([]byte, buf.Len())
	copy(colsB, buf.Bytes())

	if found < 3 {
		err = fmt.Errorf("cannot find the columns x, y, and z")
		return
	}

	// Check PBC for each atom in each molecule
	for m := 0; m < c.Mol; m++ {
		var lastXYZMol [3]float64

		for a := 0; a < c.At; a++ {
			b, _ := r.ReadSlice('\n')

			fields := strings.Fields(string(b))
			if len(fields) != colsTot {
				err = fmt.Errorf("number of columns don't match")
				return
			}

			for k := 0; k < 3; k++ {
				if a == 0 {
					lastXYZMol[k], _ = strconv.ParseFloat(fields[cols[k]], 64)
					continue
				}

				xyz, _ := strconv.ParseFloat(fields[cols[k]], 64)
				dist := lastXYZMol[k] - xyz
				if dist > c.Dist[k] {
					xyz += box[k]
				} else if dist < -c.Dist[k] {
					xyz -= box[k]
				}

				lastXYZMol[k] = xyz
			}

			lastXYZ = append(lastXYZ, lastXYZMol)

			var bytes []byte
			for k, v := range fields {
				switch k {
				case cols[0]:
					bytes = strconv.AppendFloat(bytes, lastXYZMol[0], 'g', -1, 64)
				case cols[1]:
					bytes = strconv.AppendFloat(bytes, lastXYZMol[1], 'g', -1, 64)
				case cols[2]:
					bytes = strconv.AppendFloat(bytes, lastXYZMol[2], 'g', -1, 64)
				default:
					bytes = append(bytes, []byte(v)...)
				}
				bytes = append(bytes, ' ')
			}
			bytes = append(bytes, '\n')
			w.Write(bytes)
		}
	}

	return
}
