package main

import (
	"bytes"
	"encoding/json"
	"fmt"
	"io/ioutil"
	"math"
	"math/cmplx"
	"path"
	"runtime"
	"time"
	"unsafe"

	"github.com/gonum/blas/cblas128"
	"github.com/mjibson/go-dsp/fft"
)

const minusI complex128 = complex(0, -1)

func floatToComplex(s []float64) []complex128 {
	result := make([]complex128, len(s))
	// Don't use range to avoid copying data.
	for i := 0; i < len(s); i++ {
		result[i] = complex(s[i], 0.0)
	}
	return result
}

func LoadWhale() ([]complex128, error) {
	// H/T: andrewbrookins.com/tech/golang-get-directory-of-the-current-file/
	_, filename, _, _ := runtime.Caller(1)
	whale_path := path.Join(path.Dir(filename), "bluewhale.json")

	contents, err := ioutil.ReadFile(whale_path)
	if err != nil {
		return nil, err
	}

	var float64_call []float64
	json.Unmarshal(contents, &float64_call)
	return floatToComplex(float64_call), nil
}

// This is like numpy.ones. Styled after `bytes.Repeat`.
func Repeat(val complex128, count int) []complex128 {
	result := make([]complex128, count)
	result[0] = val
	slice_ptr := 1
	for slice_ptr < len(result) {
		copy(result[slice_ptr:], result[:slice_ptr])
		slice_ptr *= 2
	}
	return result
}

// See src/runtime/slice.go (as of df027ac)
type sliceStruct struct {
	Array unsafe.Pointer
	Len   int
	Cap   int
}

func RepeatAsBytes(val complex128, count int) []complex128 {
	b := *(*[16]byte)(unsafe.Pointer(&val))
	b_repeat := bytes.Repeat(b[:], count)
	b_struct := *(*sliceStruct)(unsafe.Pointer(&b_repeat))

	complex_struct := sliceStruct{
		Array: b_struct.Array,
		Len:   count,
		Cap:   count,
	}
	return *(*[]complex128)(unsafe.Pointer(&complex_struct))
}

func dftKernel(t, s complex128) complex128 {
	return cmplx.Exp(minusI * t * s)
}

func dftKernelVectorized(t complex128, s cblas128.Vector) cblas128.Vector {
	N := len(s.Data)
	result := cblas128.Vector{Inc: 1, Data: make([]complex128, N)}
	cblas128.Axpy(N, t*minusI, s, result)
	// TODO(dhermes): Determine if there is a way to vectorize this
	//                operation.
	// Don't use range to avoid copying data.
	for i := 0; i < N; i++ {
		result.Data[i] = cmplx.Exp(result.Data[i])
	}
	return result
}

func computeFHat(f, t, s []complex128) []complex128 {
	N := len(f)
	f_hat := make([]complex128, N)
	f_vec := cblas128.Vector{Inc: 1, Data: f}
	s_vec := cblas128.Vector{Inc: 1, Data: s}
	for k := 0; k < N; k++ {
		v := dftKernelVectorized(t[k], s_vec)
		f_hat[k] = cblas128.Dotu(N, v, f_vec)
	}
	return f_hat
}

func getDFTData(N int) ([]complex128, []complex128) {
	t := make([]complex128, N)
	s := make([]complex128, N)
	for k := 0; k < N; k++ {
		t[k] = 2.0 * math.Pi * complex(float64(k), 0)
		s[k] = complex(float64(k)/float64(N), 0)
	}
	return t, s
}

func PrintMatrix(m cblas128.General) {
	for row := 0; row < m.Rows; row++ {
		index := row
		fmt.Printf("[")
		for col := 0; col < m.Cols; col++ {
			fmt.Printf("%v", m.Data[index])
			fmt.Printf(" ")
			index += m.Stride
		}
		fmt.Println("]")
	}
}

func CopyColumn(m cblas128.General, col int, v cblas128.Vector) {
	v_tmp := cblas128.Vector{Inc: 1, Data: m.Data[col*m.Stride:]}
	cblas128.Copy(m.Rows, v, v_tmp)
}

func CopyRow(m cblas128.General, row int, v cblas128.Vector) {
	v_tmp := cblas128.Vector{Inc: m.Stride, Data: m.Data[row:]}
	cblas128.Copy(m.Cols, v, v_tmp)
}

func CopyDiagonal(m cblas128.General, diag int, v cblas128.Vector) {
	diag_last_row := m.Cols - diag
	if diag_last_row > m.Rows {
		diag_last_row = m.Rows
	}

	var diag_start, diag_len int
	if diag >= 0 {
		diag_start = diag * m.Stride
		diag_len = diag_last_row
	} else {
		diag_start = -diag
		diag_len = diag_last_row + diag // Really - (-diag)
	}
	v_tmp := cblas128.Vector{Inc: 1 + m.Stride, Data: m.Data[diag_start:]}
	cblas128.Copy(diag_len, v, v_tmp)
}

func compareNaive() {
	blue_whale_call, err := LoadWhale()
	if err != nil {
		fmt.Printf("File error: %v\n", err)
	}

	N := len(blue_whale_call)
	t, s := getDFTData(N)

	start := time.Now()
	dft_whale_call := fft.FFT(blue_whale_call)
	elapsed_fft := time.Since(start)

	start = time.Now()
	computed_f_hat := computeFHat(blue_whale_call, t, s)
	elapsed_naive := time.Since(start)

	fmt.Println("N:", N)
	fmt.Println("FFT took:", elapsed_fft)
	fmt.Println("Naive DFT took:", elapsed_naive)

	err_vals := cblas128.Vector{Inc: 1, Data: dft_whale_call}
	c := cblas128.Vector{Inc: 1, Data: computed_f_hat}
	// err_vals -= computed_f_hat
	cblas128.Axpy(N, -1.0, c, err_vals)

	fmt.Println("||e||_2:", cblas128.Nrm2(N, err_vals))
}

func matrixManipulation() {
	rows, cols := 3, 4
	m_data := make([]complex128, rows*cols)
	m := cblas128.General{Rows: rows, Cols: cols, Stride: rows, Data: m_data}

	fmt.Println("Before:")
	PrintMatrix(m)

	s := []complex128{1.0, 2.0, 3.0, 4.0}
	v := cblas128.Vector{Inc: 1, Data: s}
	CopyRow(m, 2, v)
	fmt.Println("After Row 2:")
	PrintMatrix(m)
}

func main() {
	// See: http://godoc.org/github.com/gonum/blas
	// sudo apt-get install libopenblas-dev
	// CGO_LDFLAGS="-L/usr/lib/libopenblas.so -lopenblas" go install github.com/gonum/blas/cgo
	val := complex(4.2, 13.37)
	fmt.Println(RepeatAsBytes(val, 10))
}
