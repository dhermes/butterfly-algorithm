package main

import (
	"bytes"
	"encoding/json"
	"fmt"
	"io/ioutil"
	"math"
	"math/cmplx"
	"path"
	"reflect"
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

func RepeatAsBytes(val complex128, count int) []complex128 {
	b := *(*[16]byte)(unsafe.Pointer(&val))
	b_repeat := bytes.Repeat(b[:], count)
	b_struct := *(*reflect.SliceHeader)(unsafe.Pointer(&b_repeat))

	complex_struct := reflect.SliceHeader{
		Data: b_struct.Data,
		Len:  count,
		Cap:  count,
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

func ComponentWiseMultiply(N int, v cblas128.Vector, w cblas128.Vector) {
	v_index, w_index := 0, 0
	for i := 0; i < N; i++ {
		v.Data[i] *= w.Data[w_index]
		v_index += v.Inc
		w_index += w.Inc
	}
}

func ComponentWiseProduct(N int, v cblas128.Vector, w cblas128.Vector) cblas128.Vector {
	result := cblas128.Vector{Inc: 1, Data: make([]complex128, N)}
	cblas128.Copy(N, v, result)
	ComponentWiseMultiply(N, result, w)
	return result
}

func getBinsAndDeltas(vals []float64, min_val, bin_width float64, num_bins int) ([]int, cblas128.Vector) {
	N := len(vals)

	bin_indices := make([]int, N)
	deltas := cblas128.Vector{Inc: 1, Data: make([]complex128, N)}
	center0 := min_val + 0.5*bin_width
	for i := 0; i < N; i++ {
		bin_index := int(math.Floor((vals[i] - min_val) / bin_width))
		if bin_index == num_bins {
			bin_index = num_bins - 1
		}
		bin_indices[i] = bin_index
		deltas.Data[i] = complex(vals[i]-(center0+float64(bin_index)*bin_width), 0.0)
	}
	return bin_indices, deltas
}

type BinData struct {
	Tau    float64
	Sigma  float64
	Coeffs cblas128.Vector
}

func createInitialData(s []float64, min_s, max_s, tau float64,
	actual_data []complex128, num_bins, M int) []BinData {
	bin_width := (max_s - min_s) / float64(num_bins)

	bin_indices, s_deltas := getBinsAndDeltas(s, min_s, bin_width, num_bins)
	N := len(s)

	sum_parts := cblas128.General{
		Rows:   N,
		Cols:   M,
		Stride: N,
		Data:   make([]complex128, N*M),
	}
	s_vec := cblas128.Vector{Inc: 1, Data: floatToComplex(s)}

	// When alpha = 0 (column 0), we have dft(tau, s) * D(s).
	col_vals := dftKernelVectorized(complex(tau, 0.0), s_vec)
	ComponentWiseMultiply(N, col_vals, cblas128.Vector{Inc: 1, Data: actual_data})
	CopyColumn(sum_parts, 0, col_vals)

	for alpha := 1; alpha < M; alpha++ {
		// Update one factor of (s - sigma) in (s - sigma)^(alpha)
		ComponentWiseMultiply(N, col_vals, s_deltas)
		// Update one term in (-i)^(alpha) / (alpha!)
		cblas128.Scal(N, complex(0.0, -1.0/float64(alpha)), col_vals)

		// Copy the current values in.
		CopyColumn(sum_parts, alpha, col_vals)
	}

	bin_coeffs := cblas128.General{
		Rows:   num_bins,
		Cols:   M,
		Stride: num_bins,
		Data:   make([]complex128, num_bins*M),
	}
	for i := 0; i < N; i++ {
		row_for_s := cblas128.Vector{
			Inc:  sum_parts.Stride,
			Data: sum_parts.Data[i:],
		}
		// y += alpha x; Axpy(N, alpha, x, y)
		row_for_bin := cblas128.Vector{
			Inc:  bin_coeffs.Stride,
			Data: bin_coeffs.Data[bin_indices[i]:],
		}
		cblas128.Axpy(N, 1.0, row_for_s, row_for_bin)
	}

	result := make([]BinData, num_bins)
	curr_sigma := min_s + 0.5*bin_width
	for i := 0; i < num_bins; i++ {
		row_for_bin := cblas128.Vector{
			Inc:  bin_coeffs.Stride,
			Data: bin_coeffs.Data[i:],
		}
		result[i] = BinData{
			Tau:    tau,
			Sigma:  curr_sigma,
			Coeffs: row_for_bin,
		}
		curr_sigma += bin_width
	}
	return result
}

func PrintVector(N int, v cblas128.Vector) {
	fmt.Printf("[")
	index := 0
	for i := 0; i < N; i++ {
		fmt.Printf("%v", v.Data[index])
		fmt.Printf(" ")
		index += v.Inc
	}
	fmt.Println("]")
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

func checkCreateInitialData() {
	s := []float64{0.0, 1.0}
	min_s, max_s := 0.0, 1.0
	tau := 1.0
	actual_data := []complex128{3.0, 4.0}
	num_bins := 4
	M := 2
	bin_data := createInitialData(s, min_s, max_s, tau, actual_data, num_bins, M)
	for i := 0; i < num_bins; i++ {
		fmt.Printf("%v ", bin_data[i].Tau)
		fmt.Printf("%v ", bin_data[i].Sigma)
		PrintVector(M, bin_data[i].Coeffs)
	}
}

func timeRepeatAsBytes() {
	val := complex(4.2, 13.37)
	fmt.Println(RepeatAsBytes(val, 3))

	var start time.Time
	var elapsed_repeat time.Duration
	for i := 0; i < 1000; i++ {
		start = time.Now()
		RepeatAsBytes(val, 4096)
		elapsed_repeat += time.Since(start)
	}

	fmt.Println(elapsed_repeat)
}

func main() {
	// See: http://godoc.org/github.com/gonum/blas
	// sudo apt-get install libopenblas-dev
	// CGO_LDFLAGS="-L/usr/lib/libopenblas.so -lopenblas" go install github.com/gonum/blas/cgo
	timeRepeatAsBytes()
	checkCreateInitialData()
}
