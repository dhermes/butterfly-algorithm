package main

import (
	"encoding/json"
	"fmt"
	"io/ioutil"
	"math"
	"math/cmplx"
	"path"
	"runtime"
	"time"

	"github.com/gonum/blas/cblas128"
	"github.com/mjibson/go-dsp/fft"
)

const minusI complex128 = complex(0, -1)

func LoadWhale() (*[]complex128, error) {
	// H/T: andrewbrookins.com/tech/golang-get-directory-of-the-current-file/
	_, filename, _, _ := runtime.Caller(1)
	whale_path := path.Join(path.Dir(filename), "bluewhale.json")

	contents, err := ioutil.ReadFile(whale_path)
	if err != nil {
		return nil, err
	}

	var float64_call []float64
	json.Unmarshal(contents, &float64_call)

	blue_whale_call := make([]complex128, len(float64_call))
	for i, value := range float64_call {
		blue_whale_call[i] = complex(value, 0.0)
	}
	return &blue_whale_call, nil
}

func dftKernel(t, s complex128) complex128 {
	return cmplx.Exp(minusI * t * s)
}

func dftKernelVectorized(t complex128, s []complex128) []complex128 {
	N := len(s)
	result := make([]complex128, N)
	v := cblas128.Vector{Inc: 1, Data: s}
	w := cblas128.Vector{Inc: 1, Data: result}
	cblas128.Axpy(N, t*minusI, v, w)
	// TODO(dhermes): Determine if there is a way to vectorize this
	//                operation.
	for i, value := range result {
		result[i] = cmplx.Exp(value)
	}
	return result
}

func computeFHat(f, t, s []complex128) []complex128 {
	N := len(f)
	f_hat := make([]complex128, N)
	f_vec := cblas128.Vector{Inc: 1, Data: f}
	for k := 0; k < N; k++ {
		v := cblas128.Vector{Inc: 1, Data: dftKernelVectorized(t[k], s)}
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

func compareNaive() {
	blue_whale_call, err := LoadWhale()
	if err != nil {
		fmt.Printf("File error: %v\n", err)
	}

	N := len(*blue_whale_call)
	t, s := getDFTData(N)

	start := time.Now()
	dft_whale_call := fft.FFT(*blue_whale_call)
	elapsed_fft := time.Since(start)

	start = time.Now()
	computed_f_hat := computeFHat(*blue_whale_call, t, s)
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

func main() {
	// See: http://godoc.org/github.com/gonum/blas
	// sudo apt-get install libopenblas-dev
	// CGO_LDFLAGS="-L/usr/lib/libopenblas.so -lopenblas" go install github.com/gonum/blas/cgo
	compareNaive()
}
