package main

import (
	"encoding/json"
	"fmt"
	"io/ioutil"
	"math/cmplx"
	"path"
	"runtime"

	"github.com/gonum/blas/cblas128"
)

const minusI complex128 = complex(0, -1)

func LoadWhale() (*[]float64, error) {
	// H/T: andrewbrookins.com/tech/golang-get-directory-of-the-current-file/
	_, filename, _, _ := runtime.Caller(1)
	whale_path := path.Join(path.Dir(filename), "bluewhale.json")

	contents, err := ioutil.ReadFile(whale_path)
	if err != nil {
		return nil, err
	}

	var blue_whale_call []float64
	json.Unmarshal(contents, &blue_whale_call)
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
		fmt.Println("v.Data:", v.Data)
		f_hat[k] = cblas128.Dotu(N, v, f_vec)
	}
	return f_hat
}

func main() {
	// See: http://godoc.org/github.com/gonum/blas
	// sudo apt-get install libopenblas-dev
	// CGO_LDFLAGS="-L/usr/lib/libopenblas.so -lopenblas" go install github.com/gonum/blas/cgo
	v := cblas128.Vector{Inc: 1, Data: []complex128{1, 1, 1}}
	fmt.Println("v has length:", cblas128.Nrm2(len(v.Data), v))

	t := []complex128{1.0, 2.0}
	s := []complex128{3.0, 4.0}
	f := []complex128{5.0, 6.0}
	fmt.Println("t:", t)
	fmt.Println("s:", s)
	fmt.Println("f:", f)
	fmt.Println("fhat:", computeFHat(f, t, s))

	blue_whale_call, err := LoadWhale()
	if err != nil {
		fmt.Printf("File error: %v\n", err)
	} else {
		fmt.Printf("Number of values: %v\n", len(*blue_whale_call))
	}
}
