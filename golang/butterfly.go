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

func main() {
	// See: http://godoc.org/github.com/gonum/blas
	// sudo apt-get install libopenblas-dev
	// CGO_LDFLAGS="-L/usr/lib/libopenblas.so -lopenblas" go install github.com/gonum/blas/cgo
	v := cblas128.Vector{Inc: 1, Data: []complex128{1, 1, 1}}
	fmt.Println("v has length:", cblas128.Nrm2(len(v.Data), v))

	blue_whale_call, err := LoadWhale()
	if err != nil {
		fmt.Printf("File error: %v\n", err)
	} else {
		fmt.Printf("Number of values: %v\n", len(*blue_whale_call))
	}
}
