package main

import (
	"fmt"
	"math"
	"strconv"

	"github.com/tealeg/xlsx"
)

var (
	m0       float64
	R        float64
	r0       float64
	b        float64
	RPM      float64
	z        int
	k        float64
	dPsi     float64
	Lhh      float64
	wR       float64
	N        = 300
	r        []float64
	_r_      []float64
	dFi      []float64
	AirTwist []air_twist
)

type air_twist struct {
	r    float64
	foil Foil
}

func main() {
	read_data("test.xlsx")
	//foil := Read_foil("NACA 63012")
	/*	for i := 0; i <len(foil.Mach);i++ {
		fmt.Println(foil.Mach[i])
		for k := 0; k<len(foil.Data[i].Alpha); k++{
			fmt.Printf("%f\t", foil.Data[i].Alpha[k])
			fmt.Printf("%f\t",foil.Data[i].Cy[k])
			fmt.Printf("%f\n",foil.Data[i].Cx[k])
		}
	}*/
	var ctrl control
	ctrl.Fi7 = 8
	run(0,0,0,0, ctrl)
	//output()
}

func read_data(filename string) {
	var (
		i int
	)
	xlFile, err := xlsx.OpenFile(filename)
	if err != nil {
		fmt.Print("Error")
	}
	sheet := xlFile.Sheet["HeliData"]
	row := sheet.Rows
	mp := make(map[string]float64)
	for i := 0; i < len(row); i += 1 {
		mp[row[i].Cells[0].Value], err = strconv.ParseFloat(row[i].Cells[1].Value, 64)
	}
	m0 = mp["m0"]
	R = mp["R"]
	r0 = mp["r0"]
	b = mp["b"]
	RPM = mp["RPM"]
	z = int(mp["z"])
	k = mp["k"]
	dPsi = mp["dPsi"]
	Lhh = mp["Lhh"]
	wR = math.Pi * RPM * R / 30
	r = make([]float64, N)
	_r_ = make([]float64, N)
	dFi = make([]float64, N)
	r[0] = r0
	_r_[0] = r[0] / R
	for i = 1; i < N; i++ {
		r[i] = r[i-1] + (R-r0)/(float64(N)-1)
		_r_[i] = r[i] / R
	}
	_r_[N-1] = 1
	sheet = xlFile.Sheet["Twist"]
	row = sheet.Rows
	Xt := make([]float64, len(row))
	Yt := make([]float64, len(row))
	for i := 1; i < len(row); i++ {
		Xt[i], err = strconv.ParseFloat(row[i].Cells[0].Value, 64)
		Yt[i], err = strconv.ParseFloat(row[i].Cells[1].Value, 64)
	}
	for i = 0; i < N; i++ {
		dFi[i] = aprox(r[i], Xt, Yt)
	}
	sheet = xlFile.Sheet["Foils"]
	row = sheet.Rows
	fls := make(map[float64]string)
	for i := 1; i < len(row); i += 1 {
		arg, err1 := strconv.ParseFloat(row[i].Cells[0].Value, 64)
		if err1 != nil {
			fmt.Printf(err1.Error())
		}
		fls[arg] = row[i].Cells[1].Value
	}
	AirTwist = make([]air_twist, len(fls))
	i = 0
	for key, val := range fls {
		AirTwist[i].r = key
		AirTwist[i].foil = Read_foil(val)
		i++
	}

}

func aprox(arg float64, X []float64, Y []float64) float64 {
	var (
		i      int
		j      int
		x      float64
		x1, x2 float64
		y1, y2 float64
		y      float64
	)
	i = 0
	j = 0
	for i = 0; i < len(X)-1; i++ {
		if X[j] <= arg {
			j = j + 1
		}
	}
	x = arg
	x1 = X[j-1]
	x2 = X[j]
	y1 = Y[j-1]
	y2 = Y[j]
	y = (x-x1)*(y2-y1)/(x2-x1) + y1
	return y
}

func output() {
	for i := 0; i < N; i++ {
		fmt.Printf("%d\t", i)
		fmt.Printf("%f\t", r[i])
		fmt.Printf("%f\t", _r_[i])
		fmt.Printf("%f\n", dFi[i])
	}
}

func tCy(r, alpha, m float64) float64 {
	var (
		i, j   int
		x1, x2 float64
		y1, y2 float64
		y      float64
	)
	j = 0

	for i = 0; i < len(AirTwist); i++ {
		if AirTwist[i].r <= r {
			j++
		}
	}

	x1 = AirTwist[j-1].r
	x2 = AirTwist[j].r

	y1 = FindCy(alpha, m, AirTwist[j-1].foil)
	y2 = FindCy(alpha, m, AirTwist[j].foil)

	y = (r-x1)*(y2-y1)/(x2-x1) + y1
	return y
}

func tCx(r, alpha, m float64) float64 {
	var (
		i, j   int
		x1, x2 float64
		y1, y2 float64
		y      float64
	)
	j = 0

	for i = 0; i < len(AirTwist); i++ {
		if AirTwist[i].r <= r {
			j++
		}
	}

	x1 = AirTwist[j-1].r
	x2 = AirTwist[j].r

	y1 = FindCx(alpha, m, AirTwist[j-1].foil)
	y2 = FindCx(alpha, m, AirTwist[j].foil)

	y = (r-x1)*(y2-y1)/(x2-x1) + y1
	return y
}
