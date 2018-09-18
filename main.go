package main

import (
	"fmt"
	"math"
	"sort"
	"strconv"

	"github.com/tealeg/xlsx"
)

const g  = 9.81

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
	N        = 100
	Nps 	=	 1.0
	dTime float64
	Omega float64
	r        []float64
	_r_      []float64
	dr []float64
	dFi      []float64
	AirTwist []air_twist
	EIx, EIy []float64
	mass []float64
	rm []float64
	Shh, Ihh float64
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
	Omega = wR/R
	dTime = 1/(Omega/(2*math.Pi)*(360/Nps))
	//fmt.Printf("%f\n", dTime)
	dr = make([]float64, N)
	r = make([]float64, N)
	_r_ = make([]float64, N)
	dFi = make([]float64, N)
	mass = make([]float64, N)
	EIx = make([]float64, N)
	r[0] = r0
	_r_[0] = r[0] / R
	for i = 1; i < N; i++ {
		r[i] = r[i-1] + (R-r0)/(float64(N)-1)
		_r_[i] = r[i] / R
		dr[i-1] = r[i]-r[i-1]
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
	var keys []float64
	for k:=range fls{
		keys = append(keys,k)
	}

		sort.Float64s(keys)
	//fmt.Printf("%v\n", keys)
	for _, k := range keys {
		AirTwist[i].r = k
		AirTwist[i].foil = Read_foil(fls[k])
		AirTwist[i].foil.name = fls[k]
		//fmt.Printf("%f\t", k)
		//fmt.Println(fls[k])
		i++
	}
	sheet = xlFile.Sheet["Mass"]
	row = sheet.Rows
	rm = make([]float64, len(row))

	for i := 1; i < len(row); i++ {
		rm[i], err = strconv.ParseFloat(row[i].Cells[0].Value, 64)
		mass[i], err = strconv.ParseFloat(row[i].Cells[1].Value, 64)
		EIx[i], err = strconv.ParseFloat(row[i].Cells[2].Value, 64)
		//fmt.Printf("%d\t%f\t%f\t%f\n",i,  rm[i], mass[i], EIx[i])
	}

	Shh = 0
	Ihh = 0
	for i:=0; i<len(rm)-1;i++{
		Shh = Shh + ((mass[i]*rm[i])+(mass[i+1]*rm[i+1]))/2*(rm[i+1]-rm[i])
		Ihh = Ihh + ((mass[i]*rm[i]*rm[i])+(mass[i+1]*rm[i+1]*rm[i+1]))/2*(rm[i+1]-rm[i])
	}
	 //fmt.Printf("%f\t%f\n", Shh, Ihh)
	//plot(rm, EIx, "EIx.png","Радиус, м", "Погонная жесткость, Нм")
	//plot(rm, mass, "mass.png","Радиус, м", "Погонная масса, кг/м")

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
	if arg>X[0]{
		x1 = X[j-1]
		y1 = Y[j-1]
		}
	if arg<=X[0]{
		x1 = X[0]
		y1 = Y[0]}
	x2 = X[j]

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

	for i = 1; i < len(AirTwist); i++ {
		j++
		if AirTwist[i].r <= r {
			break
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
		j++
		if AirTwist[i].r <= r {
			break
		}
	}

	x1 = AirTwist[j-1].r
	x2 = AirTwist[j].r

	y1 = FindCx(alpha, m, AirTwist[j-1].foil)
	y2 = FindCx(alpha, m, AirTwist[j].foil)

	y = (r-x1)*(y2-y1)/(x2-x1) + y1
	return y
}
