package main

import "math"

type RData struct{
	Vx, Vy, V float64
	ctrl control
	Hrho float64
	rho float64
	dVx, dVy float64
	T,H,S float64
	theta, Theta float64
	dtheta float64
	bparam []BData
	alphaNV float64
	c_mang mang
	vcr float64
}
type BData struct{
	G0 []float64
	G1 []float64
}
type control struct {
	Fi7 float64
	DeltaV float64
	DeltaK float64
}
type mang struct{
	c0 []float64
	c1 []float64
	c2 []float64
	c3 []float64
	c4 []float64
}

func run(H0, V0, Vy0, theta0 float64, ctrl control){
	var(
		rd RData
	)
	rd.Hrho = H0
	rd.V = V0
	rd.rho = 1.225*(20000-rd.Hrho)/(20000+rd.Hrho)
	rd.Vy = Vy0
	rd.Theta = radians(90)
	if rd.V != 0 {
		rd.Theta = math.Asin(rd.Vy/rd.V)
	}
	rd.Vx = rd.V*math.Cos(rd.Theta)
	rd.theta = radians(theta0)
	rd.alphaNV = rd.theta - rd.Theta
	rd.c_mang = mangler(rd.alphaNV)
	rd.bparam = make([]BData, 360)

}

func degrees(arg float64) float64{
	d := arg/math.Pi*180
	return d
}

func radians(arg float64) float64{
	r := arg*math.Pi/180
	return r
}

func initial(dt RData) RData{
	for i:=0;i<len(dt.bparam);i++{

	}
	return dt
}

func mangler(alphaNV float64) mang{
	var(
		m mang
		mu []float64
		)
	m.c0 = make([]float64, N)
	mu = make([]float64,N)
	for i:=0;i<N;i++{
		mu[i] = 1-math.Pow(_r_[i], 2)
		m.c0[i] = 15*math.Pow(mu[i],2)*(1-math.Pow(mu[i],2))/8
		m.c1[i] = -15*math.Pi*(5-9*math.Pow(mu[i],2))*math.Pow((1- math.Pow(mu[i],2)),0.5)*math.Pow((1-math.Sin(alphaNV)/(1+math.Sin(alphaNV))),0.5) / 256
		m.c2[i] = math.Pow(-1, 0)*15/8*((mu[i]+2)*(9*math.Pow(mu[i],2)+4-6)/(-15)+3*mu[i]/(-5))*math.Pow((1-mu[i])/(1+mu[i]), 1)*math.Pow((1-math.Sin(alphaNV)/(1+math.Sin(alphaNV))),1)
		m.c3[i] = 45*math.Pi/256*math.Pow((1-math.Pow(mu[i],2)), 1.5)*math.Pow((1-math.Sin(alphaNV)/(1+math.Sin(alphaNV))),1.5)
		m.c4[i] = math.Pow(-1, 1)*15/8*((mu[i]+4)*(9*math.Pow(mu[i],2)+16-6)/(15*7)+3*mu[i]/(-5))*math.Pow((1-mu[i])/(1+mu[i]), 2)*math.Pow((1-math.Sin(alphaNV)/(1+math.Sin(alphaNV))),2)
	}
		return m
	}

func v0(T, V, Rho, AlphaNV float64, vcr float64) float64{
	S := math.Pi*math.Pow(R, 2)
	Ct:=2*T/(Rho*math.Pow(wR,2)*S)
	Vrel := V/wR*math.Sin(AlphaNV)
	Lambda := Vrel-vcr
	mu := V*math.Cos(AlphaNV)/wR
	vcr1 := Ct/(4*math.Sqrt(math.Pow(Lambda,2)+math.Pow(mu,2)))
	return vcr1
}
