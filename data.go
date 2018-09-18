package main

import (
	"fmt"
	"math"
	"time"
)

type RData struct {
	Vx, Vy, V    float64
	ctrl         control
	Hrho         float64
	rho          float64
	dVx, dVy     float64
	T, H, S      float64
	theta, Theta float64
	dtheta       float64
	bparam       []BData
	alphaNV      float64
	c_mang       mang
	vcr          float64
	D1 float64
	t1, t2 float64

}
type Forces struct {
	T,H,S float64
	Mx, My, Mz float64
}
type BData struct {
	psi, dpsi  float64
	vind []float64
	fi   []float64
	Vind float64
	_u1_ []float64
	_v1_ []float64
	beta, beta1, beta2 float64
	mu float64
	_vbeta_ []float64
	F []float64
	_w1_ []float64
	W1 []float64
	al []float64
	mach []float64
	Cy, Cx []float64
	dy, dx []float64
	dt, dq, dqp,dqi []float64
	dmx, dmy, dmz []float64
	dh, ds []float64
	dfr []float64
	T, H, S, Mx, My, Mz float64

}
type control struct {
	Fi7    float64
	DeltaV float64
	DeltaK float64
}
type mang struct {
	c0 []float64
	c1 []float64
	c2 []float64
	c3 []float64
	c4 []float64
}

func run(H0, V0, Vy0, theta0 float64, ctrl control) {
	var (
		rd RData
		his []BData
	)
	rd.Hrho = H0
	rd.V = V0
	rd.rho = 1.225 * (20000 - rd.Hrho) / (20000 + rd.Hrho)
	rd.Vy = Vy0
	rd.Theta = radians(90)
	if rd.V != 0 {
		rd.Theta = math.Asin(rd.Vy / rd.V)
	}
	rd.ctrl = ctrl
	rd.Vx = rd.V * math.Cos(rd.Theta)
	rd.theta = radians(theta0)
	rd.alphaNV = rd.Theta - rd.theta
	rd.c_mang = mangler(rd.alphaNV)
	rd.bparam = make([]BData, z)
	his = make([]BData, z)
	rd.t1 = radians(rd.D1*(rd.ctrl.DeltaK+k*rd.ctrl.DeltaV))
	rd.t2 = radians(rd.D1*(-rd.ctrl.DeltaV+k*rd.ctrl.DeltaK))
	T:=m0*g
	T1:=0.0
	start := time.Now()
	vh := math.Pow((T / (2 * rd.rho * math.Pi * math.Pow(R, 2))), 0.5)
	Vind := v0(T, rd.V, rd.rho, rd.alphaNV, vh/wR)
	for k:=0;k<z;k++{
		ps := float64(k*360/z)
		rd.bparam[k] = init0(rd, ps, Vind)
		T1 = T1 + rd.bparam[k].T
		his[k] = rd.bparam[k]
		//fmt.Printf("%f\t", rd.bparam[k].T)
	}
	//fmt.Println("")
	elapsed := time.Since(start)
	fmt.Printf("Init time %s\n", elapsed)
	
//fmt.Printf("%f\n", T1)
	for j:=1;j<360;j++{
		T = (T+T1)/2
		T1 = 0
		vh = math.Pow((T / (2 * rd.rho * math.Pi * math.Pow(R, 2))), 0.5)
		Vind = v0(T, rd.V, rd.rho, rd.alphaNV, vh/wR)
		for k:=0;k<z;k++{
			rd.bparam[k] = calc(rd.bparam[k], his[k], rd, Vind, Nps*float64(j)+(360/float64(z))*float64(k))
			T1 = T1 + rd.bparam[k].T
			his[k] = rd.bparam[k]
			//fmt.Printf("%f\t", rd.bparam[k].T)
		}
		//fmt.Println("")
		//fmt.Printf("T=%f\tT1=%f\tpsi=%d\tbeta=%f\tVind=%f\n", T,T1, j, degrees(rd.bparam[0].beta), Vind)
	}
	elapsed = time.Since(start)
	fmt.Printf("Round time %s\n", elapsed)
}

func degrees(arg float64) float64 {
	d := arg / math.Pi * 180
	return d
}

func radians(arg float64) float64 {
	r := arg * math.Pi / 180
	return r
}

func init0(dt RData, psi, Vind float64) BData {
	var (
		ret BData
				i int
	)

	ret.psi = radians(psi)

	//fmt.Printf("vh = %f\n",vh)
	ret.Vind = Vind
	//fmt.Printf("%f\n", ret.Vind)
	ret.mu = dt.V*math.Cos(dt.alphaNV)/wR
	ret.vind = make([]float64, N)
	ret.fi = make([]float64, N)
	ret._u1_ = make([]float64, N)
	ret._vbeta_ = make([]float64, N)
	ret._v1_ = make([]float64, N)
	ret.F = make([]float64, N)
	ret._w1_ = make([]float64, N)
	ret.al = make([]float64, N)
	ret.mach = make([]float64, N)
	ret.Cy = make([]float64, N)
	ret.Cx = make([]float64, N)
	ret.dy = make([]float64, N)
	ret.dx = make([]float64, N)
	ret.W1 = make([]float64, N)
	ret.dt = make([]float64, N)
	ret.dq = make([]float64, N)
	ret.dqp = make([]float64, N)
	ret.dqi = make([]float64, N)
	ret.dmx = make([]float64, N)
	ret.dmy = make([]float64, N)
	ret.dmz = make([]float64, N)
	ret.dfr = make([]float64, N)
	ret.dh = make([]float64, N)
	ret.ds = make([]float64, N)
	for i = 0; i < N; i++ {
		ret.vind[i] = -4 * ret.Vind * (0.5*dt.c_mang.c0[i] - dt.c_mang.c1[i]*math.Cos(ret.psi) - dt.c_mang.c2[i]*math.Cos(2*ret.psi) - dt.c_mang.c3[i]*math.Cos(3*ret.psi) - dt.c_mang.c4[i]*math.Cos(4*ret.psi))
		ret.fi[i] = radians(dt.ctrl.Fi7)+radians(dFi[i])-k*ret.beta+dt.t1*math.Cos(ret.psi)+dt.t2*math.Sin(ret.psi)
		ret._u1_[i] = _r_[i]+ret.mu*math.Sin(ret.psi)
		ret._vbeta_[i] = -_r_[i]*ret.beta1-ret.mu*math.Cos(dt.alphaNV)*math.Sin(ret.beta)*math.Cos(ret.psi)
		ret._v1_[i] = dt.V*math.Sin(dt.alphaNV)/wR+ret.vind[i]+ret._vbeta_[i]
		ret._w1_[i] = math.Sqrt(ret._u1_[i]*ret._u1_[i]+ret._v1_[i]*ret._v1_[i])
		ret.W1[i] = ret._w1_[i]*wR
		ret.F[i] = -math.Atan2(ret._v1_[i], ret._u1_[i])
		ret.al[i] = degrees(ret.fi[i] - ret.F[i])
		ret.mach[i] = ret.W1[i]/340
		ret.Cy[i] = tCy(r[i], ret.al[i], ret.mach[i])
		ret.Cx[i] = tCx(r[i], ret.al[i], ret.mach[i])
		ret.dy[i] = ret.Cy[i]*dt.rho*ret.W1[i]*ret.W1[i]*b/2
		ret.dx[i] = ret.Cx[i]*dt.rho*ret.W1[i]*ret.W1[i]*b/2
		ret.dt[i] = (ret.dy[i]*math.Cos(ret.F[i])-ret.dx[i]*math.Sin(ret.F[i]))*math.Cos(ret.beta)
		ret.dqp[i] = ret.dx[i]*math.Cos(ret.F[i])
		ret.dqi[i] = ret.dy[i]*math.Sin(ret.F[i])
		ret.dq[i] = ret.dqp[i]+ret.dqi[i]
		ret.dmx[i] = r[i]*ret.dt[i]*math.Sin(ret.psi)
		ret.dmy[i] = ret.dq[i]*r[i]
		ret.dmz[i] = -ret.dt[i]*r[i]*math.Cos(ret.psi)
		ret.dfr[i] = ret.dt[i]*math.Tan(ret.beta)
		ret.dh[i] = ret.dq[i]*math.Sin(ret.psi)+ret.dfr[i]*math.Cos(ret.psi)
		ret.ds[i] = -ret.dq[i]*math.Cos(ret.psi)+ret.dfr[i]*math.Sin(ret.psi)
	}
	/*str := strconv.Itoa(int(psi))+" "
	plot(r, ret.vind, str+"1.png","Радиус, м", "Индуктивная скорость, м/с")
	plot(r, ret.Cx, str+"Cx.png","Радиус, м", "Cx")
	plot(r, ret.Cy, str+"Cy.png","Радиус, м", "Cy")
	plot(r, ret.al, str+"Al.png","Радиус, м", "Угол атаки")
	plot(r, ret.mach, str+"M.png","Радиус, м", "Число Маха")
	plot(r, ret.dy, str+"dy.png","Радиус, м", "Воздушная нагрузка dy, Н/м")
	plot(r, ret.dx, str+"dx.png","Радиус, м", "Воздушная нагрузка dx, Н/м")
	plot(r, ret.dt, str+"dt.png","Радиус, м", "Воздушная нагрузка dt, Н/м")
	plot(r, ret.dh, str+"dh.png","Радиус, м", "Воздушная нагрузка dh, Н/м")
	plot(r, ret.ds, str+"ds.png","Радиус, м", "Воздушная нагрузка ds, Н/м")
	plot(r, ret.dq, str+"dq.png","Радиус, м", "Воздушная нагрузка dq, Н/м")
	plot(r, ret.dfr, str+"dfr.png","Радиус, м", "Воздушная нагрузка dfr, Н/м")
	plot(r, ret.dmx, str+"dmx.png","Радиус, м", "Погонный момент dmx, Нм")
	plot(r, ret.dmy, str+"dmy.png","Радиус, м", "Погонный момент dmy, Нм")
	plot(r, ret.dmz, str+"dmz.png","Радиус, м", "Погонный момент dmz, Нм")
	*/
	ret.T = 0
	ret.H = 0
	ret.S = 0
	ret.Mx = 0
	ret.My = 0
	ret.Mz = 0
	for i=0;i<N-1;i++{
		ret.T = ret.T + (ret.dt[i]+ret.dt[i+1])*dr[i]/2
		ret.H = ret.H + (ret.dh[i]+ret.dh[i+1])*dr[i]/2
		ret.S = ret.S + (ret.ds[i]+ret.ds[i+1])*dr[i]/2
		ret.Mx = ret.Mx + (ret.dmx[i]+ret.dmx[i+1])*dr[i]/2
		ret.My = ret.My + (ret.dmy[i]+ret.dmy[i+1])*dr[i]/2
		ret.Mz = ret.Mz + (ret.dmz[i]+ret.dmz[i+1])*dr[i]/2
	}
/*	fmt.Printf("T=%f\n", ret.T)
	fmt.Printf("H=%f\n", ret.H)
	fmt.Printf("S=%f\n", ret.S)
	fmt.Printf("Mx=%f\n", ret.Mx)
	fmt.Printf("My=%f\n", ret.My)
	fmt.Printf("Mz=%f\n", ret.Mz)*/
	ST:=stat_mom(r,ret.dt)
	//fmt.Printf("%f\n", ST)
	A1:=(Ihh*math.Cos(ret.beta)-Lhh*Shh)*math.Pow(Omega,2)*math.Sin(ret.beta)
	ret.beta2 = (ST-A1-g*Shh)/(Ihh*math.Pow(Omega,2))
	//fmt.Printf("d2beta=%f\n", ret.beta2)
	ret.beta1 = ret.beta1 + ret.beta2
	//fmt.Printf("dbeta=%f\n", ret.beta1)
	ret.beta = ret.beta + ret.beta1
	//fmt.Printf("beta=%f\n", degrees(ret.beta))
	return ret
}

func mangler(alphaNV float64) mang {
	var (
		m  mang
		mu []float64
	)
	m.c0 = make([]float64, N)
	m.c1 = make([]float64, N)
	m.c2 = make([]float64, N)
	m.c3 = make([]float64, N)
	m.c4 = make([]float64, N)
	mu = make([]float64, N)
	for i := 0; i < N; i++ {
		mu[i] = math.Sqrt(1 - math.Pow(_r_[i], 2))
		//fmt.Printf("%f\t", _r_[i])
		//fmt.Printf("%f\n", mu[i])
		m.c0[i] = 15 * mu[i] * (1 - math.Pow(mu[i], 2)) / 8
		m.c1[i] = -15 * math.Pi * (5 - 9*math.Pow(mu[i], 2)) * math.Pow((1-math.Pow(mu[i], 2)), 0.5) * math.Pow((1-math.Sin(alphaNV))/(1+math.Sin(alphaNV)), 0.5) / 256
		m.c1[i] = (1 - math.Sin(alphaNV)) / (1 + math.Sin(alphaNV))
		m.c2[i] = math.Pow(-1, 0) * 15 / 8 * ((mu[i]+2)*(9*math.Pow(mu[i], 2)+4-6)/(-15) + 3*mu[i]/(-5)) * math.Pow((1-mu[i])/(1+mu[i]), 1) * math.Pow((1-math.Sin(alphaNV))/(1+math.Sin(alphaNV)), 1)
		m.c3[i] = 45 * math.Pi / 256 * math.Pow((1-math.Pow(mu[i], 2)), 1.5) * math.Pow((1-math.Sin(alphaNV))/(1+math.Sin(alphaNV)), 1.5)
		m.c4[i] = math.Pow(-1, 1) * 15 / 8 * ((mu[i]+4)*(9*math.Pow(mu[i], 2)+16-6)/(15*7) + 3*mu[i]/(-5)) * math.Pow((1-mu[i])/(1+mu[i]), 2) * math.Pow((1-math.Sin(alphaNV))/(1+math.Sin(alphaNV)), 2)
	}
	return m
}

func v0(T, V, Rho, AlphaNV float64, vcr float64) float64 {
	S := math.Pi * math.Pow(R, 2)
	Ct := 2 * T / (Rho * math.Pow(wR, 2) * S)

	Vrel := V / wR * math.Sin(AlphaNV)
	Lambda := Vrel - vcr
	//fmt.Printf("Lambda = %f\n",Lambda)
	mu := V * math.Cos(AlphaNV) / wR
	//fmt.Printf("mu = %f\n",mu)
	vcr1 := Ct / (4 * math.Sqrt(math.Pow(Lambda, 2)+math.Pow(mu, 2)))
	//fmt.Printf("vcrl = %f\n",vcr1*wR)
	return vcr1
}

func stat_mom(r1 []float64, m1 []float64) float64{
	ST:= 0.0
	for i:=0; i<len(r1)-1;i++{
		ST = ST + ((m1[i]*(r1[i]-Lhh))+(m1[i+1]*(r1[i+1]-Lhh)))/2*(r1[i+1]-r1[i])
		}
	return ST
}

func calc(ret, old BData, dt RData, Vind, psi float64) BData{
	var (
		i int
	)
	ret.psi = radians(psi)
	ret.dpsi = ret.psi - old.psi
	//fmt.Printf("d2beta=%f\n", old.beta2)
	ret.beta1 = ret.beta1 + old.beta2*ret.dpsi
	//fmt.Printf("dbeta=%f\n", old.beta1)
	ret.beta = ret.beta + old.beta1*ret.dpsi
	//fmt.Printf("beta=%f\n", degrees(ret.beta))
	//math.Pow((T / (2 * dt.rho * math.Pi * math.Pow(R, 2))), 0.5)
	//fmt.Printf("vh = %f\n",vh)
	ret.Vind = Vind
	//fmt.Printf("%f\n", ret.Vind)
	ret.mu = dt.V*math.Cos(dt.alphaNV)/wR
	//fmt.Printf("Vind=%f\n", ret.Vind)
	for i = 0; i < N; i++ {
		ret.vind[i] = -4 * ret.Vind * (0.5*dt.c_mang.c0[i] - dt.c_mang.c1[i]*math.Cos(ret.psi) - dt.c_mang.c2[i]*math.Cos(2*ret.psi) - dt.c_mang.c3[i]*math.Cos(3*ret.psi) - dt.c_mang.c4[i]*math.Cos(4*ret.psi))
		ret.fi[i] = radians(dt.ctrl.Fi7)+radians(dFi[i])-k*ret.beta+dt.t1*math.Cos(ret.psi)+dt.t2*math.Sin(ret.psi)
		ret._u1_[i] = _r_[i]+ret.mu*math.Sin(ret.psi)
		ret._vbeta_[i] = -_r_[i]*ret.beta1-ret.mu*math.Cos(dt.alphaNV)*math.Sin(ret.beta)*math.Cos(ret.psi)
		ret._v1_[i] = dt.V*math.Sin(dt.alphaNV)/wR+ret.vind[i]+ret._vbeta_[i]
		ret._w1_[i] = math.Sqrt(ret._u1_[i]*ret._u1_[i]+ret._v1_[i]*ret._v1_[i])
		ret.W1[i] = ret._w1_[i]*wR
		ret.F[i] = -math.Atan2(ret._v1_[i], ret._u1_[i])
		ret.al[i] = degrees(ret.fi[i] - ret.F[i])
		ret.mach[i] = ret.W1[i]/340
		//fmt.Printf("%f\t%f\t%f\n", r[i], ret.al[i], ret.mach[i])
		ret.Cy[i] = tCy(r[i], ret.al[i], ret.mach[i])

		ret.Cx[i] = tCx(r[i], ret.al[i], ret.mach[i])

		ret.dy[i] = ret.Cy[i]*dt.rho*ret.W1[i]*ret.W1[i]*b/2
		ret.dx[i] = ret.Cx[i]*dt.rho*ret.W1[i]*ret.W1[i]*b/2
		ret.dt[i] = (ret.dy[i]*math.Cos(ret.F[i])-ret.dx[i]*math.Sin(ret.F[i]))*math.Cos(ret.beta)
		ret.dqp[i] = ret.dx[i]*math.Cos(ret.F[i])
		ret.dqi[i] = ret.dy[i]*math.Sin(ret.F[i])
		ret.dq[i] = ret.dqp[i]+ret.dqi[i]
		ret.dmx[i] = r[i]*ret.dt[i]*math.Sin(ret.psi)
		ret.dmy[i] = ret.dq[i]*r[i]
		ret.dmz[i] = -ret.dt[i]*r[i]*math.Cos(ret.psi)
		ret.dfr[i] = ret.dt[i]*math.Tan(ret.beta)
		ret.dh[i] = ret.dq[i]*math.Sin(ret.psi)+ret.dfr[i]*math.Cos(ret.psi)
		ret.ds[i] = -ret.dq[i]*math.Cos(ret.psi)+ret.dfr[i]*math.Sin(ret.psi)
	}
	/*str := strconv.Itoa(int(psi))+" "
	plot(r, ret.vind, str+"1.png","Радиус, м", "Индуктивная скорость, м/с")
	plot(r, ret.Cx, str+"Cx.png","Радиус, м", "Cx")
	plot(r, ret.Cy, str+"Cy.png","Радиус, м", "Cy")
	plot(r, ret.al, str+"Al.png","Радиус, м", "Угол атаки")
	plot(r, ret.mach, str+"M.png","Радиус, м", "Число Маха")
	plot(r, ret.dy, str+"dy.png","Радиус, м", "Воздушная нагрузка dy, Н/м")
	plot(r, ret.dx, str+"dx.png","Радиус, м", "Воздушная нагрузка dx, Н/м")
	plot(r, ret.dt, str+"dt.png","Радиус, м", "Воздушная нагрузка dt, Н/м")
	plot(r, ret.dh, str+"dh.png","Радиус, м", "Воздушная нагрузка dh, Н/м")
	plot(r, ret.ds, str+"ds.png","Радиус, м", "Воздушная нагрузка ds, Н/м")
	plot(r, ret.dq, str+"dq.png","Радиус, м", "Воздушная нагрузка dq, Н/м")
	plot(r, ret.dfr, str+"dfr.png","Радиус, м", "Воздушная нагрузка dfr, Н/м")
	plot(r, ret.dmx, str+"dmx.png","Радиус, м", "Погонный момент dmx, Нм")
	plot(r, ret.dmy, str+"dmy.png","Радиус, м", "Погонный момент dmy, Нм")
	plot(r, ret.dmz, str+"dmz.png","Радиус, м", "Погонный момент dmz, Нм")
	*/
	ret.T = 0
	ret.H = 0
	ret.S = 0
	ret.Mx = 0
	ret.My = 0
	ret.Mz = 0
	for i=0;i<N-1;i++{
		ret.T = ret.T + (ret.dt[i]+ret.dt[i+1])*dr[i]/2
		ret.H = ret.H + (ret.dh[i]+ret.dh[i+1])*dr[i]/2
		ret.S = ret.S + (ret.ds[i]+ret.ds[i+1])*dr[i]/2
		ret.Mx = ret.Mx + (ret.dmx[i]+ret.dmx[i+1])*dr[i]/2
		ret.My = ret.My + (ret.dmy[i]+ret.dmy[i+1])*dr[i]/2
		ret.Mz = ret.Mz + (ret.dmz[i]+ret.dmz[i+1])*dr[i]/2
	}
	/*
	fmt.Printf("psi=%f\t", degrees(ret.psi))
	fmt.Printf("T=%f\t", ret.T)
	fmt.Printf("H=%f\t", ret.H)
	fmt.Printf("S=%f\t", ret.S)
	fmt.Printf("Mx=%f\t", ret.Mx)
	fmt.Printf("My=%f\t", ret.My)
	fmt.Printf("Mz=%f\n", ret.Mz)
*/
	ST:=stat_mom(r,ret.dt)
	//fmt.Printf("%f\n", ST)
	A1:=(Ihh*math.Cos(ret.beta)-Lhh*Shh)*math.Pow(Omega,2)*math.Sin(ret.beta)
	ret.beta2 = (ST-A1-g*Shh)/(Ihh*math.Pow(Omega,2))
	return ret
}