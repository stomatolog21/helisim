package main

import (
	"io/ioutil"
	"log"
	"fmt"
	"os"
	"bufio"
	"strings"
	"strconv"
	"sort"
)


type Foil struct {
	name string
	Mach []float64
	Data []foildata
}
type foildata struct {
	Alpha []float64
	Cy []float64
	Cx []float64
}

func Read_foil(name string) Foil {
	var(
		dir string
		foil Foil
		M []float64
		line string
		cy  []float64
		cx  []float64
		alpha  []float64
	)
	dir = ".\\"+name+"\\"
	files, err := ioutil.ReadDir(dir)
	if err != nil{
		log.Fatal(err)
	}

	file_name := make([]string, len(files))
	M = make([]float64, len(files))
	for i,f := range files{
		file_name[i] = dir + f.Name()
		//fmt.Println(file_name[i])
		file, err1 := os.Open(file_name[i])
		defer file.Close()
		if err1 != nil {
			fmt.Printf("Error")
		}
		reader := bufio.NewReader(file)
		for{
			line,err = reader.ReadString('\n')
			if strings.Contains(line, "Mach = "){
				runes :=[]rune(line)
				M[i], err = strconv.ParseFloat(string(runes[10:15]), 64)
			}
			if err != nil{
				break
			}
		}
		}
	type FLS struct{
		Mach float64
		File_Name string
	}
	fls := make([]FLS, len(M))
	for i := 0; i< len(M); i++ {
		fls[i].File_Name = file_name[i]
		fls[i].Mach = M[i]
	}
	sort.Slice(fls, func(i,j int) bool{return fls[i].Mach<fls[j].Mach})
	foil.Mach = make([]float64, len(M))
	for i:=0; i<len(M); i++ {
		foil.Mach[i] = fls[i].Mach
		}

	foil.Data = make([]foildata, len(M))

	for i := 0; i <len(fls); i++{
		file, err1 := os.Open(fls[i].File_Name)
		defer file.Close()
		if err1 != nil {
			fmt.Printf("Error")
		}
		reader := bufio.NewReader(file)
		for{
			line,err = reader.ReadString('\n')
			if strings.Contains(line, "-------"){
				j := 1
				for{

					line,err = reader.ReadString('\n')
					if line == "\r\n"{
						break
					}
					runes :=[]rune(line)
					a, err1 := strconv.ParseFloat(strings.Replace(string(runes[0:8]), " ", "", -1),64)
					b, err2 := strconv.ParseFloat(strings.Replace(string(runes[10:17]), " ", "", -1),64)
					c, err3 := strconv.ParseFloat(strings.Replace(string(runes[20:27]), " ", "", -1),64)
					alpha = append(alpha, 0)
					cy = append(cy,0)
					cx = append(cx, 0)
					alpha[len(alpha)-1] = a
					cy[len(cy)-1] = b
					cx[len(cx)-1] = c
					runes = runes[:0]
					if err1 != nil{
						fmt.Print("a ")
						fmt.Printf("%d ", j)
						fmt.Println(err1.Error())

					}
					if err2 != nil{
						fmt.Print("b ")
						fmt.Println(err2.Error())

					}
					if err3 != nil{
						fmt.Print("c ")
						fmt.Println(err3.Error())

					}
					if err != nil{
						break
					}
					j++
					}





			}

			if err != nil{
				break
			}
		}

		foil.Data[i].Alpha = make([]float64, len(alpha))
		foil.Data[i].Cy = make([]float64, len(alpha))
		foil.Data[i].Cx = make([]float64, len(alpha))

		for q := 0; q<len(alpha); q++{
			foil.Data[i].Alpha[q] = alpha[q]
			foil.Data[i].Cy[q] = cy[q]
			foil.Data[i].Cx[q] = cx[q]
		}

		alpha = alpha[:0]
		cy  = cy[:0]
		cx = cx[:0]
	}
	return foil
}

func FindCy(Alpha, Mach float64, foil Foil) float64{
	var(
		i, j int
		x1, x2 float64
		y1, y2 float64
		y float64
	)


	j = 0
	for i=1;i<len(foil.Mach);i++{
		j++
		if Mach <= foil.Mach[i]{
			break
		}

	}
	x1 = foil.Mach[j-1]
	x2 = foil.Mach[j]
	y1 = aprox(Alpha, foil.Data[j-1].Alpha, foil.Data[j-1].Cy)
	y2 = aprox(Alpha, foil.Data[j].Alpha, foil.Data[j].Cy)
	y = (Mach - x1) * (y2 - y1) / (x2 - x1) + y1
	return y
}

func FindCx(Alpha, Mach float64, foil Foil) float64{
	var(
		i, j int
		x1, x2 float64
		y1, y2 float64
		y float64
	)
	j = 0
	for i=1;i<len(foil.Mach);i++{
		j++
		if Mach <= foil.Mach[i]{
			break
		}
	}
	x1 = foil.Mach[j-1]
	x2 = foil.Mach[j]
	y1 = aprox(Alpha, foil.Data[j-1].Alpha, foil.Data[j-1].Cx)
	y2 = aprox(Alpha, foil.Data[j].Alpha, foil.Data[j].Cx)
	y = (Mach - x1) * (y2 - y1) / (x2 - x1) + y1
	return y
}

