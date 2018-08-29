package main

import (
	"bytes"
	"fmt"
	"github.com/wcharczuk/go-chart"
	"github.com/wcharczuk/go-chart/drawing"
	"io/ioutil"
)

func plot(X, Y []float64, filename string){
	graph := chart.Chart{
		XAxis: chart.XAxis{
			Name:      "Радиус, м",
			NameStyle: chart.StyleShow(),
			Style:     chart.StyleShow(),
		},
		YAxis: chart.YAxis{
			Name:      "Индуктивная скорость, м/с",
			NameStyle: chart.StyleShow(),
			Style:     chart.StyleShow(),
		},

		Series: []chart.Series{
			chart.ContinuousSeries{
				Style: chart.Style{
					Show:        true,                           //note; if we set ANY other properties, we must set this to true.
					StrokeColor: drawing.ColorRed,               // will supercede defaults
					FillColor:   drawing.ColorRed.WithAlpha(64), // will supercede defaults
				},
				XValues: X,
				YValues: Y,
			},
		},
	}
	buffer := bytes.NewBuffer([]byte{})

	err := graph.Render(chart.PNG, buffer)
	readbuf, err1 := ioutil.ReadAll(buffer)
	err = ioutil.WriteFile(filename, readbuf, 0644)
	fmt.Println(err)
	fmt.Println(err1)
	fmt.Println(buffer)
}
