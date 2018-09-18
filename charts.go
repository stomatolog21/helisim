package main

import (
	"bytes"
	"io/ioutil"
	"log"

	"github.com/wcharczuk/go-chart"
	"github.com/wcharczuk/go-chart/drawing"
)

func plot(X, Y []float64, filename, xaxis, yaxis string) {
	filename = ".\\pic\\" + filename
	graph := chart.Chart{
		XAxis: chart.XAxis{
			Name:      xaxis,
			NameStyle: chart.StyleShow(),
			Style:     chart.StyleShow(),
		},
		YAxis: chart.YAxis{
			Name:      yaxis,
			AxisType:  1,
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
	if err != nil {
		log.Println(err)
	}
	readbuf, err1 := ioutil.ReadAll(buffer)
	if err1 != nil {
		log.Println(err1)
	}
	err = ioutil.WriteFile(filename, readbuf, 0644)

}
