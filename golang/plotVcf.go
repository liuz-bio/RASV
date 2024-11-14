package main

import (
        "C"
        "fmt"
        "gonum.org/v1/plot"
        "gonum.org/v1/plot/plotter"
        "gonum.org/v1/plot/plotutil"
        "gonum.org/v1/plot/vg"
        gn_draw "gonum.org/v1/plot/vg/draw"
        "gonum.org/v1/plot/vg/vgimg"
        "gonum.org/v1/plot/vg/vgpdf"
        "gonum.org/v1/plot/vg/vgsvg"
        ifont "golang.org/x/image/font"
        "gonum.org/v1/plot/font"
        "gonum.org/v1/plot/text"
        "image"
        "image/color"
        "log"
        "math"
        "os"
        "time"
        "sort"
        "strconv"
        "strings"
)

//////////////////////////////////////////////////////////

const (
        // dlamchE is the machine epsilon. For IEEE this is 2^{-53}.
        dlamchE = 1.0 / (1 << 53)

        // dlamchB is the radix of the machine (the base of the number system).
        dlamchB = 2

        // dlamchP is base * eps.
        dlamchP = dlamchB * dlamchE
)

const (
        // free indicates no restriction on label containment.
        free = iota
        // containData specifies that all the data range lies
        // within the interval [label_min, label_max].
        containData
        // withinData specifies that all labels lie within the
        // interval [dMin, dMax].
        withinData
)

func talbotLinHanrahan(dMin, dMax float64, want int, containment int, Q []float64, w *weights, legibility func(lMin, lMax, lStep float64) float64) (values []float64, step, q float64, magnitude int) {
        const eps = dlamchP * 100

        if dMin > dMax {
                panic("labelling: invalid data range: min greater than max")
        }

        if Q == nil {
                Q = []float64{1, 5, 2, 2.5, 4, 3}
        }
        if w == nil {
                w = &weights{
                        simplicity: 0.25,
                        coverage:   0.2,
                        density:    0.5,
                        legibility: 0.05,
                }
        }
        if legibility == nil {
                legibility = unitLegibility
        }

        if r := dMax - dMin; r < eps {
                l := make([]float64, want)
                step := r / float64(want-1)
                for i := range l {
                        l[i] = dMin + float64(i)*step
                }
                magnitude = minAbsMag(dMin, dMax)
                return l, step, 0, magnitude
        }

        type selection struct {
                // n is the number of labels selected.
                n int
                // lMin and lMax are the selected min
                // and max label values. lq is the q
                // chosen.
                lMin, lMax, lStep, lq float64
                // score is the score for the selection.
                score float64
                // magnitude is the magnitude of the
                // label step distance.
                magnitude int
        }
        best := selection{score: -2}

outer:
        for skip := 1; ; skip++ {
                for _, q := range Q {
                        sm := maxSimplicity(q, Q, skip)
                        if w.score(sm, 1, 1, 1) < best.score {
                                break outer
                        }

                        for have := 2; ; have++ {
                                dm := maxDensity(have, want)
                                if w.score(sm, 1, dm, 1) < best.score {
                                        break
                                }

                                delta := (dMax - dMin) / float64(have+1) / float64(skip) / q

                                const maxExp = 309
                                for mag := int(math.Ceil(math.Log10(delta))); mag < maxExp; mag++ {
                                        step := float64(skip) * q * math.Pow10(mag)

                                        cm := maxCoverage(dMin, dMax, step*float64(have-1))
                                        if w.score(sm, cm, dm, 1) < best.score {
                                                break
                                        }

                                        fracStep := step / float64(skip)
                                        kStep := step * float64(have-1)

                                        minStart := (math.Floor(dMax/step) - float64(have-1)) * float64(skip)
                                        maxStart := math.Ceil(dMax/step) * float64(skip)
                                        for start := minStart; start <= maxStart && start != start-1; start++ {
                                                lMin := start * fracStep
                                                lMax := lMin + kStep

                                                switch containment {
                                                case containData:
                                                        if dMin < lMin || lMax < dMax {
                                                                continue
                                                        }
                                                case withinData:
                                                        if lMin < dMin || dMax < lMax {
                                                                continue
                                                        }
                                                case free:
                                                        // Free choice.
                                                }

                                                score := w.score(
                                                        simplicity(q, Q, skip, lMin, lMax, step),
                                                        coverage(dMin, dMax, lMin, lMax),
                                                        density(have, want, dMin, dMax, lMin, lMax),
                                                        legibility(lMin, lMax, step),
                                                )
                                                if score > best.score {
                                                        best = selection{
                                                                n:         have,
                                                                lMin:      lMin,
                                                                lMax:      lMax,
                                                                lStep:     float64(skip) * q,
                                                                lq:        q,
                                                                score:     score,
                                                                magnitude: mag,
                                                        }
                                                }
                                        }
                                }
                        }
                }
        }

        if best.score == -2 {
                l := make([]float64, want)
                step := (dMax - dMin) / float64(want-1)
                for i := range l {
                        l[i] = dMin + float64(i)*step
                }
                magnitude = minAbsMag(dMin, dMax)
                return l, step, 0, magnitude
        }

        l := make([]float64, best.n)
        step = best.lStep * math.Pow10(best.magnitude)
        for i := range l {
                l[i] = best.lMin + float64(i)*step
        }
        return l, best.lStep, best.lq, best.magnitude
}

// minAbsMag returns the minumum magnitude of the absolute values of a and b.
func minAbsMag(a, b float64) int {
        return int(math.Min(math.Floor(math.Log10(math.Abs(a))), (math.Floor(math.Log10(math.Abs(b))))))
}

// simplicity returns the simplicity score for how will the curent q, lMin, lMax,
// lStep and skip match the given nice numbers, Q.
func simplicity(q float64, Q []float64, skip int, lMin, lMax, lStep float64) float64 {
        const eps = dlamchP * 100

        for i, v := range Q {
                if v == q {
                        m := math.Mod(lMin, lStep)
                        v = 0
                        if (m < eps || lStep-m < eps) && lMin <= 0 && 0 <= lMax {
                                v = 1
                        }
                        return 1 - float64(i)/(float64(len(Q))-1) - float64(skip) + v
                }
        }
        panic("labelling: invalid q for Q")
}

// maxSimplicity returns the maximum simplicity for q, Q and skip.
func maxSimplicity(q float64, Q []float64, skip int) float64 {
        for i, v := range Q {
                if v == q {
                        return 1 - float64(i)/(float64(len(Q))-1) - float64(skip) + 1
                }
        }
        panic("labelling: invalid q for Q")
}

// coverage returns the coverage score for based on the average
// squared distance between the extreme labels, lMin and lMax, and
// the extreme data points, dMin and dMax.
func coverage(dMin, dMax, lMin, lMax float64) float64 {
        r := 0.1 * (dMax - dMin)
        max := dMax - lMax
        min := dMin - lMin
        return 1 - 0.5*(max*max+min*min)/(r*r)
}

// maxCoverage returns the maximum coverage achievable for the data
// range.
func maxCoverage(dMin, dMax, span float64) float64 {
        r := dMax - dMin
        if span <= r {
                return 1
        }
        h := 0.5 * (span - r)
        r *= 0.1
        return 1 - (h*h)/(r*r)
}

// density returns the density score which measures the goodness of
// the labelling density compared to the user defined target
// based on the want parameter given to talbotLinHanrahan.
func density(have, want int, dMin, dMax, lMin, lMax float64) float64 {
        rho := float64(have-1) / (lMax - lMin)
        rhot := float64(want-1) / (math.Max(lMax, dMax) - math.Min(dMin, lMin))
        if d := rho / rhot; d >= 1 {
                return 2 - d
        }
        return 2 - rhot/rho
}

// maxDensity returns the maximum density score achievable for have and want.
func maxDensity(have, want int) float64 {
        if have < want {
                return 1
        }
        return 2 - float64(have-1)/float64(want-1)
}

// unitLegibility returns a default legibility score ignoring label
// spacing.
func unitLegibility(_, _, _ float64) float64 {
        return 1
}

// weights is a helper type to calcuate the labelling scheme's total score.
type weights struct {
        simplicity, coverage, density, legibility float64
}

// score returns the score for a labelling scheme with simplicity, s,
// coverage, c, density, d and legibility l.
func (w *weights) score(s, c, d, l float64) float64 {
        return w.simplicity*s + w.coverage*c + w.density*d + w.legibility*l
}

// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
func getSortedValues(input map[int]*plot.Plot) []*plot.Plot {
        keys := make([]int, 0, len(input))
        for k := range input {
                keys = append(keys, k)
        }
        sort.Ints(keys)
        values := make([]*plot.Plot, len(keys))
        for i, k := range keys {
                values[i] = input[k]
        }
        return values
}

type integerTicks struct{}

func (integerTicks) Ticks(min, max float64) []plot.Tick {
        if max <= min {
                panic("illegal range")
        }

        const suggestedTicks = 3

        labels, step, q, mag := talbotLinHanrahan(min, max, suggestedTicks, withinData, nil, nil, nil)
        majorDelta := step * math.Pow10(mag)
        if q == 0 {
                // Simple fall back was chosen, so
                // majorDelta is the label distance.
                majorDelta = labels[1] - labels[0]
        }

        ticks := make([]plot.Tick, len(labels))
        for i, v := range labels {
                //ticks[i] = plot.Tick{Value: v, Label: strconv.FormatFloat(v, fc, prec, 64)}
                ticks[i] = plot.Tick{Value: float64(int(v)), Label: strconv.Itoa(int(v))}
        }

        var minorDelta float64
        // See talbotLinHanrahan for the values used here.
        switch step {
        case 1, 2.5:
                minorDelta = majorDelta / 5
        case 2, 3, 4, 5:
                minorDelta = majorDelta / step
        default:
                if majorDelta/2 < dlamchP {
                        return ticks
                }
                minorDelta = majorDelta / 2
        }

        // Find the first minor tick not greater
        // than the lowest data value.
        var i float64
        for labels[0]+(i-1)*minorDelta > min {
                i--
        }
        // Add ticks at minorDelta intervals when
        // they are not within minorDelta/2 of a
        // labelled tick.
        for {
                val := labels[0] + i*minorDelta
                if val > max {
                        break
                }
                found := false
                for _, t := range ticks {
                        if math.Abs(t.Value-val) < minorDelta/2 {
                                found = true
                        }
                }
                if !found {
                        ticks = append(ticks, plot.Tick{Value: val})
                }
                i++
        }
        return ticks
}
//////////////////////////////////////////////////////////////////////////////////////////////////////
func creatBar(hp1 int, hp2 int)(*plot.Plot){
    hpa1 := plotter.Values{float64(hp1)}
    hpa2 := plotter.Values{float64(hp2)}

    p := plot.New()

    p.Y.Label.Text = "Coverage" //"Support reads"

    w := vg.Points(20)
    hp1Bar, err := plotter.NewBarChart(hpa1, w)
    if err != nil {
        log.Fatal(err)
    }
    hp1Bar.LineStyle.Width = vg.Length(0)
    hp1Bar.Color = plotutil.Color(0)
    hp1Bar.Offset = -w

    hp2Bar, err := plotter.NewBarChart(hpa2, w)
    if err != nil {
        log.Fatal(err)
    }
    hp2Bar.LineStyle.Width = vg.Length(0)
    hp2Bar.Color = plotutil.Color(2)
    hp2Bar.Offset = w

    p.Add(hp1Bar,hp2Bar)
    p.Legend.Add("hp1", hp1Bar)
    p.Legend.Add("hp2", hp2Bar)
    //p.Legend.Top = true
    p.Legend.Top = true
    p.Legend.Left = true
    //p.NominalX("GT")
    return p
}


func randomDataFloat64(n [][]float64) plotter.XYs {
        pts := make(plotter.XYs, len(n))
        for i := 0; i < len(n); i++ {
                x := n[i][0]
                y := n[i][1]
                pts[i].X = x
                pts[i].Y = y
        }
        return pts
}

func arrowDat(rawPos [][]int, strand string)([][][]float64){
    part_arrow := make([][][]float64,0)
    start_pos := float64(rawPos[0][0])
    geneSize := float64(rawPos[1][0] - rawPos[0][0])
    arrow_win := geneSize/6
    //angle := 10.0 // 角度（以度为单位）
    //c := float64(5) // \
    //radian := angle * (math.Pi / 180.0)
    //sinValue := math.Sin(radian)
    //cosValue := math.Cos(radian)
    a := 0.2 //h |
    //b := a*cosValue/sinValue //w _
    b := geneSize*0.002

    for i := 1; i <= 5; i++ {
        start :=  start_pos + arrow_win*float64(i)
        if strand == "+" { //rigth
            part_arrow = append(part_arrow, [][]float64{{start+b,float64(rawPos[0][1])}, {start,float64(rawPos[0][1])+a}})
            part_arrow = append(part_arrow, [][]float64{{start+b,float64(rawPos[0][1])}, {start,float64(rawPos[0][1])-a}})
        }else{
            part_arrow = append(part_arrow, [][]float64{{start,float64(rawPos[0][1])}, {start+b,float64(rawPos[0][1])+a}})
            part_arrow = append(part_arrow, [][]float64{{start,float64(rawPos[0][1])}, {start+b,float64(rawPos[0][1])-a}})
        }
    }
    return part_arrow
}

func plotArrow(rawPos [][]int,
               p *plot.Plot,
               arrow_color color.RGBA,
               element_size float64,
               strand string)(*plot.Plot){

    part_arrow := arrowDat(rawPos, strand)
    for _, pos := range part_arrow{
        line_xy := randomDataFloat64(pos)
        arrow, _ := plotter.NewLine(line_xy)
        arrow.Color = arrow_color
        arrow.Width = vg.Points(element_size)
        p.Add(arrow)
    }
    return p
}

func plotLine(rawPos [][]int,
              p *plot.Plot,
              legend plot.Legend,
              element_color color.RGBA,
              element_size float64,
              strand string,
              legendAdded  map[string]bool,
              lineType string)(*plot.Plot, plot.Legend, map[string]bool){

    line_xy := randomData(rawPos)
    element, _ := plotter.NewLine(line_xy)
    element.Color = element_color
    element.Width = vg.Points(element_size)
    p.Add(element)
    if !legendAdded[lineType]{
        legend.Add(fmt.Sprintf("%s",  lineType), element)
        legendAdded[lineType] = true
        //fmt.Println("plotLine:", lineType)
    }

    if lineType == "gene"{
        p = plotArrow(rawPos, p, element_color, element_size, strand)
    }
    return p, legend, legendAdded
}

func plotSnp(rawPos [][]int,
             p *plot.Plot,
             legend plot.Legend,
             scatter_color color.RGBA,
             scatter_size float64,
             legendAdded  map[string]bool,
             scatterType string)(*plot.Plot, plot.Legend, map[string]bool){

    scatter_xy := randomData(rawPos)
    scatter, err := plotter.NewScatter(scatter_xy)
    if err != nil {
        log.Fatal(err)
    }
    // 设置竖线点的形状为竖线
    scatter.GlyphStyle.Shape = gn_draw.BoxGlyph{} // PlusGlyph{}
    scatter.GlyphStyle.Color = scatter_color
    //scatter.GlyphStyle.Radius  = scatter_size
    p.Add(scatter)
    if !legendAdded[scatterType]{
        legend.Add(fmt.Sprintf("%s",  scatterType), scatter)
        legendAdded[scatterType] = true
    }
    return p, legend, legendAdded
}

func snpDat(snp_str string, h int)([][]int){
    snps := strings.Split(snp_str, ",")
    out_snps := make([][]int,0)
    for _, pos_x := range snps{
        posX, err := strconv.Atoi(pos_x)
        if err != nil {
            return out_snps
            //fmt.Println("无法将字符串转换为整数:snpDat", err)
        }
        out_snps = append(out_snps, []int{posX, h})
    }
    return out_snps
}

type Gene struct {
    geneID string
    gene_strand string
    genePos []int
    tranIDs []string
    tran_exonPos map[string][][]int
    tran_tranPos map[string][]int
    tran_fivePos map[string][][]int
    tran_threePos map[string][][]int
}

func strTolist(str string)([]int){
    tmp := strings.Split(str, ",")
    out := make([]int, 0)
    for _, tmp_pos := range tmp{
        pos, err := strconv.Atoi(tmp_pos)
        if err != nil {
            return out
            //fmt.Println("无法将字符串转换为整数:strTolist", err)
        }
        out = append(out, pos)
    }
    return out
}

func dealGenesInfo(geneIDstr string,
                   tranIDstr string,
                   strandstr string,
                   genePos_str string,
                   tranPos_str string,
                   exonPos_str string,
                   fivePos_str string,
                   threePos_str string)([]Gene){

    allGenesInfo := make([]Gene,0)
    if geneIDstr == "" {
        return allGenesInfo
    }
    geneIDsList := strings.Split(geneIDstr,";")
    strandsList  := strings.Split(strandstr,";")
    genePosList := strings.Split(genePos_str,";")
    tranIDsList := strings.Split(tranIDstr,";")
    tranPosList := strings.Split(tranPos_str,";")
    exonPosList := strings.Split(exonPos_str,";")
    fivePosList := strings.Split(fivePos_str,";")
    threePosList := strings.Split(threePos_str,";")

    for i, geneID := range geneIDsList{
        fmt.Println(geneID,strandsList, i)
        strand := strandsList[i]
        fmt.Println("genePosList:",genePosList[i])
        genePos := strTolist(genePosList[i])
        tranIDs := strings.Split(tranIDsList[i], ",")
        tmp_tran_Pos := strings.Split(tranPosList[i], "|")
        tmp_exon_Pos := strings.Split(exonPosList[i], "|")
        tmp_five_Pos := strings.Split(fivePosList[i], "|")
        tmp_three_Pos := strings.Split(threePosList[i], "|")

        tran_tranPos := make(map[string][]int)
        tran_fivePos := make(map[string][][]int)
        tran_threePos := make(map[string][][]int)
        tran_exonPos := make(map[string][][]int)
        for ii, tranid := range tranIDs{
            tranPos := strTolist(tmp_tran_Pos[ii])
            tran_tranPos[tranid] = tranPos

            five_Pos := strings.Split(tmp_five_Pos[ii], ":")
            fivePos := make([][]int,0)
            for _, tmp_pos := range five_Pos{
                fivePos = append(fivePos, strTolist(tmp_pos))
            }
            tran_fivePos[tranid] = fivePos

            three_Pos := strings.Split(tmp_three_Pos[ii], ":")
            threePos := make([][]int,0)
            for _, tmp_pos := range three_Pos{
                threePos = append(threePos, strTolist(tmp_pos))
            }
            tran_threePos[tranid] = threePos

            exon_Pos := strings.Split(tmp_exon_Pos[ii], ":")
            exonPos := make([][]int,0)
            for _, tmp_pos := range exon_Pos{
                exonPos = append(exonPos, strTolist(tmp_pos))
            }
            tran_exonPos[tranid] = exonPos
        }

        newGene := Gene {geneID: geneID,
                         gene_strand: strand,
                         genePos: genePos,
                         tranIDs: tranIDs,
                         tran_exonPos: tran_exonPos,
                         tran_tranPos: tran_tranPos,
                         tran_fivePos: tran_fivePos,
                         tran_threePos: tran_threePos,
                         }
        allGenesInfo = append(allGenesInfo, newGene)
    }

    return allGenesInfo
}


func plotGeneSnp(allGenesInfo []Gene, snp_str string)(*plot.Plot, int, bool){
    colors := make(map[string]color.RGBA)
    colors["gene"] = color.RGBA{R: 166, G: 166, B: 166, A: 255}      //166,166,166
    colors["transcript"] = color.RGBA{R: 252, G: 239, B: 156, A: 255} //252,239,156
    colors["exon"] = color.RGBA{R: 255, G: 190, B: 78, A: 255}      //255,190,78
    colors["five_utr"] = color.RGBA{R: 148, G: 230, B: 236, A: 255} //148,230,236
    colors["three_utr"] = color.RGBA{R: 83, G: 170, B: 240, A: 255} //83,170,240
    colors["snp"] = color.RGBA{R: 0, G: 0, B: 0, A: 50}

    elements_size := make(map[string]float64)
    elements_size["gene"] = float64(1)
    elements_size["transcript"] = float64(2)
    elements_size["exon"] = float64(6)
    elements_size["five_utr"] = float64(6)
    elements_size["three_utr"] = float64(6)
    elements_size["snp"] = float64(0.1)

    flag := false
    h := 0
    idList := make([]string, 0)

    p := plot.New()
    legendAdded := make(map[string]bool)
    legend := plot.NewLegend()
    legend.ThumbnailWidth = 0.5 * vg.Inch

    snps := snpDat(snp_str, h)
    if len(snps)>0{
        idList = append(idList, "clinvar_snp")
        p, legend, legendAdded  = plotSnp(snps, p, legend,
                                          colors["snp"],
                                          elements_size["snp"],
                                          legendAdded, "clinvar_snp")
        h = h +1
        flag = true
    }
    if len(allGenesInfo) >0 {
        for _, gene := range allGenesInfo{
            idList = append(idList, gene.geneID)
            strand := gene.gene_strand
            gene_pos := gene.genePos
            fmt.Println("gene_pos:",gene_pos)
            p, legend, legendAdded = plotLine([][]int{{gene_pos[0],h},{gene_pos[1],h}}, p, legend,
                                              colors["gene"],
                                              elements_size["gene"],
                                              strand, legendAdded,  "gene")
            flag = true
            h = h +1
            for _, tranID := range gene.tranIDs{
                if tranID == "" {
                     continue
                }
                idList = append(idList, tranID)
                tran_pos := gene.tran_tranPos[tranID]
                five_pos := gene.tran_fivePos[tranID]
                three_pos := gene.tran_threePos[tranID]
                p, legend, legendAdded = plotLine([][]int{{tran_pos[0],h},{tran_pos[1],h}}, p, legend,
                                         colors["transcript"],
                                         elements_size["transcript"],
                                         strand, legendAdded, "Intron") //"transcript"

                for _, exon_pos := range gene.tran_exonPos[tranID]{
                    p, legend, legendAdded = plotLine([][]int{{exon_pos[0],h},{exon_pos[1],h}}, p, legend,
                                             colors["exon"],
                                             elements_size["exon"],
                                             strand, legendAdded, "exon")
                }

                for _, five_p := range five_pos{
                    if five_p[0] > 0{
                        p, legend, legendAdded = plotLine([][]int{{five_p[0],h},{five_p[1],h}}, p, legend,
                                                 colors["five_utr"],
                                                 elements_size["five_utr"],
                                                 strand, legendAdded, "5'UTR")
                    }
                }

                for _, three_p := range three_pos{
                    if three_p[0] > 0{
                        p, legend, legendAdded = plotLine([][]int{{three_p[0],h},{three_p[1],h}}, p, legend,
                                                 colors["three_utr"],
                                                 elements_size["three_utr"],
                                                 strand, legendAdded, "3'UTR")
                    }
                }
                h = h +1
            }
        }
    }

    if flag{
        p.NominalY(idList...)
        legend.Top = true
        legend.Left = true
        p.Legend = legend

    }

    return p, h, flag
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
func createPlots(allSignal map[int]map[int][][]int,
                 chrom string,
                 svtype string,
                 start int,
                 stop int,
                 flanks []int,
                 flanksInvalid map[int]int,
                 flankSize int,
                 legendFlag bool) ([]*plot.Plot, *plot.Plot) {

        figDat := make(map[int]map[int][][][]int)
        flanksInvalidSize := len(flanksInvalid)
        hp1, hp2 := 0, 0
        //fmt.Println("flanksInvalid:", flanksInvalid)
        if flanksInvalidSize == 0 {
                figDat, hp1, hp2 = figData(allSignal, flanks, start, stop)
        } else if flanksInvalidSize < flankSize {
                figDat, hp1, hp2 = figData(allSignal, flanks, start, stop)
                for _, flank := range flanksInvalid {
                        flankData := make(map[int][][][]int)
                        flankData[6] = [][][]int{{{start, 1}, {start, 25},{-1}}, {{stop, 1}, {stop, 25},{-1}}} //[][][]int{{{start, 1}, {start, nu}, {-1}}, {{stop, 1}, {stop, nu}, {-1}}}
                        figDat[flank] = flankData
                }
        } else {
                for _, flank := range flanksInvalid {
                        flankData := make(map[int][][][]int)
                        flankData[6] = [][][]int{{{start, 1}, {start, 25},{-1}}, {{stop, 1}, {stop, 25},{-1}}} //[][][]int{{{x1, y}, {x2, y}, {mpQ}}}
                        figDat[flank] = flankData
                }
        }

        pBar := creatBar(hp1, hp2)

       // colors := []color.RGBA{
       //         {R: 1, G: 193, B: 110, A: 255},   //1 L 
       //         {R: 3, G: 126, B: 243, A: 255},   //2 Del
       //         {R: 248, G: 90, B: 64, A: 255},   //3 Ins
       //         {R: 255, G: 200, B: 69, A: 255},  //4 R 
       //         {R: 255, G: 255, B: 255, A: 255},  //alignment 5 gray
       //         {R: 100, G: 100, B: 100, A: 255}, //breakpoint 6 
       // }

        colors := map[string][]color.RGBA{
                  "DEL":{ 
                          {R: 1, G: 193, B: 110, A: 255},   //1 L
                          {R: 3, G: 126, B: 243, A: 255},   //2 Del
                          {R: 248, G: 90, B: 64, A: 255},   //3 Ins
                          {R: 255, G: 200, B: 69, A: 255},  //4 R
                          {R: 255, G: 255, B: 255, A: 255},  //alignment 5 gray
                          {R: 100, G: 100, B: 100, A: 255}, //breakpoint 6
                        },
                  "INS":{
                          {R: 1, G: 193, B: 110, A: 255},   //1 L
                          {R: 3, G: 126, B: 243, A: 255},   //2 Del
                          {R: 248, G: 90, B: 64, A: 255},   //3 Ins
                          {R: 255, G: 200, B: 69, A: 255},  //4 R
                          {R: 255, G: 255, B: 255, A: 255},  //alignment 5 gray
                          {R: 100, G: 100, B: 100, A: 255}, //breakpoint 6
                        },
                  "DUP":{
                          {R: 170, G:221, B:34,  A: 255},  //1 L
                          {R: 3, G: 126, B: 243, A: 255},   //2 Del
                          {R: 248, G: 90, B: 64, A: 255},   //3 Ins
                          {R: 145, G:70,  B:255,   A: 255},  //4 R
                          {R: 255, G: 255, B: 255, A: 255}, //alignment 5 gray
                          {R: 100, G: 100, B: 100, A: 255}, //breakpoint 6
                        },
                  "INV":{
                          {R: 255, G:120,  B:203, A: 255},  //1 L
                          {R: 3, G: 126, B: 243, A: 255},   //2 Del
                          {R: 248, G: 90, B: 64, A: 255},   //3 Ins
                          {R: 0,   G:194,  B:224, A: 255},  //4 R
                          {R: 255, G: 255, B: 255, A: 255}, //alignment 5 gray
                          {R: 100, G: 100, B: 100, A: 255}, //breakpoint 6
                        },
                   "BND":{
                          {R: 153, G:0,    B:0,   A: 255},  //1 L
                          {R: 3, G: 126, B: 243, A: 255},   //2 Del
                          {R: 248, G: 90, B: 64, A: 255},   //3 Ins
                          {R: 173, G:132,  B:31,  A: 255},  //4 R
                          {R: 255, G: 255, B: 255, A: 255}, //alignment 5 gray
                          {R: 100, G: 100, B: 100, A: 255}, //breakpoint 6
                        },
                  }

        // 自定义渐变的颜色区间
        //f3f6f6                        243,246,246
        //eaeded                        234,237,237
        //d5dbdb                        213,219,219
        //bfc9ca                        191,201,202
        //a9b7b8                        169,183,184
        //95a5a6                        149,165,166
        alignColors := []color.RGBA{
                       {R: 243, G: 246, B: 246, A: 255}, //mapQ 0-10
                       {R: 234, G: 237, B: 237, A: 255}, //mapQ 10-20
                       {R: 213, G: 219, B: 219, A: 255}, //mapQ 20-30
                       {R: 191, G: 201, B: 202, A: 255}, //mapQ 30-40
                       {R: 169, G: 183, B: 184, A: 255}, //mapQ 40-50
                       {R: 149, G: 165, B: 166, A: 255}, //mapQ 50-60
        }

        //               {R: 97,  G:111,  B:103, A: 255} //1 L
        //               {R: 198, G:129,  B:67,  A: 255} //4 R

       //                {R: 255, G:120,  B:203, A: 255} //1 L
       //                {R: 0,   G:194,  B:224, A: 255} //4 R

       //                {R: 153, G:0,    B:0,   A: 255} //1 L
       //                {R: 173, G:132,  B:31,  A: 255} //4 R

       //                {R: 169, G:60,   B:59,  A: 255} //1 L
       //                {R: 168, G:150,  B:160, A: 255} //4 R

        legendTitle := map[int]string{
                1: "L skip segment",
                2: "DEL signal",
                3: "INS signal",
                4: "R skip segment",
                5: "Align read",
                6: "BND",
        }

        imgFiles := make(map[int]*plot.Plot, len(figDat))
        colorLinesMap := make(map[int]map[int][][][]int, len(figDat))
        for flank, colorLines := range figDat {
                colorLinesMap[flank] = make(map[int][][][]int)
                for colorId, lines := range colorLines {
                        colorLinesMap[flank][colorId] = lines
                }
        }
        //var figName string
        legendAdded := make(map[color.RGBA]bool) // 用于跟踪已添加的颜色

        colorIds := []int{5, 6, 1, 4, 3, 2}
        switch svtype{
        case "DEL":
            colorIds = []int{5, 6, 2, 1, 4}
        case "INS":
            colorIds = []int{5, 6, 3, 1, 4}    
        case "DUP":
            colorIds = []int{5, 6, 1, 4, 2, 3}
        case "INV":
            colorIds = []int{5, 6, 1, 4, 2, 3}
        case "BND":
            colorIds = []int{5, 6, 1, 4, 2, 3}
        case "TRA":
            colorIds = []int{5, 6, 1, 4, 2, 3}
        }

        for i, flank := range flanks {
                colorLines := colorLinesMap[flank]
                p := plot.New()
                legend := plot.NewLegend()
                legend.ThumbnailWidth = 0.5 * vg.Inch
                //colorIds := []int{5, 6, 1, 4, 3, 2}
                //colorLines
                for _, colorId := range colorIds {
                        lines := colorLines[colorId]
                        for _, line := range lines {
                                l_size := len(line)-1
                                lin := line[:l_size]
                                //fmt.Println("line:", line)
                                mapQ:= line[l_size][0]
                                xy := randomData(lin)
                                line, err := plotter.NewLine(xy)

                                if err != nil {
                                        log.Fatal(err)
                                }

                                if mapQ < 0{
                                    line.Color = colors[svtype][colorId-1]
                                }else{
                                    // 计算线条颜色
                                    //int(math.Abs(float64(mapQ)))
                                    //fmt.Println("alignColors:", int(int(math.Abs(float64(mapQ-1))) / 10), mapQ)
                                    mapQ_c := alignColors[int(int(math.Abs(float64(mapQ-1))) / 10)]
                                    line.Color = mapQ_c
                                    //line.LineStyle.Color = mapQ_c
                                }

                                line.Width = vg.Points(2)
                                p.Add(line)
                                // 只为最后一张图片添加图例
                                if i == flankSize-1 && !legendAdded[colors[svtype][colorId-1]] && legendFlag {
                                        if colorId == 6 {
                                                legend.Add(fmt.Sprintf("%s %s", svtype, legendTitle[colorId]), line)
                                                legendAdded[colors[svtype][colorId-1]] = true
                                        } else {
                                                legend.Add(fmt.Sprintf("%s", legendTitle[colorId]), line)
                                                legendAdded[colors[svtype][colorId-1]] = true
                                        }
                                }
                        }
                }
                p.X.Tick.Marker = integerTicks{}
                p.Y.Tick.Marker = integerTicks{}
                legend.Top = true
                legend.Left = true
                p.Legend = legend
                // 添加X和Y轴标题
                imgFiles[flank] = p
        }
        return getSortedValues(imgFiles), pBar
}

func flipImageXAxis(img image.Image) image.Image {
        width := img.Bounds().Dx()
        height := img.Bounds().Dy()

        flippedImg := image.NewRGBA(image.Rect(0, 0, width, height))

        for y := 0; y < height; y++ {
                for x := 0; x < width; x++ {
                        flippedImg.Set(x, y, img.At(x, height-y-1))
                }
        }

        return flippedImg
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
func CombineSV(plots map[string][]*plot.Plot,
               dcc gn_draw.Canvas,
               canvas *vgimg.Canvas,
               imgwith vg.Length,
               imgheight vg.Length,
               geneHeight vg.Length,
               Xdistance vg.Length,
               Ydistance vg.Length,
               textStyleFlank text.Style,
               textStyleSample text.Style,
               xlabel string,
               infos string,
               flankSize int,
               dpi int,
               img_Y int,
               flanksID []string,
               samples []string)(gn_draw.Canvas, *vgimg.Canvas){

        var tmp_legend plot.Legend
        ii := img_Y

        geneHeight_tmp := geneHeight*vg.Length(img_Y)

        for _, sample := range samples {
                ps := plots[sample]
                for i, p := range ps {
                        if ii == img_Y {
                                p.X.Label.Text = xlabel
                        }
                        if i == 0 {
                                p.Y.Label.Text = "Coverage" //"Reads of Number"
                                if ii == img_Y+len(plots)-1 {
                                       var x2 vg.Length = vg.Length(i)*(imgwith+Xdistance) + 0.5*vg.Inch
                                       var y2 vg.Length = vg.Length(ii+1-img_Y)*(imgheight+Ydistance) + 0.5*vg.Inch + geneHeight_tmp
                                       dcc.FillText(textStyleFlank, vg.Point{X: x2, Y: y2 + 0.2*vg.Inch}, infos)
                                }
                        }


                        if ii == img_Y+len(plots)-1 {
                                var x0 vg.Length = vg.Length(i)*(imgwith+Xdistance) + 0.5*vg.Inch
                                var y0 vg.Length = vg.Length(ii+1-img_Y)*(imgheight+Ydistance) + 0.5*vg.Inch + geneHeight_tmp
                                dcc.FillText(textStyleFlank, vg.Point{X: x0 + imgwith/2, Y: y0 + 0.2*vg.Inch}, flanksID[i])
                        }

                        if i == flankSize-1 {
                                tmp_legend = p.Legend
                                p.Legend = plot.NewLegend()
                        } else {
                                p.Legend = plot.NewLegend()
                        }

                        c := vgimg.NewWith(
                                vgimg.UseWH(imgwith, imgheight),
                                vgimg.UseDPI(dpi),
                        )

                        dc := gn_draw.New(c)
                        p.Draw(dc)
                        if i == 0 {
                                var x1 vg.Length = vg.Length(i)*(imgwith+Xdistance) + 0.5*vg.Inch
                                var y1 vg.Length = vg.Length(ii+1-img_Y)*(imgheight+Ydistance) + 0.5*vg.Inch + geneHeight_tmp
                                dcc.FillText(textStyleSample, vg.Point{X: x1 + 0.8*vg.Inch, Y: y1 - 0.22*vg.Inch}, sample)
                        }

                        if ii == img_Y+len(plots)-1 && i == flankSize-1 {
                                p.Legend = tmp_legend
                                p.Legend.Draw(gn_draw.Crop(dcc, dcc.Size().X-2*vg.Inch, 0, 0, -0.5*vg.Inch))
                        }

                        img := vgimg.PngCanvas{Canvas: c}
                        var x vg.Length = vg.Length(i)*(imgwith+Xdistance) + 0.5*vg.Inch
                        var y vg.Length = vg.Length(ii-img_Y)*(imgheight+Ydistance) + 0.5*vg.Inch + geneHeight_tmp
                        rect := vg.Rectangle{
                                Min: vg.Point{X: x, Y: y},
                                Max: vg.Point{X: x + imgwith, Y: y + imgheight},
                        }
                        canvas.DrawImage(rect, img.Image())
                }
                ii += 1
        }

        return dcc, canvas
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
func CombineGT(plotsGT map[string]*plot.Plot,
               dcc gn_draw.Canvas,
               canvas *vgimg.Canvas,
               hpWith vg.Length,
               imgwith vg.Length,
               imgheight vg.Length,
               geneHeight vg.Length,
               Xdistance vg.Length,
               Ydistance vg.Length,
               textStyleFlank text.Style,
               flankSize int,
               dpi int,
               img_Y int,
               samples []string,
               samplesGT map[string]string)(gn_draw.Canvas, *vgimg.Canvas){

        var tmp_legend plot.Legend
        ii := flankSize
        geneHeight_tmp := geneHeight*vg.Length(img_Y)
        for i, sample := range samples {
                //fmt.Println("sample:", sample, "Genotype")
                p := plotsGT[sample]
                p.HideX()
                if i == 0 {
                    p.X.Label.Text = "GT"
                }

                if i == len(samples)-1 {
                    tmp_legend = p.Legend
                    p.Legend = plot.NewLegend()
                    var x0 vg.Length = vg.Length(ii)*(imgwith+Xdistance) + 0.3*vg.Inch
                    var y0 vg.Length = vg.Length(i+1)*(imgheight+Ydistance) + 0.5*vg.Inch + geneHeight_tmp
                    dcc.FillText(textStyleFlank, vg.Point{X: x0 + hpWith/2, Y: y0 + 0.2*vg.Inch}, "Genotype")
                }else{
                    p.Legend = plot.NewLegend()
                }

                c := vgimg.NewWith(
                                   vgimg.UseWH(hpWith, imgheight),
                                   vgimg.UseDPI(dpi),
                                  )

                dc := gn_draw.New(c)
                p.Draw(dc)

                if i == len(samples)-1 {
                        p.Legend = tmp_legend
                        p.Legend.Draw(gn_draw.Crop(dcc, dcc.Size().X-2*vg.Inch, 0, 0, -1.5*vg.Inch))
                }

                img := vgimg.PngCanvas{Canvas: c}
                var x vg.Length = vg.Length(ii)*(imgwith+Xdistance) + 0.5*vg.Inch
                var y vg.Length = vg.Length(i)*(imgheight+Ydistance) + 0.5*vg.Inch + geneHeight_tmp+ vg.Length(0.2)*vg.Inch
                rect := vg.Rectangle{
                                     Min: vg.Point{X: x, Y: y},
                                     Max: vg.Point{X: x + hpWith, Y: y + imgheight},
                        }
                canvas.DrawImage(rect, img.Image())

                var x2 vg.Length = vg.Length(ii)*(imgwith+Xdistance) + 0.4*vg.Inch + hpWith
                var y2 vg.Length = vg.Length(i)*(imgheight+Ydistance) + 0.5*vg.Inch + geneHeight_tmp + imgheight/2 + Ydistance
                dcc.FillText(textStyleFlank, vg.Point{X: x2, Y: y2 }, "GT: "+samplesGT[sample])
        }

        return dcc, canvas
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
func CombineGene(plotsGene []*plot.Plot,
                 dcc gn_draw.Canvas,
                 canvas *vgimg.Canvas,
                 imgwith vg.Length,
                 imgheight vg.Length,
                 textStyleFlank text.Style,
                 dpi int,
                 samples []string,)(gn_draw.Canvas, *vgimg.Canvas){

        ii := 0
        for i, p := range plotsGene{
            c := vgimg.NewWith(
                               vgimg.UseWH(imgwith, imgheight),
                               vgimg.UseDPI(dpi),
                              )

            tmp_legend := p.Legend
            p.Legend = plot.NewLegend()

            dc := gn_draw.New(c)
            p.Draw(dc)

            if i == 0 {
                       p.Legend = tmp_legend
                       p.Legend.Draw(gn_draw.Crop(dcc, dcc.Size().X-2*vg.Inch, 0, 0, -2*vg.Inch))
            }

            img := vgimg.PngCanvas{Canvas: c}
            var x vg.Length = vg.Length(i) + 0.5*vg.Inch
            var y vg.Length = vg.Length(ii) + 0.5*vg.Inch
            rect := vg.Rectangle{
                                 Min: vg.Point{X: x, Y: y},
                                 Max: vg.Point{X: x + imgwith, Y: y + imgheight},
                    }
            canvas.DrawImage(rect, img.Image())
        }

        return dcc, canvas
}
// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
func CombineAndSave(plots map[string][]*plot.Plot,
                    plotsGT map[string]*plot.Plot,
                    plotsGene []*plot.Plot,
                    gene_h int,
                    gene_plot *plot.Plot,
                    gene_flag bool,
                    filename string,
                    with int,
                    height int,
                    flankSize int,
                    xpad float32,
                    ypad float32,
                    xlabel string,
                    flanksID []string,
                    samples []string,
                    samplesGT map[string]string,
                    infos string,
                    dpi int) {

        Ydistance := vg.Length(ypad) * vg.Inch
        Xdistance := vg.Length(xpad) * vg.Inch

        hpWith    := vg.Length(with/2) * vg.Inch
        geneHeight := vg.Length(float64(gene_h)*0.15 +1.0) * vg.Inch

        imgwith := vg.Length(with) * vg.Inch
        geneWith := vg.Length(flankSize)*(imgwith+Xdistance) + hpWith
        imgheight := vg.Length(height) * vg.Inch

        sampleSize := vg.Length(len(plots))
        canvasWith := vg.Length(flankSize)*(imgwith+Xdistance) + 2*vg.Inch + vg.Length(0.8)*(imgwith+Xdistance)
        canvasHeight := sampleSize*(imgheight+Ydistance) + 2*vg.Inch


        textStyleSample := text.Style{
                Color:   color.Black,
                Font:    font.Font{Typeface: "Liberation", Variant: "Serif", Size: 14, Weight: ifont.WeightBold},
                Handler: plot.DefaultTextHandler,
        }

        textStyleFlank := text.Style{
                Color:   color.Black,
                Font:    font.Font{Typeface: "Liberation", Variant: "Serif", Size: 14, Weight: ifont.WeightNormal},
                Handler: plot.DefaultTextHandler,
        }

        //dcc := gn_draw.New(canvas)
        var dcc gn_draw.Canvas
        var canvas *vgimg.Canvas

        img_Y := 0
        if len(plotsGene) >0 {
            canvasHeight := canvasHeight + geneHeight+Ydistance
            canvas = vgimg.NewWith(
                                    vgimg.UseWH(canvasWith, canvasHeight),
                                    vgimg.UseDPI(dpi),
                                   )
            dcc = gn_draw.New(canvas)
            dcc, canvas = CombineGene(plotsGene, dcc, canvas,
                                      geneWith, geneHeight,
                                      textStyleFlank, dpi,
                                      samples)
            img_Y = img_Y + 1
        }else{
            canvas = vgimg.NewWith(
                                    vgimg.UseWH(canvasWith, canvasHeight),
                                    vgimg.UseDPI(dpi),
                                   )
            dcc = gn_draw.New(canvas)
        }

        geneH := vg.Length(img_Y)*geneHeight+Ydistance
        dcc, canvas = CombineSV(plots, dcc, canvas,
                                imgwith, imgheight, geneH,
                                Xdistance, Ydistance,
                                textStyleFlank, textStyleSample,
                                xlabel, infos, flankSize,
                                dpi, img_Y, flanksID, samples)

        dcc, canvas = CombineGT(plotsGT, dcc, canvas,
                                hpWith, imgwith, imgheight, geneH,
                                Xdistance, Ydistance,
                                textStyleFlank,
                                flankSize, dpi, img_Y, samples, samplesGT)

        file, err := os.OpenFile(filename, os.O_WRONLY|os.O_CREATE|os.O_TRUNC, 0644)
        if err != nil {
                log.Println("Error creating file:", err)
        }
        defer file.Close()

        fileParts := strings.Split(filename, ".")
        imgformat := fileParts[len(fileParts)-1]

        switch imgformat {
        case "png":
                var png vgimg.PngCanvas = vgimg.PngCanvas{Canvas: canvas}
                if _, err := png.WriteTo(file); err != nil {
                        log.Println("Error saving image:", err)
                }
        case "pdf":
                pdfcanvas := vgpdf.New(canvasWith, canvasHeight)
                img := vgimg.PngCanvas{Canvas: canvas}
                flippedImage := flipImageXAxis(img.Image())
                rect := vg.Rectangle{
                        Min: vg.Point{X: vg.Length(0), Y: vg.Length(0)}, //w, h := img.Bounds().Size().X, img.Bounds().Size().Y
                        Max: vg.Point{X: dcc.Size().X, Y: dcc.Size().Y}, //dcc.Size().X vg.Length(0)
                }
                pdfcanvas.DrawImage(rect, flippedImage)
                if _, err := pdfcanvas.WriteTo(file); err != nil {
                        log.Println("Error saving image:", err)
                }
        case "svg":
                svgcanvas := vgsvg.New(canvasWith, canvasHeight)
                img := vgimg.PngCanvas{Canvas: canvas}
                rect := vg.Rectangle{
                        Min: vg.Point{X: vg.Length(0), Y: vg.Length(0)}, //w, h := img.Bounds().Size().X, img.Bounds().Size().Y
                        Max: vg.Point{X: dcc.Size().X, Y: dcc.Size().Y}, //dcc.Size().X vg.Length(0)
                }
                svgcanvas.DrawImage(rect, img.Image())
                if _, err := svgcanvas.WriteTo(file); err != nil {
                        log.Println("Error saving image:", err)
                }
        case "jpeg":
                var jpeg vgimg.JpegCanvas = vgimg.JpegCanvas{Canvas: canvas}
                if _, err := jpeg.WriteTo(file); err != nil {
                        log.Println("Error saving image:", err)
                }
        case "tiff":
                var tiff vgimg.TiffCanvas = vgimg.TiffCanvas{Canvas: canvas}
                if _, err := tiff.WriteTo(file); err != nil {
                        log.Println("Error saving image:", err)
                }
        }
}

func randomData(n [][]int) plotter.XYs {
        pts := make(plotter.XYs, len(n))
        for i := 0; i < len(n); i++ {
                x := float64(n[i][0])
                y := float64(n[i][1])
                pts[i].X = x
                pts[i].Y = y
        }
        return pts
}

// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
func dealHP(flankSignals map[int][][]int, start int, stop int)(int, int){
    hp1, hp2 := 0, 0
    start_tmp, stop_tmp := float64(start-20), float64(stop+20)
    half_size := (stop_tmp - start_tmp)/2.0
    for _, signals := range flankSignals{
        hg_flag := -1
        label1:
        for _, signal := range signals{
            x1 := float64(signal[0])
            x2 := float64(signal[1])
            //fmt.Println("signal[2]:", signal)
            if signal[2] >= 500{
                if (x1 <= start_tmp && start_tmp <= x2) || (x1 <= stop_tmp && stop_tmp <= x2){
                    hg_flag = 0
                    continue
                }
            }

            switch signal[2]{
            case 1:
                if start_tmp + half_size <= x2 && x2  <= stop_tmp {
                    hg_flag = 1
                    break label1
                }
            case 2:
                if (start_tmp <= x1 && x1 <= stop_tmp) || (start_tmp <= x2 && x2 <= stop_tmp) {
                    hg_flag = 1
                    break label1
                }
            case 3:
                if (start_tmp <= x1 && x1 <= stop_tmp) || (start_tmp <= x2 && x2 <= stop_tmp) {
                    hg_flag = 1
                    break label1
                }
            case 4:
                if start_tmp <= x1 && x1  <= stop_tmp - half_size {
                    hg_flag = 1
                    break label1
                }
            }
        }
        if hg_flag == 1{
            hp2 = hp2 + 1
        }else if hg_flag ==0{
            hp1 = hp1 + 1
        }
        //fmt.Println("num:",num, signals)
    }
    return hp1, hp2
}

func figData(allSignal map[int]map[int][][]int, flanks []int,
             start int, stop int) (map[int]map[int][][][]int, int, int) {
        figDat := make(map[int]map[int][][][]int)
        // 外部转换和排序

        hp1, hp2 := 0, 0
        for flank_i, flank := range flanks {
                reads := allSignal[flank]
                if flank_i == 0{
                    hp1, hp2 = dealHP(reads, start, stop)
                    fmt.Println("sample:", hp1, hp2)
                }
                nu := 1
                nums := make([]int, len(reads))
                i := 0
                for key := range reads {
                        nums[i] = key
                        i++
                }
                sort.Ints(nums)
                flankData := make(map[int][][][]int)
                for _, num := range nums {
                        signals := reads[num]
                        for _, signal := range signals {
                                colorId := signal[2]
                                mpQ := -1

                                if colorId >= 500{
                                    mpQ = colorId % 500
                                    colorId = 5
                                }
                                x1 := signal[0]
                                x2 := signal[1]
                                y := nu
                                colorData, exists := flankData[colorId]

                                if exists {
                                        colorData = append(colorData, [][]int{{x1, y}, {x2, y}, {mpQ}})
                                } else {
                                        colorData = [][][]int{{{x1, y}, {x2, y}, {mpQ}}}
                                }
                                flankData[colorId] = colorData
                        }
                        nu += 1
                }
                flankData[6] = [][][]int{{{start, 1}, {start, nu}, {-1}}, {{stop, 1}, {stop, nu}, {-1}}}
                /*[[start ,1],[start, nu]],[[stop ,1],[stop, nu]]
                [start ,1] x, y
                [start, nu] x, y

                */
                figDat[flank] = flankData
        }
        return figDat, hp1, hp2
}

func DealPosition(rawSignals string) ([][]int, int, int, int) {
        signals := make([][]int,0)
        readSignals := strings.Split(strings.TrimRight(rawSignals, "\t"), "\t")
        //fmt.Println("readSignals:",readSignals)
        //fmt.Println("rawSignals:",rawSignals)
        //chr1    10000   15290   500     1       1       c1h1fh_5481b6ac-a312-0216-af8f-c39ec2d6c82e     8865,10000,1    10360,10365,2   10809,10820,2   10851,10858,2   14455,14460,3,CTCTTAAGAG
        alStart, _  := strconv.Atoi(readSignals[1])
        alStop , _  := strconv.Atoi(readSignals[2])
        mapQ   , _  := strconv.Atoi(readSignals[3])
        //fmt.Println("readSignals:",readSignals)
        for _, sig := range readSignals[7:]{
            sigs := strings.Split(sig, ",")
            start, _ := strconv.Atoi(sigs[0])
            end  , _ := strconv.Atoi(sigs[1])
            typ  , _ := strconv.Atoi(sigs[2])
            newArray := []int{start, end, typ}
            signals = append(signals, newArray)
        }
        return signals, alStart, alStop, mapQ
}

func KeyInMap(allSignals map[int]map[int][][]int, flank int, nums int, tmp []int) map[int]map[int][][]int {
        if _, ok1 := allSignals[flank]; ok1 {
                if _, ok2 := allSignals[flank][nums]; ok2 {
                        allSignals[flank][nums] = append(allSignals[flank][nums], tmp)
                } else {
                        allSignals[flank][nums] = [][]int{tmp}
                }
        } else {
                allSignals[flank] = map[int][][]int{nums: {tmp}}
        }
        return allSignals
}

func caluFlank(x1 float64, x2 float64, minPos float64, maxPos float64)(int, int){
    start := math.Max(minPos, x1)
    end := math.Min(maxPos, x2)
    return int(start), int(end)
}
func CutFlank(tmp []int, minSignal int, minPos int, maxPos int, flank int, nums int) ([]int, bool) {
        x1 := tmp[0]
        x2 := tmp[1]
        flag := tmp[2]
        tmp1 := make([]int, 5)
        var start int
        var end int
        //a b minPos maxPos
        //c d x1 x2
        // Check if the intervals intersect
        if maxPos < x1 || x2 < minPos {
                    return tmp1, false
        }else if flag == 1{
            if x2 <= maxPos && x2 >= minPos{
                start, end = caluFlank(float64(x1), float64(x2), float64(minPos), float64(maxPos))
            }else{
                return tmp1, false
            }
        }else if flag == 4 {
            if x1 <= maxPos && x1 >= minPos{
                start, end = caluFlank(float64(x1), float64(x2), float64(minPos), float64(maxPos))
            }else{
                return tmp1, false
            }
        }else{
            start, end = caluFlank(float64(x1), float64(x2), float64(minPos), float64(maxPos))

        }
 
        if Abs(end - start) < minSignal{
            return tmp1, false
        }

        tmp1[0] = start
        tmp1[1] = end
        tmp1[2] = flag
        tmp1[3] = flank
        tmp1[4] = nums
        return tmp1, true
}

func CutFlankDeal(fStart int, fStop int, ailgnStart int, alignStop int, minSignal int, flank int, nums int, mapQ int, signals [][]int, allSignals map[int]map[int][][]int) map[int]map[int][][]int {
        for i := 0; i < len(signals); i++ {
                //fmt.Println("signals:",signals)
                tmp, flg := CutFlank(signals[i], minSignal, fStart, fStop, flank, nums)
                if flg {
                        allSignals = KeyInMap(allSignals, flank, nums, tmp)
                }
        }
        ailgnPosition := []int{ailgnStart, alignStop, mapQ}
        tmp1, flg := CutFlank(ailgnPosition, minSignal, fStart, fStop, flank, nums)
        if flg {
                allSignals = KeyInMap(allSignals, flank, nums, tmp1)
        }
        return allSignals
}

func Abs(n int) int {
     if n < 0 {
         return -n
     }
     return n
}

func dealReads(rawSignals string, min_quality int, minSignal int, flanks []int, flanktofStart map[int]int, flanktofStop map[int]int, num int, flanksInvalid map[int]int, allSignals map[int]map[int][][]int) (map[int]map[int][][]int, map[int]int) {
        // 将结果发送到通道
        signals, alStart, alStop, mapQ := DealPosition(rawSignals)
        if mapQ < min_quality{
            return allSignals, flanksInvalid
        }

        rawSiz := len(allSignals)
        if len(signals)>0{
            for _, flank := range flanks {
                    fStart := flanktofStart[flank]
                    fStop  := flanktofStop[flank]
                    allSignals = CutFlankDeal(fStart, fStop, alStart, alStop, minSignal, flank, num, mapQ, signals, allSignals)
                    newSize := len(allSignals)
                    if newSize > rawSiz {
                            delete(flanksInvalid, flank)
                        }
                    rawSiz = newSize
            }
        }

        return allSignals, flanksInvalid
}

func flanktoids(flanks []int) []string {
        arraySize := len(flanks)
        keytoflanks := make(map[int]string, arraySize)
        for _, k := range flanks {
                var idx string
                if k > 1000 {
                        idx = strconv.FormatFloat(float64(k)/1000, 'f', -2, 64) + "kb"
                } else {
                        idx = strconv.Itoa(k) + "bp"
                }
                keytoflanks[k] = idx

        }
        sort.Ints(flanks)
        flanksID := make([]string, arraySize)
        for i, k := range flanks {
                flanksID[i] = keytoflanks[k]
        }
        return flanksID
}

func stringtomap(str_values []string, str_keys []string) map[string]string {
        //str_values := strings.Split(str_value, ",")
        data := make(map[string]string)
        for key, value := range str_values {
                //values := strings.Split(value, ":")
                data[str_keys[key]] = value
        }
        return data
}

func stringtointmap(flanks []string, startEnd string) map[int]int {
        data := make(map[int]int)
        startEnd_values := strings.Split(startEnd, ",")
        //fmt.Println("startEnd_values:", startEnd_values)
        for key, value := range flanks{
            value_int, _ := strconv.Atoi(value)
            startEnd_tmp, _ := strconv.Atoi(startEnd_values[key])
            data[value_int] = startEnd_tmp
        }
        return data
}

//export DealWith
func DealWith(
              c_rawSignals *C.char,
              c_chrom *C.char,
              c_svtype *C.char,
              Rflanks *C.char,
              c_fStart *C.char, 
              c_fStop *C.char,
              c_dpi *C.int,
              c_imgwith *C.int,
              c_imgheight *C.int,
              c_min_quality *C.int,
              c_minSignal *C.int,
              c_start *C.int,
              c_stop *C.int,
              c_outpathimg *C.char,
              c_samples *C.char,
              c_gt *C.char, 
              c_svInfo *C.char, 
              c_geneIDstr *C.char,
              c_tranIDstr *C.char,
              c_strandstr *C.char,
              c_genePos_str *C.char,
              c_tranPos_str *C.char,
              c_exonPos_str *C.char,
              c_fivePos_str *C.char,
              c_threePos_str *C.char,
              c_snp_str *C.char,
              legendflag bool,){
    startTime := time.Now()			  
    ///////////////////////////////////////////////////		 
    geneIDstr := C.GoString(c_geneIDstr)
    tranIDstr := C.GoString(c_tranIDstr)
    strandstr := C.GoString(c_strandstr)
    genePos_str := C.GoString(c_genePos_str)
    tranPos_str := C.GoString(c_tranPos_str)
    exonPos_str := C.GoString(c_exonPos_str)
    fivePos_str := C.GoString(c_fivePos_str)
    threePos_str := C.GoString(c_threePos_str)
    snp_str := C.GoString(c_snp_str)

    allGenesInfo := dealGenesInfo(geneIDstr,
                                  tranIDstr,
                                  strandstr,
                                  genePos_str,
                                  tranPos_str,
                                  exonPos_str,
                                  fivePos_str,
                                  threePos_str)
    fmt.Println("---------------------------------------------------")
    samplePlotGene := make([]*plot.Plot, 0)
    gene_plot, gene_h, gene_flag := plotGeneSnp(allGenesInfo, snp_str)
    if gene_flag{
        samplePlotGene = append(samplePlotGene, gene_plot)
    }
		
    //////////////////////////////////////////////////////
    chrom := C.GoString(c_chrom)
    svInfo := C.GoString(c_svInfo)
    rfls := strings.Split(C.GoString(Rflanks), ",")
    flanktofStart := stringtointmap(rfls, C.GoString(c_fStart))
    flanktofStop := stringtointmap(rfls, C.GoString(c_fStop))
    samplesGT := make(map[string]string)
    samples := strings.Split(C.GoString(c_samples), ",")
    gts := strings.Split(C.GoString(c_gt), ",")
    for i, sample := range samples {
            samplesGT[sample] = gts[i]
    }

    rawSignals := strings.Split(C.GoString(c_rawSignals), ":")
    signalsmap := stringtomap(rawSignals, samples)

    svtype := C.GoString(c_svtype)
    start := int(*c_start)
    stop := int(*c_stop)
    dpi  := int(*c_dpi)
    min_quality := int(*c_min_quality) + 500
    minSignal := int(*c_minSignal)
    outpathimg := C.GoString(c_outpathimg)
    //fmt.Println(chrom, samples, svtype, start, stop, outpathimg)
    flankSize := len(rfls)
    flanks := make([]int, flankSize)
    for i, str := range rfls {
            num, _ := strconv.Atoi(str)
            flanks[i] = num
    }
    sort.Ints(flanks)
    ///////////////////////////////////////////////////
    samplesplots := make(map[string][]*plot.Plot)
    samplePlotGT  := make(map[string]*plot.Plot)
    for sample, signal := range signalsmap {
        signals := strings.Split(signal, "|")
        flanksInvalid := map[int]int{}
        allSignals := make(map[int]map[int][][]int)
        for _, fk := range flanks {
                flanksInvalid[fk] = fk
        }
        for i, signal := range signals {
            if signal != "" {
                allSignals, flanksInvalid = dealReads(signal, min_quality, minSignal, flanks, flanktofStart, flanktofStop, i+1, flanksInvalid, allSignals)
            }
        }
		
		sampleplots, pBar := createPlots(
                                       allSignals,
                                       chrom,
                                       svtype,
                                       start,
                                       stop,
                                       flanks,
                                       flanksInvalid,
                                       flankSize,
                                       legendflag,
                                      )

		samplePlotGT[sample] = pBar
		samplesplots[sample] = sampleplots
    }
		
    ///////////////////////////////////////////////////	
    fmt.Println(outpathimg)
    imgwith := int(*c_imgwith) //4
    imgheight := int(*c_imgheight) //4
    //flanksize := 5
    var xpad float32 = 0.0
    var ypad float32 = 0.3
    xlabel := chrom 
    flanksID := flanktoids(flanks)
    CombineAndSave(
                   samplesplots,
                   samplePlotGT,
                   samplePlotGene,
                   gene_h,
                   gene_plot,
                   gene_flag,
                   outpathimg,
                   imgwith,
                   imgheight,
                   flankSize,
                   xpad,
                   ypad,
                   xlabel,
                   flanksID,
                   samples,
                   samplesGT,
                   svInfo,
                   dpi,
                  )
    // 记录函数结束时间
    endTime := time.Now()

    // 计算函数运行时间
    executionTime := endTime.Sub(startTime)

    fmt.Printf("函数运行时间：%v\n", executionTime)
}
func main() {}

