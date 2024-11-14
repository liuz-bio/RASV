package main

import (
        "C"
        //"fmt"
        "regexp"
        "strconv"
        "strings"
        //"sync"
        "compress/gzip"
        //"encoding/json"
        "unicode/utf8"
        "os"
)

func SplitString(input string) ([]int, []string) {
        //use regular expressions to match numbers and letters
        re := regexp.MustCompile(`(\d+)([A-Za-z])`)
        matches := re.FindAllStringSubmatch(input, -1)

        //results of storing numbers and letters respectively
        numbers := make([]int, len(matches))
        letters := make([]string, len(matches))

        //extract numbers and letters and store them in the corresponding array
        for i := range matches {
                match := matches[i]
                number, _ := strconv.Atoi(match[1])
                letter := match[2]
                numbers[i] = number
                letters[i] = letter
        }

        return numbers, letters
}

func cutRead(query_start int, seq string, sig_len int) (string){
    seq_len := utf8.RuneCountInString(seq)
    tmp_s := query_start
    tmp_e := query_start+sig_len
    if tmp_s < 0 {
        tmp_s = 0
    }
    if tmp_e > seq_len {
        tmp_e = seq_len
    }
    return seq[tmp_s:tmp_e]
}


func DealPosition(numbers []int, letters []string, query_seq string, start int, minSignal int) []string {
        refs_start := start
        query_start := 0
        signals := make([]string, 0)
        if strings.EqualFold(letters[0], "S") { //1 
                query_start = query_start + numbers[0]
                if numbers[0] >= minSignal{
                    tmp_start := refs_start - numbers[0]
                    if tmp_start > 0 {
                        newArray := []string{strconv.Itoa(tmp_start), strconv.Itoa(refs_start), strconv.Itoa(1), strconv.Itoa(query_start)}
                        //seq_s1, quality_s1 := cutRead(query_start, query_seq, query_quality, 150,  seq_len)
                        signals = append(signals,strings.Join(newArray[:], ","))
                    }
                }
        }
        //letters_Size := len(letters)
        for i := range letters{
        //for i := 0; i < letters_Size; i++ {
                length := numbers[i]
                if strings.EqualFold(letters[i], "M") {
                        query_start = query_start + length
                        refs_start = refs_start + length
                } else if strings.EqualFold(letters[i], "D") { //2
                        refs_stop := refs_start + length
                        if length >= minSignal {
                                newArray := []string{strconv.Itoa(refs_start), strconv.Itoa(refs_stop), strconv.Itoa(2)}
                                signals = append(signals, strings.Join(newArray[:], ","))
                                //seq_d, quality_d := cutRead(query_start, query_seq, query_quality, length*2, seq_len)
                        }
                        refs_start = refs_stop
                } else if strings.EqualFold(letters[i], "I") { //3
                        reads_stop := refs_start + length
                        if length >= minSignal {
                                //newArray := []string{strconv.Itoa(refs_start), strconv.Itoa(reads_stop), strconv.Itoa(3)}
                                //signals = append(signals, strings.Join(newArray[:], ","))
                                seq_i  := cutRead(query_start, query_seq, length)
                                newArray := []string{strconv.Itoa(refs_start), strconv.Itoa(reads_stop), strconv.Itoa(3), seq_i}
                                signals = append(signals, strings.Join(newArray[:], ","))
                        }
                        query_start = query_start + length

                }
        }

        last_index := len(letters)-1
        if strings.EqualFold(letters[last_index], "S") { //4
                s_len := numbers[last_index]
                if s_len >= minSignal {
                    newArray := []string{strconv.Itoa(refs_start), strconv.Itoa(refs_start + s_len), strconv.Itoa(4), strconv.Itoa(query_start)}
                    signals = append(signals, strings.Join(newArray[:], ","))
                    //seq_s2, quality_s2 := cutRead(query_start, query_seq, query_quality, 150, seq_len)
                }
        }
        return signals
}


func processReads(input string, 
                  query_name string, 
                  mapQ string, 
                  read_is_reverse string, 
                  read_is_secondary string, 
                  query_seq string,
                  alStart int, 
                  alStop int,
                  queryAlStart int,
                  queryAlStop  int,
                  queryLen int,
                  minSignal int)(string) {
        // 将结果发送到通道
        numbers, letters := SplitString(input)
        signals := DealPosition(numbers, letters, query_seq, alStart, minSignal)
        var result string
        if len(signals)>0 {
            result =  strconv.Itoa(alStart)+"\t"+strconv.Itoa(alStop)+"\t"+strconv.Itoa(queryAlStart)+"\t"+strconv.Itoa(queryAlStop)+"\t"+strconv.Itoa(queryLen)+"\t"+mapQ+"\t"+read_is_reverse+"\t"+read_is_secondary+"\t"+query_name+"\t"+strings.Join(signals[:], "\t")
        }else{
            result =  strconv.Itoa(alStart)+"\t"+strconv.Itoa(alStop)+"\t"+strconv.Itoa(queryAlStart)+"\t"+strconv.Itoa(queryAlStop)+"\t"+strconv.Itoa(queryLen)+"\t"+mapQ+"\t"+read_is_reverse+"\t"+read_is_secondary+"\t"+query_name
        }
        return result
}

func writerData(allSignals []string, outPut string){
    file, _ := os.Create(outPut)
    defer file.Close()

    //for _, signal := range allSignals {
    for i := range allSignals {
        signal := allSignals[i]
        if signal != "" {
            _, _ = file.WriteString(signal + "\n")
        }
    }
}

func writerData1(allSignals []string, outPut string){
    gzipFile, err := os.Create(outPut)
    if err != nil {
        panic(err)
    }
    defer gzipFile.Close()

    gzipWriter := gzip.NewWriter(gzipFile)
    defer gzipWriter.Close()
  
    //for _, signal := range allSignals {
    for i := range allSignals {
        signal := allSignals[i]
        if signal != "" {
            _, _ = gzipWriter.Write([]byte(signal + "\n"))
        }
    }
}


//','.join(cigs).encode(), c_start, c_stop, outPaths.encode(), sample.encode(),
//                         prefix.encode(), ','.join(map(str, alignStarts)).encode(),
//                         ','.join(map(str, alignStops)).encode()
//export DealWith
func DealWith(cigas *C.char, 
              c_outpaths *C.char, 
              c_query_names *C.char,
              c_alignStart *C.char, 
              c_alignStop *C.char, 
              c_queryAlignStart *C.char,
              c_queryAlignStop *C.char,
              c_queryLengths *C.char,
              c_map_quality *C.char,             
              c_query_is_reverse *C.char,
              c_query_is_secondary *C.char,
              c_query_seqs *C.char,
              c_minSignal *C.int) *C.char {

        inputs := strings.Split(C.GoString(cigas), ",")
        outpaths := C.GoString(c_outpaths)
        alignStart := strings.Split(C.GoString(c_alignStart), ",")
        alignStop := strings.Split(C.GoString(c_alignStop), ",")
        queryAlignStart := strings.Split(C.GoString(c_queryAlignStart), ",")
        queryAlignStop := strings.Split(C.GoString(c_queryAlignStop), ",")
        queryLengths := strings.Split(C.GoString(c_queryLengths), ",")
        query_names := strings.Split(C.GoString(c_query_names), "|")
        map_quality := strings.Split(C.GoString(c_map_quality), ",")
        query_is_reverse := strings.Split(C.GoString(c_query_is_reverse), ",")
        query_is_secondary := strings.Split(C.GoString(c_query_is_secondary), ",")
        query_seqs := strings.Split(C.GoString(c_query_seqs), ",")
        minSignal := int(*c_minSignal)

        allSignals := make([]string,len(inputs))
        // 并行处理每个项目
        //for i, input := range inputs {
        for i := range inputs {
                input := inputs[i]
                alStart, _ := strconv.Atoi(alignStart[i])
                alStop, _ := strconv.Atoi(alignStop[i])
                queryAlStart, _ := strconv.Atoi(queryAlignStart[i])
                queryAlStop, _  := strconv.Atoi(queryAlignStop[i])
                queryLen, _ := strconv.Atoi(queryLengths[i])
                result := processReads(input, 
                                       query_names[i], 
                                       map_quality[i], 
                                       query_is_reverse[i], 
                                       query_is_secondary[i],  
                                       query_seqs[i],
                                       alStart, 
                                       alStop, 
                                       queryAlStart,
                                       queryAlStop,
                                       queryLen,
                                       minSignal)
                allSignals = append(allSignals, result)
        }

        writerData(allSignals, outpaths)
        return C.CString("asasadsasdasd")
        //C.CString(strings.Join(allFlankSignals[:], ","))
}

func main() {}
