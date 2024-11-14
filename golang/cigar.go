package main
import (
    "C"
    "strconv"
)

func cutRead(queryStart int, seq string, sigLen int) string {
    seqLen := len(seq)
    tmpS := queryStart
    tmpE := queryStart + sigLen
    if tmpS < 0 {
        tmpS = 0
    }
    if tmpE > seqLen {
        tmpE = seqLen
    }
    return seq[tmpS:tmpE]
}

func dealCigar(input string, querySeq string, start int, minSignal int) string {
    signals := make([]byte, 0, len(input)) // 预分配信号缓冲区

    refsStart := start
    queryStart := 0

    for i := 0; i < len(input); {
        number := 0
        for ; i < len(input) && input[i] >= '0' && input[i] <= '9'; i++ {
            number = number*10 + int(input[i]-'0')
        }
        letter := input[i]
        i++

        switch letter {
        case 'M':
            queryStart += number
            refsStart += number
        case 'I':
            readsStop := refsStart + number
            if number >= minSignal {
                appendSignal(&signals, refsStart, readsStop, queryStart, 3, querySeq)
            }
            queryStart += number
        case 'S':
            if refsStart == start {
                queryStart += number
                if number >= minSignal {
                    tmpStart := refsStart - number
                    if tmpStart > 0 {
                        appendSignal(&signals, tmpStart, refsStart, queryStart, 1, querySeq)
                    }
                }
            } else {
                if number >= minSignal {
                    appendSignal(&signals, refsStart, refsStart+number, queryStart, 4, querySeq)
                }
            }

	
        case 'D':
            refsStop := refsStart + number
            if number >= minSignal {
                appendSignal(&signals, refsStart, refsStop, queryStart, 2, "")
            }
            refsStart = refsStop
        }
    }

    return string(signals)
}

func appendSignal(signals *[]byte, start, stop, queryStart, signalType int, querySeq string) {
    *signals = strconv.AppendInt(*signals, int64(start), 10)
    *signals = append(*signals, ',')
    *signals = strconv.AppendInt(*signals, int64(stop), 10)
    *signals = append(*signals, ',')
    *signals = strconv.AppendInt(*signals, int64(signalType), 10)
    *signals = append(*signals, ',')
    if signalType == 3 {
        *signals = append(*signals, cutRead(queryStart, querySeq, stop-start)...)
	*signals = append(*signals, ',')
	*signals = strconv.AppendInt(*signals, int64(queryStart), 10)
    } else {
        *signals = strconv.AppendInt(*signals, int64(queryStart), 10)
    }
    *signals = append(*signals, '\t')
}

//export DealWith
func DealWith(cigar *C.char, alignStart *C.int, querySeq *C.char, minSignal *C.int) *C.char {
    start := int(*alignStart)
    min := int(*minSignal)
    query := C.GoString(querySeq)
    signals := dealCigar(C.GoString(cigar), query, start, min)
    return C.CString(signals)
}

func main() {}
