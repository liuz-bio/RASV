package main

import (
	"fmt"
	"log"
	"os"
        //"io"
        "github.com/biogo/hts/fai"
)


func faidx(m *fai.File, chrom string, start int, stop int)(string){
        s, err := m.SeqRange(chrom, start, stop)
        buf := make([]byte, stop-start)
        var got []byte
        n, err := s.Read(buf)
        got = append(got, buf[:n]...)
        if err != nil {
            log.Fatal("unexpected error: %v", err)
        }
        seq := string(got)
        s.Reset()
        return seq
    }

func genomeIndex(genomeFile string) (*fai.File, *os.File){//io.ReaderAt){
        //genomeFile := "/home/lz/work_space/Database/VISOR_test/Genome_old/Genome.fa"

        f, err := os.Open(genomeFile + ".fai")
        if err != nil {
                log.Fatal("unexpected error: %v", err)
        }

        idx, err := fai.ReadFrom(f)
        f.Close()
        if err != nil {
                log.Fatal("unexpected error: %v", err)
        }

        f, err = os.Open(genomeFile)
        if err != nil {
                log.Fatal("unexpected error: %v", err)
        }

        m := fai.NewFile(f, idx)
        //f.Close()
        fmt.Println("seq:", "##########################")
        return m, f
        }

func main() {
        genomeFile := "/home/lz/work_space/Database/VISOR_test/Genome_old/Genome.fa"
        genome, fasta := genomeIndex(genomeFile)
        chrom := "chr1"
        start, stop :=  0, 100
        seq := faidx(genome, chrom, start, stop)
        fmt.Println("seq1:", seq)        

        chrom = "chr1"
        start, stop = 1, 100
        seq = faidx(genome, chrom, start, stop)
        fmt.Println("seq2:", seq)
        fasta.Close()
}

/*
func main() {
	// 替换为你的基因组文件路径
	genomeFile := "/home/lz/work_space/Database/VISOR_test/Genome_old/Genome.fa"

        f, err := os.Open(genomeFile + ".fai")
	if err != nil {
		log.Fatal("unexpected error: %v", err)
	}

	idx, err := fai.ReadFrom(f)
	f.Close()
	if err != nil {
		log.Fatal("unexpected error: %v", err)
	}

        f, err = os.Open(genomeFile)
	if err != nil {
	        log.Fatal("unexpected error: %v", err)
	}
	m := fai.NewFile(f, idx)

        s, err := m.SeqRange("chr1", 1234567, 1234667)
        buf := make([]byte, 1234667-1234567)
        var got []byte
        n, err := s.Read(buf)
        got = append(got, buf[:n]...)
        if err != nil {
            log.Fatal("unexpected error: %v", err)
        }
        a := string(got)
        //startIndex := len("%!(EXTRA string=")
        //endIndex := len(a) - 1
        fmt.Println("seq:", a, got)
        s.Reset()
}
*/

