go mod init golang
go mod tidy

go build -o extractSignal.so -buildmode=c-shared extractSignal.go
go build -o plotSample.so -buildmode=c-shared plotSample.go
go build -o plotVcf.so -buildmode=c-shared plotVcf.go
go build -o cigar.so -buildmode=c-shared cigar.go
