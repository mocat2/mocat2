all: fastq_trim_filter fastx_quality_stats

LDFLAGS = -static

fastx_quality_stats: faststats.cpp fastq.cpp
	g++ -g -O2 -Wall -o $@ faststats.cpp fastq.cpp 

fastq_trim_filter: fastrtf.cpp fastq.cpp popenstream.cpp
	g++ -g -O2 -Wall -o $@ fastrtf.cpp fastq.cpp popenstream.cpp $(CFLAGS) -L /usr/local/lib -lboost_program_options-mt -lz

clean:
	rm -f fastq_trim_filter fastx_quality_stats

.PHONY: all
