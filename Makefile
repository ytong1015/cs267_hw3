CC = gcc
UPCC = upcc

KMER_LENGTH 		= 19
KMER_PACKED_LENGTH 	= $(shell echo $$((($(KMER_LENGTH)+3)/4)))

# Add -std=gnu99 to CFLAGS if use gnu compiler
CFLAGS 	= -O3 -std=gnu99 -march=haswell
CFLAGSUPC = -O3 -std=gnu99 -march=haswell
DEFINE 	= -DKMER_LENGTH=$(KMER_LENGTH) -DKMER_PACKED_LENGTH=$(KMER_PACKED_LENGTH)
HEADERS	= contig_generation.h kmer_hash.h packingDNAseq.h contig_generation.upc kmer_hash.upc
LIBS	=

TARGETS	= serial pgen

all: 	$(TARGETS)

serial: serial.c $(HEADERS)
		$(CC) $(CFLAGS) -o $@ $< -DKMER_LENGTH=$(KMER_LENGTH) -DKMER_PACKED_LENGTH=$(KMER_PACKED_LENGTH) $(LIBS)

pgen:	pgen.upc $(HEADERS)
		$(UPCC) $(UPCFLAGS) -Wc,"$(CFLAGSUPC)" -o $@ $< $(DEFINE) $(LIBS)

clean :
	rm -f *.o
	rm -rf $(TARGETS)
