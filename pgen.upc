#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <time.h>
#include <upc.h>

#include "packingDNAseq.h"
#include "kmer_hash.upc"

// typedef shared [] hash_table_t* shtptr;
// shared shtptr all_hash_tables[THREADS];

int main(int argc, char *argv[]){

	/** Declarations **/
	double inputTime=0.0, constrTime=0.0, traversalTime=0.0;

	char cur_contig[MAXIMUM_CONTIG_SIZE], unpackedKmer[KMER_LENGTH+1], left_ext, right_ext, *input_UFX_name;
	int64_t posInContig, contigID = 0, totBases = 0, ptr = 0, nKmers, cur_chars_read, total_chars_to_read;
	unpackedKmer[KMER_LENGTH] = '\0';
	kmer_t *cur_kmer_ptr;
	start_kmer_t *startKmersList = NULL, *curStartNode;
	unsigned char *working_buffer;
	FILE *inputFile, *serialOutputFile;
	/* Read the input file name */
	input_UFX_name = argv[1];
	/* Initialize lookup table that will be used for the DNA packing routines */
	init_LookupTable();
	nKmers = getNumKmersInUFX(input_UFX_name);

	upc_lock_t *l;
	l = upc_all_lock_alloc(); 

	/** Read input **/
	upc_barrier;
	inputTime -= gettime();
	///////////////////////////////////////////
	// Your code for input file reading here //
	int64_t kmers_per_proc = (nKmers+THREADS-1)/THREADS;
	int64_t mykmers = kmers_per_proc;
	if (THREADS - 1 == MYTHREAD) mykmers = nKmers - (THREADS-1)*kmers_per_proc;
	total_chars_to_read = mykmers * LINE_SIZE;
	working_buffer = (unsigned char*) malloc(total_chars_to_read * sizeof(unsigned char));
	inputFile = fopen(input_UFX_name, "r");
	fseek(inputFile, MYTHREAD*kmers_per_proc*LINE_SIZE*sizeof(unsigned char), SEEK_SET);
	cur_chars_read = fread(working_buffer, sizeof(unsigned char), total_chars_to_read , inputFile);
	// upc_lock(l);
	// printf("Thread %d of %d: hello UPC world, nKmers is %d, line size is %d, read size is %d, Buffer is \n%*s\n", MYTHREAD, THREADS, nKmers, LINE_SIZE, sizeof(unsigned char)* total_chars_to_read , sizeof(unsigned char)* total_chars_to_read ,working_buffer); 
	// upc_unlock(l);
	///////////////////////////////////////////
	upc_barrier;
	inputTime += gettime();

	/** Graph construction **/
	constrTime -= gettime();
	///////////////////////////////////////////
	// Your code for graph construction here //
	shared memory_heap_t* memory_heap;
	shared hash_table_t* hashtable;
	hashtable = (shared hash_table_t*) create_hash_table(nKmers, memory_heap);

	upc_barrier;
	printf("yolo0.0, size = %d\n", hashtable->size);
	// printf("yolo0.1, size = %d\n", all_hash_tables[(MYTHREAD+1)%4]->size);
	add_kmer(hashtable, memory_heap, &working_buffer[0], working_buffer[KMER_LENGTH+1], working_buffer[KMER_LENGTH+2], l);
	upc_barrier;
	if (MYTHREAD == 0) {
		for (int i = 25; i < 26;i++) {//hashtable->size; i++) {
			printf("%d\n", KMER_PACKED_LENGTH, hashtable->table[i].head->kmer);
		}
	}

	// upc_barrier;
	// kmer_t * curmer = lookup_kmer(all_hash_tables[MYTHREAD], "CACAAAGTCAGCTGTGCTC");
	// printf("THREAD: %d, curmer is %c\n", MYTHREAD, curmer->r_ext);
	// curmer = lookup_kmer(all_hash_tables[(MYTHREAD+1)%4], "CACAAAGTCAGCTGTGCTC");
	// printf("THREAD: %d, curmer is %c\n", MYTHREAD, curmer->r_ext);
	// all_hash_tables[MYTHREAD] = (hash_table_t*)create_hash_table(mykmers, &memory_heap);
	///////////////////////////////////////////
	upc_barrier;
	constrTime += gettime();

	/** Graph traversal **/
	traversalTime -= gettime();
	////////////////////////////////////////////////////////////
	// Your code for graph traversal and output printing here //
	// Save your output to "pgen.out"                         //
	////////////////////////////////////////////////////////////
	upc_barrier;
	traversalTime += gettime();

	/** Print timing and output info **/
	/***** DO NOT CHANGE THIS PART ****/
	if(MYTHREAD==0){
		printf("%s: Input set: %s\n", argv[0], argv[1]);
		printf("Number of UPC threads: %d\n", THREADS);
		printf("Input reading time: %f seconds\n", inputTime);
		printf("Graph construction time: %f seconds\n", constrTime);
		printf("Graph traversal time: %f seconds\n", traversalTime);
	}
	return 0;
}
