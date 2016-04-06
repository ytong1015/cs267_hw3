#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <time.h>
#include <upc.h>
#include <upc_collective.h>


#include "packingDNAseq.h"
#include "kmer_hash.upc"

typedef shared [] hash_table_t* shtptr;
shared shtptr all_hash_tables[THREADS];

 
 // Allocate an array of 'count' locks and cache it in private arrays.
 //
 // The actual array is replicated THREADS time in the shared heap, but
 // this function returns a pointer to a private piece of the array.
 //
 // The lock allocation calls are distributed over the threads as if
 // implemented by the following:
 //      upc_lock_t* shared [*] table[count];
 //      upc_forall(int i=0; i<count; ++i; &table[count])
 //        table[i] = upc_global_lock_alloc();
 // followed by local replication of the table.
 //
 // This code works for any non-zero value of 'count'.
 // This function must be called collectively.
 //
 upc_lock_t **allocate_lock_array(unsigned int count) {
   const unsigned int blksize = ((count + THREADS - 1) / THREADS); // Round up for "[*]"
   const unsigned int padded_count = blksize * THREADS;
   upc_lock_t* shared *tmp = upc_all_alloc(padded_count, sizeof(upc_lock_t*));
   upc_lock_t* shared *table = upc_all_alloc(padded_count*THREADS, sizeof(upc_lock_t*));
 
   // Allocate lock pointers into a temporary array.
   // This code overlays an array of blocksize [*] on the cyclic one.
   upc_lock_t** ptmp = (upc_lock_t**)(&tmp[MYTHREAD]); // Private array "slice"
   const int my_count = upc_affinitysize(count,blksize,MYTHREAD);
   for (int i=0; i<my_count; ++i) ptmp[i] = upc_global_lock_alloc();
 
   // Replicate the temporary array THREADS times into the shared table array.
   // IN_MYSYNC:   Since each thread generates its local portion of input.
   // OUT_ALLSYNC: Ensures upc_free() occurs only after tmp is unneeded.
   upc_all_gather_all(table, tmp, blksize * sizeof(upc_lock_t*),
                      UPC_IN_MYSYNC|UPC_OUT_ALLSYNC);
 
   if (!MYTHREAD) upc_free(tmp);  // Free the temporary array exactly once
 
   // Return a pointer-to-private for local piece of replicated table
   return (upc_lock_t**)(&table[MYTHREAD]);
 } 


int main(int argc, char *argv[]){

	/** Declarations **/
	double inputTime=0.0, constrTime=0.0, traversalTime=0.0;

	unsigned char cur_contig[MAXIMUM_CONTIG_SIZE+1], unpackedKmer[KMER_LENGTH+1], left_ext, right_ext, *input_UFX_name;
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

	/** Read input **/
	upc_barrier;
	inputTime -= gettime();
	///////////////////////////////////////////
	// Your code for input file reading here //
	int64_t kmers_per_proc = (nKmers+THREADS-1)/THREADS;
	int64_t tablesize = nKmers*LOAD_FACTOR;
	int64_t mykmers = kmers_per_proc;
	printf(" nkmers = %d mykmers = %d\n", nKmers, mykmers);
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
	printf("The next print statement won't print..... ");
	shared int64_t *next_index, *hash_table;
	shared kmer_t *memory_heap;
	create_hash_table(nKmers, &memory_heap, &next_index, &hash_table); 

	upc_barrier;

	int64_t k = MYTHREAD*kmers_per_proc;
	ptr = 0;

//	printf("The next print statement won't print..... ");
//	upc_lock_t *l = ( upc_lock_t*) upc_all_alloc(tablesize, sizeof(upc_lock_t*));
//	upc_barrier;
	upc_lock_t ** lock_array = allocate_lock_array(tablesize);
        upc_lock_t ** lock_next_index = allocate_lock_array(nKmers);
//	for (int64_t i = MYTHREAD; i < tablesize; i+=THREADS) {
//		lock_array[i] = upc_global_lock_alloc();
//	}
//	upc_barrier;

	while(ptr < cur_chars_read)
	{

		left_ext = working_buffer[ptr+KMER_LENGTH+1];
		right_ext = working_buffer[ptr+KMER_LENGTH+2];

		// if (k%10000 == 0)
		// printf("GETS TO THE FIRST ONE, BITCH %d\t%d\t%d\t%d\n", MYTHREAD, k, ptr, mykmers);
		// add kmer

		// upc_lock(l);
		add_kmer(hash_table, memory_heap, &working_buffer[ptr],left_ext, right_ext, next_index, k, tablesize, lock_array, lock_next_index);
		// upc_unlock(l);
		ptr += LINE_SIZE;
		k++;

	}

	printf("see i told you it got stuck\n");
	///////////////////////////////////////////
	upc_barrier;
	constrTime += gettime();

	// Graph traversal //
	traversalTime -= gettime();
	////////////////////////////////////////////////////////////
	// Your code for graph traversal and output printing here //
	// Save your output to "pgen.out"                         //
	////////////////////////////////////////////////////////////
	char output_file_name[50];
	sprintf(output_file_name, "pgen%d.out",MYTHREAD);
	FILE *out_file = fopen(output_file_name, "w");
	int64_t i = 0;
	ptr = 0;

	for (; i<mykmers; i++, ptr += LINE_SIZE)
	{
		unsigned char left_ext = working_buffer[ptr+KMER_LENGTH+1];
		if (left_ext != 'F') continue;

		memcpy(cur_contig, &working_buffer[ptr], KMER_LENGTH*sizeof(unsigned char)) ;
		posInContig = KMER_LENGTH;
		right_ext = working_buffer[ptr+KMER_LENGTH+2];

		while (right_ext != 'F' && right_ext != 0)
		{

			cur_contig[posInContig] = right_ext;
			posInContig += 1;

			right_ext = lookup_kmer(memory_heap, hash_table, tablesize,&cur_contig[posInContig-KMER_LENGTH], next_index);

		}
		cur_contig[posInContig] = '\0';
		fprintf(out_file, "%s\n", cur_contig);


	}

	fclose(out_file);

	if (MYTHREAD == 0)
	{
		upc_free(memory_heap);
		upc_free(next_index);
		upc_free(hash_table);
	}






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
