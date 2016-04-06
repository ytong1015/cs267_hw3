#ifndef KMER_HASH_H
#define KMER_HASH_H

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h> 
#include <string.h>
#include "contig_generation.upc"

/* Creates a hash table and (pre)allocates memory for the memory heap */
int64_t create_hash_table(int64_t nEntries, shared kmer_t **memory_heap, shared int64_t** next_pointers, shared int64_t** table_pointers)
{
   int64_t n_buckets = nEntries * LOAD_FACTOR;

   *table_pointers = (shared int64_t*) upc_all_alloc(n_buckets, sizeof(int64_t));
   // upc_memset(*table_pointers, DUMMY, n_buckets*sizeof(int64_t));
   upc_forall(int64_t i = 0; i < n_buckets; i++; i)
   {
         (*table_pointers)[i] = DUMMY;
   }

   *next_pointers = (shared int64_t*) upc_all_alloc(nEntries, sizeof(int64_t));
   // upc_memset(next_pointers, DUMMY, n_buckets*sizeof(int64_t));
   upc_forall(int64_t i = 0; i < nEntries; i++; i)
   {
         (*next_pointers)[i] = DUMMY;
   }
   
   if (*table_pointers == NULL) {
      fprintf(stderr, "ERROR: Could not allocate memory for the hash table: %lld buckets of %lu bytes\n", n_buckets, sizeof(int64_t));
      exit(1);
   }
   
   *memory_heap= (shared kmer_t *) upc_all_alloc(nEntries, sizeof(kmer_t));
   if (*memory_heap == NULL) {
      fprintf(stderr, "ERROR: Could not allocate memory for the heap!\n");
      exit(1);
   }
   return 0;
}

/* Auxiliary function for computing hash values */
int64_t hashseq(int64_t  hashtable_size, char *seq, int size)
{
   unsigned long hashval;
   hashval = 5381;
   for(int i = 0; i < size; i++) {
      hashval = seq[i] +  (hashval << 5) + hashval;
   }
   
   return hashval % hashtable_size;
}

/* Returns the hash value of a kmer */
int64_t hashkmer(int64_t  hashtable_size, char *seq)
{
   return hashseq(hashtable_size, seq, KMER_PACKED_LENGTH);
}

/* Looks up a kmer in the hash table and returns a pointer to that entry */
kmer_t lookup_kmer(shared kmer_t* memory_heap, shared int64_t* hash_table, int64_t tablesize, const unsigned char *kmer, shared int64_t* next_index)
{
   //printf("swag0\n");
   char packedKmer[KMER_PACKED_LENGTH];
   //printf("swag1\n");
   packSequence(kmer, (unsigned char*) packedKmer, KMER_LENGTH);
   //printf("swag2\n");
   int64_t hashval = hashkmer(tablesize, (char*) packedKmer);
   //printf("swag3\n");
   int64_t pos = hash_table[hashval];
   kmer_t temp;
   while (pos != DUMMY)
   {
      upc_memget(&temp, &memory_heap[pos], sizeof(kmer_t));
      if (memcmp(packedKmer, temp.kmer, KMER_PACKED_LENGTH*sizeof(char)) == 0) return temp;
      pos = next_index[pos];
   }  
   //printf("swag4 %d\n", hashval); 
   return temp;
}

/* Adds a kmer and its extensions in the hash tablknows e (note that a memory heap should be preallocated. ) */
int add_kmer( shared int64_t* hash_table, shared kmer_t* memory_heap, const unsigned char *kmer, char left_ext, char right_ext, shared int64_t* next_index, int64_t posInHeap, int tablesize, upc_lock_t ** lock_array, upc_lock_t ** lock_next_index)
{
   /* Pack a k-mer sequence appropriately */
   char packedKmer[KMER_PACKED_LENGTH];
   //printf("yolo0.1\n");
   packSequence(kmer, (unsigned char*) packedKmer, KMER_LENGTH);
   //printf("yolo0.2\n");
   int64_t hashval = hashkmer(tablesize, (char*) packedKmer);
   //printf("yolo0.3\n");
   kmer_t temp;
   //printf("yolo0.4\n");
   memcpy(temp.kmer, packedKmer, KMER_PACKED_LENGTH*sizeof(char));
   temp.l_ext = left_ext;
   temp.r_ext = right_ext;

   upc_memput(&memory_heap[posInHeap], &temp, sizeof(kmer_t));

   int64_t ptr, ptr2;
    upc_lock(lock_array[hashval]);
    ptr = hash_table[hashval];
    if (ptr == DUMMY)
    {
       hash_table[hashval] = posInHeap;
   // //    upc_unlock(lock_array[hashval]);
    //   return 0;
    }
    upc_unlock(lock_array[hashval]);
   while (ptr!=DUMMY)
   {
      upc_lock(lock_next_index[ptr]);
      ptr2 = next_index[ptr];
      if (ptr2 == DUMMY)
      {
         next_index[ptr] = posInHeap;
      // //    upc_unlock(lock_array[ptr]);
      //    break;
      }
      upc_unlock(lock_next_index[ptr]);
      ptr = ptr2;
   }
   //upc_unlock(lock_next_index[ptr]);
   return 0;
}


/* Deallocation functions */
int dealloc_heap(shared kmer_t *memory_heap)
{
   upc_free(memory_heap);
   return 0;
}

int dealloc_hashtable(shared int64_t *hashtable)
{
   upc_free(hashtable);
   return 0;
}

/* Adds a k-mer in the start list by using the memory heap (the k-mer was "just added" in the memory heap at position posInHeap - 1) */
void addKmerToStartList(shared kmer_t *memory_heap, start_kmer_t **startKmersList, int64_t posInHeap)
{
   start_kmer_t *new_entry;
   kmer_t *ptrToKmer;

//   int64_t prevPosInHeap = memory_heap->posInHeap - 1;
   //upc_memget(ptrToKmer, &(memory_heap[posInHeap]), sizeof(kmer_t));
   ptrToKmer = (kmer_t*) &(memory_heap[posInHeap]);
   new_entry = (start_kmer_t*) malloc(sizeof(start_kmer_t));
   new_entry->next = (*startKmersList);
   new_entry->kmerPtr = ptrToKmer;
   (*startKmersList) = new_entry;
}



#endif // KMER_HASH_H
