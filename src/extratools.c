#include "extratools.h"

//This file includes the implementations of all the extra tools adding to all steps
//e.g. Look Up Table
//     Hash Table
//     All things like those

void LoadLookupTable(LOOKUPTABLE * lookupTable, const char * fileName, const int tableSize)  {
    (*lookupTable).tableSize = tableSize;
	unsigned long long NR_TOP = 1 << (tableSize * 2);
	(*lookupTable).table = malloc(sizeof(unsigned) * NR_TOP);
	int fin = open(fileName, O_RDONLY);
	unsigned step = 1048576;
	unsigned int i;
	for (i = 0; i < NR_TOP; i += step) {
		read(fin, (*lookupTable).table + i, step * sizeof(*(*lookupTable).table));
	}
	close(fin);
}

unsigned int LookupSafe(LOOKUPTABLE lookupTable, BWT * bwt,
                        unsigned long long lKey, unsigned long long rKey,
                        unsigned int *l, unsigned int *r) {

	*l = lKey ? lookupTable.table[lKey-1]+1 : 1;
	*r = lookupTable.table[rKey];

	if (*l == bwt->inverseSa0) {
		(*l)++;
	}

	return *r-*l+1;
}

unsigned int retrieveSA=0,retrieveHASH=0;
double textPositionTime, textPositionTimeTotal = 0;
unsigned int writeQIndex;

double getTextPositionTime() {return textPositionTimeTotal;}
unsigned int getSARetrieved() {return retrieveSA;}
unsigned int getHASHRetrieved() {return retrieveHASH;}

void FreeLookupTable(LOOKUPTABLE * lookupTable) {
     free((*lookupTable).table);
}

void LoadHashTable(HASHTABLE * hashTable, const char * fileName) {
unsigned int ttlOccurrence=0;
unsigned int ttlItem=0;

    FILE *inFile;
    if(!(inFile = fopen(fileName, "r"))) return;
    fread((unsigned int *)&((*hashTable).tableSize),sizeof(unsigned int),1,inFile);
    fread((unsigned int *)&((*hashTable).a),sizeof(unsigned int),1,inFile);
    fread((unsigned int *)&((*hashTable).b),sizeof(unsigned int),1,inFile);
    fread((unsigned int *)&((*hashTable).prime),sizeof(unsigned int),1,inFile);
    fread((unsigned int *)&ttlItem,sizeof(unsigned int),1,inFile);
    fread((unsigned int *)&ttlOccurrence,sizeof(unsigned int),1,inFile);

    //printf("Initializing the hash table..(n=%u)\n",(*hashTable).tableSize);
    (*hashTable).table = (HASHCELL*) malloc(sizeof(HASHCELL)*((*hashTable).tableSize));
    (*hashTable).itemList = (HASHITEM*) malloc(sizeof(HASHITEM)*ttlItem);
    (*hashTable).occList = (OCC*) malloc(sizeof(OCC)*ttlOccurrence);
    //printf("Initialized the hash table..\n");
    unsigned int i;
    for (i=0;i<((*hashTable).tableSize);i++) {
        char mk;
        fread((char *) &mk,1,1,inFile);
        if (mk==0) {
            //Empty cell
            (*hashTable).table[i].index=0;
            (*hashTable).table[i].count=0;
        } else {
            fread((unsigned int *)&((*hashTable).table[i].index),sizeof(unsigned int),1,inFile);
            fread((unsigned int *)&((*hashTable).table[i].count),sizeof(unsigned int),1,inFile);
        }
    }
    for (i=0;i<ttlItem;i++) {
        fread((unsigned int *)&((*hashTable).itemList[i].l),sizeof(unsigned int),1,inFile);
        fread((unsigned int *)&((*hashTable).itemList[i].r),sizeof(unsigned int),1,inFile);
        fread((unsigned int *)&((*hashTable).itemList[i].occIndex),sizeof(unsigned int),1,inFile);
    }
    for (i=0;i<ttlOccurrence;i++) {
        fread((unsigned int *)&((*hashTable).occList[i]),sizeof(unsigned int),1,inFile);
    }
    fclose(inFile);

}

void FreeHashTable(HASHTABLE * hashTable) {
     free((*hashTable).table);
     free((*hashTable).itemList);
     free((*hashTable).occList);
}

unsigned int Hash(HASHTABLE * hashTable, unsigned int key) {
    //g(x)=(ax+b) mod p
    unsigned long long multipleA=(key * (*hashTable).a)% (*hashTable).prime;
    unsigned int g=(unsigned int) (( multipleA + (*hashTable).b ) % (*hashTable).prime);

    //f(x)=g(x) mod 2n
    unsigned int f=g % ((*hashTable).tableSize);

    return f;
}

HASHITEM * HashFind(HASHTABLE * hashTable, unsigned int l,unsigned int r) {
    unsigned int hashedIndex=Hash(hashTable,l);
    unsigned int index=(*hashTable).table[hashedIndex].index;
    unsigned int count=(*hashTable).table[hashedIndex].count;

    if (hashedIndex>=0 && hashedIndex<(*hashTable).tableSize) {
        unsigned int i=0;
        while (i<count && ((*hashTable).itemList[index+i].l!=l || (*hashTable).itemList[index+i].r!=r)) {
            i++;
        }
        if (i>=count) return NULL;
        return &((*hashTable).itemList[index+i]);
    }
    return NULL;
}

void RegisterDecoder(BWT * bwt,HASHTABLE * hashTable) {
     occBwt=bwt;
     occHashtable=hashTable;
     occCollected=0;
     retrieveSA=0;
     retrieveHASH=0;
     textPositionTimeTotal = 0;
     writeQIndex = 0;
     //occCollector = malloc(sizeof (unsigned int) * 1024*1024);
}

unsigned int allOne = 0;
unsigned int OCCSection=0;

inline int CalMismatch(const char *seq, const unsigned int *ref, const unsigned int occPosCord, const unsigned int seqLen, const unsigned int dnaLength){
	unsigned int i, l;
	int match = 0;
//	fprintf(stderr, "%u\t%u\t%u\n", occPosCord, seqLen, dnaLength);
//	fprintf(stderr, "%u\n", ref[0]);
	for(i =0, l=occPosCord; i < seqLen && l < dnaLength; ++i, ++l){
//		fprintf(stderr, "%d,%u ", i, l);
		if(!(((*(seq+i))&0x3) ^ ((((*(ref+(l>>4)))>>(((~l)&0xf)<<1)))&0x3)))
			++match;
	}
//	fprintf(stderr, "match %d\n", match);
	return (seqLen-match);
}

#include <assert.h>

int OCCProcess(const unsigned int l,const unsigned int r, const BWTOPT *bo, const unsigned int info, HITTABLE *hits) {
#if TRUE
//	fprintf(stderr, "OCC Process, n occ %d\n", r-l+1);
#endif
	const unsigned int cutoff = bo->cutoff;
	if(hits->n >= bo->cutoff) return 0;
	unsigned int n = hits->n;
	HITITEM *hit = hits->itemList+n;
	ChrBlock *blockList   = bo->blockList;
	const unsigned int nblock   = bo->nblock;
	const unsigned int seqLen   = bo->seqLen;
	const unsigned int alnLen  = bo->alnLen;
	const unsigned int extLen = seqLen-alnLen;
	const unsigned int max_mm   = bo->max_mm+(info>>25 & 0x7);
	const unsigned int *pacRef  = bo->pacRef;
	const unsigned int dnaLength = bo->dnaLen;
	const unsigned int strain = (info>>24)&1;
	char *seq    =  strain?bo->rc:bo->fw;
//	for(i=0; i<extLen; ++i) fprintf(stderr,"%c","ACGT"[*(seq+i)]);
//	fprintf(stderr, "%d\n", extLen);
	unsigned int occ_pos  = 0;
	if (r-l+1 < 4) {
			//SA
		unsigned int k;
		for (k=l;k<=r && n < cutoff;k++) {
			occ_pos = BWTSaValue(occBwt,k);
			HitInc(n);
		}
		int inc = n - hits->n;
		hits->n = n;
		return inc;
	} else {
			//Hash
		HASHITEM *item = HashFind(occHashtable,l,r);
		if (item != NULL) {
			unsigned int k;
			for (k=0;k<item->r-item->l+1 && n < cutoff;k++) {
				occ_pos = occHashtable->occList[item->occIndex+k];
				HitInc(n);
			}
			int inc = n - hits->n;
			hits->n = n;
			return inc;
		} else {
			unsigned int k;
			for (k=l; k<=r && n < cutoff; k++) {
				occ_pos = BWTSaValue(occBwt,k);
				HitInc(n);
			}
			int inc = n - hits->n;
			hits->n = n;
			return inc;
		}
	}
}

void registerTPFile(FILE * filePtr,unsigned int searchMode) {
    textPositionFile=filePtr;
    fwrite(&searchMode,sizeof(unsigned int),1,textPositionFile);
	allOne=(1U<<31)-1;
	allOne<<=1;
	allOne+=1;
}

void registerQIndex(unsigned int index) {
    writeQIndex=index;
    OCCSection=0;
}
void registerQSection() {
    if (writeQIndex==0) {
        fwrite(&allOne,sizeof(unsigned int),1,textPositionFile);
    } else {
        OCCSection++;
    }
}
