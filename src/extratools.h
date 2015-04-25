#ifndef _EXTRATOOLS_H_
#define _EXTRATOOLS_H_

#include <stdio.h>
#include <stdlib.h>
#include "MiscUtilities.h" 
#include "MemManager.h"
#include "TextConverter.h"
#include "Timing.h"
#include "BWT.h"
#include "kstring.h"
#include <fcntl.h>
#include <unistd.h>
#define MAX_DIFF 32

typedef struct LOOKUPTABLE_TYPE {
       unsigned int tableSize;
       unsigned int * table;
}LOOKUPTABLE;

typedef struct HASHCELL_TYPE {
    unsigned int count;
    unsigned int index;
}HASHCELL;

typedef struct HASHITEM_TYPE {
    unsigned int l;
    unsigned int r;
    unsigned int occIndex;
}HASHITEM;

typedef unsigned int OCC;

typedef struct HASHTABLE_TYPE {
       unsigned int prime;
       unsigned int a;
       unsigned int b;
       unsigned int tableSize;
       HASHCELL * table;
       HASHITEM * itemList;
       OCC * occList;
}HASHTABLE;

typedef struct _HITITEM_TYPE_{
	int info;
	int strain;
	int chr;
	unsigned int occ_pos, pos, blockid;
	int n_diff;
	int n_mm;
	int n_gapo, n_gape, gap_beg;
	int n_cigar;
	char *md;
	unsigned short *cigar;
}HITITEM;

typedef struct _HITTABLE_TYPE_{
	int n;
	HITITEM *itemList;
}HITTABLE;

BWT * occBwt;
HASHTABLE * occHashtable;
unsigned int * occCollector;
unsigned int occCollected;

FILE * textPositionFile;
void registerTPFile(FILE * filePtr,unsigned int searchMode);

void registerQIndex(unsigned int queryIndex);
void registerQSection();

void LoadLookupTable(LOOKUPTABLE * lookupTable, const char * fileName, const int tableSize);
void FreeLookupTable(LOOKUPTABLE * lookupTable);
unsigned int LookupSafe(LOOKUPTABLE lookupTable, BWT * bwt,unsigned long long lKey, unsigned long long rKey,unsigned int *l, unsigned int *r);
void LoadHashTable(HASHTABLE * hashTable, const char * fileName);
HASHITEM * HashFind(HASHTABLE * hashTable, unsigned int l,unsigned int r);
void FreeHashTable(HASHTABLE * hashTable);
void RegisterDecoder(BWT * bwt,HASHTABLE * hashTable);

//void OCCClean();
//void OCCProcess(unsigned int l,unsigned int r);
inline int altCalMM(unsigned int x);
inline int CalMismatch(const char *,const unsigned int *,const unsigned int , const unsigned int, const unsigned int);
int OCCProcess(unsigned int l, unsigned int r, const BWTOPT *bo, const unsigned int info, HITTABLE *hits);

#define GenOCCArr(arr) do{ \
	int occ = 0;			\
	if (r-l+1 >= 4) {	\
		HASHITEM *item = HashFind(occHashtable,(l),(r));	\
		if (item==NULL) {	\
			unsigned int k;	\
			for (k=l;k<=r;++k) {	\
				arr[occ++] = BWTSaValue(occBwt,k);	\
			}\
		} else {\
			unsigned int k;\
			for (k=0;k<item->r-item->l+1;++k) {\
				arr[occ++] = occHashtable->occList[item->occIndex+k];\
			}\
		}\
	} else {\
		unsigned int k;\
		for (k=l;k<=r;++k) {\
			arr[occ++] = BWTSaValue(occBwt,k);\
		}\
	}\
}while(0);

#define OrientPacPos(){\
	int start, end;		\
	start = end = 0;			\
	int l, m, h;		\
	l = 0; h = nblock; m = nblock/2;	\
	/*			\
	fprintf(stderr, "pacPos %u\n",occ_pos);			\
	*/			\
	while(l<=h){		\
		m = (h+l)>>1;	\
		if((start=blockList[m].blockStart)>occ_pos){ 	\
			h = m - 1;				\
		}		\
		else if((end = blockList[m].blockEnd)<occ_pos){	\
			l = m + 1;		\
		}			\
		else if ((start <= (occ_pos-(strain?extLen:0))) && (end >= (occ_pos + (strain?alnLen:seqLen)))){	\
			chr = blockList[m].chrID; 			\
			pos = occ_pos - start + (blockList+m)->ori + 1;			\
			blockid = m;			\
			break;						\
		}else break;					\
	}				\
}

#define MAX_MD_LEN 1024
#define HitInc(n) {	\
	/*		\
	if(cutoff == n)		\
		fprintf(stderr, "max %d->n %d\n", cutoff, n);		\
		*/		\
	int chr = -1;int pos = -1;int blockid = 0;int mm = 0;		\
	OrientPacPos();		\
/*	if(chr > -1 && (pos-(strain?extLen:0)) > 0 && (!extLen || max_mm >= ((info>>24)&0x7>3?0:((info>>24)&0x7))+(mm=CalMismatch(seq, pacRef, strain?(occ_pos-extLen):(occ_pos+alnLen), extLen, errTmp, strain?0:alnLen, allErr)))) {		*/ \
	if(chr > -1 && (pos-(strain?extLen:0)) > 0 && (max_mm >= (mm=CalMismatch(seq, pacRef, strain?(occ_pos-extLen):(occ_pos), seqLen, dnaLength)))) {		\
		hit->strain = strain;			\
		hit->chr = chr;		\
		hit->pos = pos-(strain?extLen:0);		\
		if (hit->pos < 0 )  printf("%d\t%d\n", pos, extLen);				\
		hit->blockid = blockid;			\
		hit->occ_pos = occ_pos-(strain?extLen:0);	\
		hit->info = info;		\
		hit->n_cigar = 0;			\
		hit->n_mm = mm + ((info>>24)&0x7);			\
		n++; hit++;			\
	}				\
}

//void CleanDecoder();

double getTextPositionTime();
unsigned int getSARetrieved();
unsigned int getHASHRetrieved();

#endif /*_EXTRATOOLS_H_*/

