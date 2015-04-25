/*
 * =============================================================================
 *
 *       Filename:  soapio.h
 *
 *    Description:  
 *
 *       Revision:  none
 *       Compiler:  gcc 4.3.2 or above
 *
 *         Author:  Chang Yu (yc), yuchang@genomics.org.cn
 *        Company:  BGI Shenzhen
 *      CopyRight:  Copyright (c) 2009, BGI Shenzhen
 *
 * =============================================================================
 */
#ifndef  _SOAPIO_H_
#define  _SOAPIO_H_
#include <stdio.h>
#include <stdlib.h>
#include "extratools.h"
#include "SeqIO.h"
#define MAX_MULTI_READS 0x20000

typedef struct _ALNSEQ_TYPE_{
	int tid;
	int id, len, ns;
	char *name, *seq, *rc, *qual, *rcqual;
	unsigned int flag;
	int nhits;
	struct {
		int H0;
		int H1;
		int H2;
	}top;
	int report;
	HITITEM *itemList;
}ALNSEQ;

typedef struct _MULTISEQ_TYPE_{
	int n, max;
	ALNSEQ *seqList;
}MULTISEQ;

typedef struct _OUTAUX_TYPE_{
	int id, un, chrNum;
	char **chrName;
	int allErr;
}OUTAUX;

typedef struct _INFILELIST_{
	FILE *ifpA, *ifpB;
        int id;
	int lock;
}InFileList;

typedef struct _OUTFILELIST_{
	FILE *ofpAln, *ofpSe, *ofpUn;
        int id;
	int lock;
}OutFileList;

void FreeMultiSeq(MULTISEQ *);
int GenMultiReads(const HSP *, MULTISEQ *, const int , const int , unsigned int *, int *);
int GetMultiSeq (InFileList *, MULTISEQ *, const int , int (*)(FILE *, seq_t *, const int));
void DumpAln(MULTISEQ *, OUTAUX *, OutFileList *,unsigned int *, unsigned int *);

#endif   /* ----- #ifndef SOAPIO_INC  ----- */
