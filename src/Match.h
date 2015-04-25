/*
 * =============================================================================
 *
 *       Filename:  Match.h
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


#ifndef  _MATCH_H_
#define  _MATCH_H_

#include "SeqIO.h"
#include "BWTAln.h"
#include "BWT.h"
#include "extratools.h"
#include "soapio.h"
#include "stdaln.h"

#ifdef PTHREADS
#include <pthread.h>
#define NSEQ_PER_THREAD 0xF00
static pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
#define SEQ_ALLOC() {\
	if (opt->nthreads > 1) {		\
		pthread_mutex_lock(&lock);	\
		if (alnSeq->tid < 0) { 		\
			int j;			\
			for (j = i; j < mseqs->n && j < i + NSEQ_PER_THREAD; ++j)	\
				alnSeq[j].tid = tid;		\
		} else if (alnSeq->tid != tid) {			\
			pthread_mutex_unlock(&lock);	\
			continue;				\
		}						\
		pthread_mutex_unlock(&lock);			\
	}		\
}
#else 
#define SEQ_ALLOC() 
#endif

#define MULTI_SEQ 0x100000
#ifndef MAX_MISMATCH
#define MAX_MISMATCH 20
#endif

#ifndef MAX_GAP_LEN
#define MAX_GAP_LEN 10
#endif
#define MAX_SEQ_LEN 256
#define MAX_ALN 10000
#define FORWARD	0
#define REVERSE	1
#define ALN_MAT 0
#define ALN_MIS 0x11
#define ALN_INS 0x22
#define ALN_DEL 0x33 
#include <assert.h>

#define HITCPY(dest, ori) {	\
	(dest)->info    = (ori)->info;		\
	(dest)->strain  = (ori)->strain;		\
	(dest)->chr     = (ori)->chr;		\
	(dest)->occ_pos = (ori)->occ_pos;		\
	(dest)->pos     = (ori)->pos;		\
	(dest)->n_mm    = (ori)->n_mm;		\
	(dest)->n_gapo  = (ori)->n_gapo;		\
	(dest)->n_gape  = (ori)->n_gape;		\
	(dest)->gap_beg  = (ori)->gap_beg;		\
	(dest)->n_diff  = (ori)->n_diff;		\
	(dest)->n_cigar = (ori)->n_cigar;			\
}

#define PacReadExt(fw, rc, start, len, seqPac, rcPac)  {\
	int j;			\
	for(j=0; j<len; ++j){				\
		seqPac[j>>4] <<= 2;			\
		seqPac[j>>4] |= *(fw+j+start);			\
		rcPac[j>>4] <<= 2;			\
		rcPac[j>>4] |= *(rc+j);			\
	}									\
}

typedef struct _SOAPOPT_{
	int fast, o_format, chain;
	int aln_len, ns, max_mm, gap_len, gap_fb;
	int mode, cutoff; 
	int pe;
	int zero_qual;
	int min_ins, max_ins, FR;
	int rr;
	int unmapped;
	int nthreads;	//number of pthreads
	int id;
	int bisulfite;
	int allErr;
	int min_len;
	int uniq;
}SOAPOPT;

typedef struct _FILEDS_{
	int ifdA, ifdB;
	int ofdAln, ofdSe, ofdUn;
}FILEDS;

typedef struct _MATCHAUX_TYPE_{
	int max_mm;
	int len, ext;
	unsigned int *pac;
	unsigned int dnaLen;
	int allErr;
}MATCHAUX;

typedef struct _PEAUX_TYPE_{
	int min_ins, max_ins;
	int FR;
        int cutoff, len;
	int allErr;
}PEAUX;

inline int CheckIns(HITITEM *, HITITEM *, PEAUX *);
void MatchProcess (FILEDS *, BWT *, BWT *, LOOKUPTABLE *, LOOKUPTABLE *, HSP *, SOAPOPT * const );
inline void PickupHit(ALNSEQ *, const int ,int *, HITTABLE *,const unsigned int *, const unsigned int, unsigned short * );
void SEAlnCore(int , MULTISEQ *, BWT *, BWT *, LOOKUPTABLE *, LOOKUPTABLE *, HSP *, const SOAPOPT *);
void PEAlnCore(int , MULTISEQ *, BWT *, BWT *, LOOKUPTABLE *, LOOKUPTABLE *, HSP *, const SOAPOPT *);
int HITCMP(const void *a, const void *b);

#endif   /* ----- #ifndef _MATCH_H_INC  ----- */

