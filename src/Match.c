/*
 * =============================================================================
 *
 *       Filename:  Match.c
 *
 *    Description:  *
 *       Revision:  none
 *       Compiler:  gcc 4.3.2 or above
 *
 *         Author:  Chang Yu (yc), yuchang@genomics.org.cn
 *        Company:  BGI Shenzhen
 *      CopyRight:  Copyright (c) 2009, BGI Shenzhen
 *
 * =============================================================================
 */

#include "Match.h"

#define GenCigarMD() {		\
	if (hits->itemList[i].n_cigar == 0){						\
		int j, match;				\
		match = 0;					\
		for(j=0; j<len && occPos < dnaLength; ++j, ++occPos){			\
			unsigned c = (((*(pacRef+(occPos>>4)))>>(((~occPos)&0xf)<<1)))&0x3;				\
			if(!(((*(seq+j))&0x3) ^ c)){				\
				++match;					\
			} else {					\
				if(match||!j) ksprintf(str, "%d", match);				\
				kputc("ACGTN"[c], str);					\
				match = 0;						\
			}								\
		}				\
		ksprintf(str, "%d", match);				\
		(alnSeq->itemList+i)->n_cigar = 1;						\
		(alnSeq->itemList+i)->cigar = (unsigned short *)malloc(sizeof(unsigned short)*(alnSeq->itemList->n_cigar));	\
		(alnSeq->itemList+i)->cigar[0] = (FROM_M << 14) | (len & 0x3ff);		\
		(alnSeq->itemList+i)->md = strdup(str->s);				\
	} else {				\
		int n_cigar = hits->itemList[i].n_cigar;			\
		hits->itemList[i].n_cigar = 0;			\
		(alnSeq->itemList+i)->cigar = (unsigned short *)malloc(sizeof(unsigned short)*(1+n_cigar));	\
		memcpy((alnSeq->itemList+i)->cigar, cigar, n_cigar*sizeof(unsigned short));	\
		unsigned int x = (alnSeq->itemList+i)->occ_pos;				\
		unsigned int y, z;							\
		y = z = 0;						\
		int k, l, u;				\
		unsigned char c;						\
		for (k = u = 0; k < n_cigar; ++k) {					\
			l = cigar[k]&0x3fff;					\
			if (cigar[k]>>14 == FROM_M) {					\
				for (z = 0; z < l && x+z < dnaLength; ++z) {		\
					c = (((*(pacRef+((x+z)>>4)))>>(((~(x+z))&0xf)<<1))) & 0x3;	\
					if (c > 3 || seq[y+z] > 3 || c != seq[y+z]) {	\
						if(u||!(y+z)) ksprintf(str, "%d", u);			\
						kputc("ACGTN"[c], str);			\
						u = 0;					\
					} else ++u;					\
				}																					\
				x += l; y += l;																		\
			} else if (cigar[k]>>14 == FROM_I || cigar[k]>>14 == 3) {								\
				y += l;																				\
			} else if (cigar[k]>>14 == FROM_D) {													\
				ksprintf(str, "%d", u);																\
				kputc('D', str);																	\
				for (z = 0; z < l && x+z < dnaLength; ++z)											\
					kputc("ACGTN"[(((*(pacRef+((x+z)>>4)))>>(((~(x+z))&0xf)<<1))) & 0x3], str);		\
				u = 0;			\
				x += l;			\
			}			\
		}			\
/*		free(cigar); cigar = NULL;	*/		\
		ksprintf(str, "%d", u);\
		(alnSeq->itemList+i)->md = strdup(str->s);			\
	}			\
}

inline void PickupHit(ALNSEQ *alnSeq, const int rr,int *site, HITTABLE *hits, const unsigned int *pacRef, const unsigned int dnaLength, unsigned short *cigar){
#ifdef DEBUG
//		fprintf(stderr, "Pick up for output\n");
#endif
	int i = *site;
	kstring_t *str = (kstring_t *)calloc(1, sizeof(kstring_t));
	str->l = 0; str->m = 0;
	if (!hits->n || (hits->n > 1 && !rr)) {alnSeq->report = 0; alnSeq->nhits =0; return;}
	else {
		int n;
		n = hits->n;
		if(rr == 1 || rr == 0) {
			alnSeq->report = 1;
			alnSeq->itemList = (HITITEM *)malloc(sizeof(HITITEM) *1);
//			assert(i<hits->n);
			HITCPY(alnSeq->itemList, hits->itemList+i);
			unsigned int occPos = alnSeq->itemList->occ_pos;
			unsigned int len = alnSeq->len;	
			char *seq = alnSeq->itemList->strain?alnSeq->rc:alnSeq->seq;	
			if (hits->itemList[i].n_cigar == 0){						
				int j, match;	
				match = 0;			
				for(j=0; j<len && occPos < dnaLength; ++j, ++occPos){	
					unsigned c = (((*(pacRef+(occPos>>4)))>>(((~occPos)&0xf)<<1)))&0x3;
					if(!(((*(seq+j))&0x3) ^ c)){
						++match;	
					} else {
						if(match || !j)ksprintf(str, "%d", match);
						kputc("ACGTN"[c], str);
						match = 0;
					}				
				}
				if(match)ksprintf(str, "%d", match);
				alnSeq->itemList->n_cigar = 1;
				alnSeq->itemList->cigar = (unsigned short *)malloc(sizeof(unsigned short)*(alnSeq->itemList->n_cigar));		
				alnSeq->itemList->cigar[0] = (FROM_M << 14) | (len & 0x3fff);
				alnSeq->itemList->md = strdup(str->s);
//		fprintf(stderr, "%d%c\n", alnSeq->itemList->cigar[0]&0x3ff, "MIDS"[alnSeq->itemList->cigar[0]>>14]);
			} else {
				int n_cigar = hits->itemList[i].n_cigar;
				hits->itemList[i].n_cigar = 0;
				alnSeq->itemList->cigar = (unsigned short *)malloc(sizeof(unsigned short)*(1+n_cigar));
				memcpy(alnSeq->itemList->cigar, cigar, n_cigar*sizeof(unsigned short));
				unsigned int x = alnSeq->itemList->occ_pos;
				unsigned int y, z;
				y = z = 0;
				int k, l, u;	
				k = l = u = 0;
				unsigned char c;
				for (k = u = 0; k < n_cigar; ++k) {
					l = cigar[k]&0x3fff;
					if (cigar[k]>>14 == FROM_M) {
						for (z = 0; z < l && x+z < dnaLength; ++z) {
							c = (((*(pacRef+((x+z)>>4)))>>(((~(x+z))&0xf)<<1))) & 0x3;
							if (c > 3 || seq[y+z] > 3 || c != seq[y+z]) {
								if(u||!(y+z))ksprintf(str, "%d", u);
								kputc("ACGTN"[c], str);
								u = 0;
							} else ++u;
						}
						x += l; y += l;
					} else if (cigar[k]>>14 == FROM_I || cigar[k]>>14 == 3) {
						y += l;
					} else if (cigar[k]>>14 == FROM_D) {
						ksprintf(str, "%d", u);
						kputc('D', str);
						for (z = 0; z < l && x+z < dnaLength; ++z)
							kputc("ACGTN"[(((*(pacRef+((x+z)>>4)))>>(((~(x+z))&0xf)<<1))) & 0x3], str);		
						u = 0;
						x += l;
					}
				}
				if (u) ksprintf(str, "%d", u);
				alnSeq->itemList->md = strdup(str->s);
			}
//			GenCigarMD();
			alnSeq->nhits = n;
		} else {
			alnSeq->report = n;
			alnSeq->itemList = (HITITEM *)malloc(sizeof(HITITEM) * n);
			for (i = 0; i < n; ++i){
				str->l = 0;
				HITCPY(alnSeq->itemList+i, hits->itemList+i);
				unsigned int occPos = (alnSeq->itemList+i)->occ_pos;
				unsigned int len = alnSeq->len;	
				char *seq = (alnSeq->itemList+i)->strain?alnSeq->rc:alnSeq->seq;
				GenCigarMD();
			}
			alnSeq->nhits = n;
		}
	}
	free(str->s);
	free(str);
}

void SEAlnCore(int tid, MULTISEQ *mseqs, BWT *bwt, BWT *rev_bwt, LOOKUPTABLE *lookup, LOOKUPTABLE *rev_lookup, HSP *hsp, const SOAPOPT *opt) {
	int i;
	ALNSEQ *alnSeq;
	HITTABLE *hits;
	hits = (HITTABLE *) malloc (sizeof(HITTABLE) * 1);
	hits->itemList = (HITITEM *) malloc (sizeof(HITITEM) * MAX_ALN);
	/*
	   for(i=0; i<MAX_ALN; ++p,++i)
	   p->path = (unsigned short *) malloc (sizeof(unsigned short) * MAX_DIFF);
	//*/
	int mode, seedLen, ns, rr, cutoff;
	mode = opt->mode; rr = opt->rr;ns = opt->ns; seedLen = opt->aln_len; cutoff = opt->cutoff;

	if (opt->uniq)  mode = 4;
	BWTOPT bo;
	//*
	bo.nblock = hsp->numOfBlock; bo.blockList = hsp->blockList;
	bo.cutoff = MAX_ALN;
	bo.max_mm = opt->max_mm; bo.gap_len = opt->gap_len; bo.gap_fb = opt->gap_fb;
	bo.pacRef = hsp->packedDNA; bo.dnaLen = hsp->dnaLength;
//	int count = 0;
	//*/
#ifdef DEBUG
	//	fprintf(stderr, "%d\n", mseqs->n);
#endif
	for(i=0; i < mseqs->n; i+=1){
		//		fprintf(stderr, "n reads %d\n", i);
		alnSeq = mseqs->seqList+i;
		//*

#ifdef  PTHREADS
		ALNSEQ *p = mseqs->seqList + i;
		if (opt->nthreads > 1) {		
			pthread_mutex_lock(&lock);	
			if (alnSeq->tid < 0) { 	
				int j;	
				for (j = i; j < mseqs->n && j < i + NSEQ_PER_THREAD; ++j)
					p++->tid = tid;
			} else if (alnSeq->tid != tid) {
				pthread_mutex_unlock(&lock);
				continue;
			}
			pthread_mutex_unlock(&lock);	
		}
#endif    ///* -----  not PTHREADS  -----*/

		if (alnSeq->ns <= ns){
			int h0, h1, h2, h3;
			h0 = h1 = h2 = h3 = 0;
			hits->n = 0;
			bo.seqLen = bo.alnLen = alnSeq->len;
			unsigned int extLen = 0;
			bo.fw = alnSeq->seq;
			bo.rc = alnSeq->rc;
ALIGN:
			bo.h = bo.alnLen>>1;
			bo.x = bo.alnLen>39?bo.alnLen/3:13;
			bo.y = bo.alnLen-13-bo.x;
			if (bo.y <= 0) {fprintf(stderr, "length y < 0, countinue as 13\n");}
			switch (mode) {
				case 5:
				case 4:
				case 0:
					h0  = BWTExactMatching((unsigned char *)alnSeq->seq, &bo, FORWARD, bwt, lookup, hits);
					h0 += BWTExactMatching((unsigned char *)alnSeq->rc+extLen,  &bo, REVERSE, bwt, lookup, hits);
					if(hits->n >= cutoff || mode == 0) break;
				case 1:
					h1  = BWT1ErrorMatching((unsigned char *)alnSeq->seq, &bo, FORWARD, bwt, rev_bwt, lookup, rev_lookup, hits);
					h1 += BWT1ErrorMatching((unsigned char *)alnSeq->rc+extLen,  &bo, REVERSE, bwt, rev_bwt, lookup, rev_lookup, hits);
					if(hits->n >= cutoff || mode == 1) break;
				case 2:
					h2  = BWT2ErrorMatching((unsigned char *)alnSeq->seq, &bo, FORWARD, bwt, rev_bwt, lookup, rev_lookup, hits);
					h2 += BWT2ErrorMatching((unsigned char *)alnSeq->rc+extLen,  &bo, REVERSE, bwt, rev_bwt, lookup, rev_lookup, hits);
					if(mode == 4 || hits->n >= cutoff || mode == 2) break;
			}

			if (!hits->n && seedLen < bo.alnLen) {
				bo.alnLen = seedLen;
				extLen = alnSeq->len-seedLen;
				goto ALIGN;
			}

			if (hits->n) {
				alnSeq->flag = 0;
				int site = hits->n?(hits->n == 1?0:rand()%hits->n):-1;
				PickupHit(alnSeq, rr, &site, hits, hsp->packedDNA, hsp->dnaLength, NULL);
			}else{
				alnSeq->flag = 0;
				alnSeq->report = 0;
			}

		}else{
			alnSeq->flag = 0;
			alnSeq->report = 0;
		}
		
		//*/
	}
//	fprintf(stderr, "Alignment Time: %2.7f\n", getElapsedTime(startTime));
	free(hits->itemList);free(hits);
}

#ifdef PTHREADS

typedef struct _THREADAUX_TYPE_{
	int tid;
	BWT *bwt;
	BWT *rev_bwt;
	LOOKUPTABLE *lookup;
	LOOKUPTABLE *rev_lookup;
	HSP *hsp;
	MULTISEQ *mseqs;
	SOAPOPT *o;
}THREADAUX;

static void *Workers(void *threadAux){
	THREADAUX *aux = (THREADAUX *)threadAux;
	aux->o->pe? PEAlnCore(aux->tid, aux->mseqs, aux->bwt, aux->rev_bwt, aux->lookup,aux->rev_lookup, aux->hsp, aux->o)
		:SEAlnCore(aux->tid, aux->mseqs, aux->bwt, aux->rev_bwt, aux->lookup, aux->rev_lookup, aux->hsp, aux->o);
}
#endif


void MatchProcess (FILEDS *fds, BWT *bwt, BWT *rev_bwt, LOOKUPTABLE *lookup, LOOKUPTABLE *rev_lookup, HSP *hsp, SOAPOPT * const opt) {

	InFileList *ifp;
	OutFileList *ofp;
	ifp = (InFileList *) malloc(sizeof(InFileList) * 1);
	ofp = (OutFileList *) malloc(sizeof(OutFileList) * 1);
	ifp->ifpA   = fdopen(fds->ifdA, "r");
	ofp->ofpAln = fdopen(fds->ofdAln, "w");
	if (opt->pe) {
		ifp->ifpB   = fdopen(fds->ifdB, "r");
		ofp->ofpSe  = fdopen(fds->ofdSe, "w");
	}
	if (opt->unmapped)
		ofp->ofpUn  = fdopen(fds->ofdUn, "w");
	int fast = CheckFast(fds->ifdA);
	ifp->id = 0;
	
	MULTISEQ mseqs;
	mseqs.n = mseqs.max = 0;
	mseqs.seqList = (ALNSEQ *)malloc(sizeof(ALNSEQ) * MAX_MULTI_READS);

#define INITALN(aln) {					\
	int j;						\
	for(j=0;j<MAX_MULTI_READS;j++){				\
		aln[j].tid = aln[j].id = aln[j].len = aln[j].ns = aln[j].flag = aln[j].nhits = aln[j].report = 0;			\
	}		\
}
	INITALN(mseqs.seqList);
	unsigned int i, nseq, nAln, nSE;
	i = nseq = nAln = nSE = 0;
//	double startAlnTime = setStartTime();
	double startAlnTime; 
	
	OUTAUX o_aux;

//*
	o_aux.allErr = opt->allErr;
	o_aux.un = opt->unmapped;
	o_aux.id = opt->id;
	o_aux.chrName = hsp->chrName;
	o_aux.chrNum = hsp->chrNum;
	//*/

#ifdef DEBUG
//	fprintf(stderr, "Begin Aln process\n");
#endif
	while (GetMultiSeq(ifp, &mseqs,opt->pe, fast?fastq:fasta) != 0) {
		nseq += mseqs.n;
		startAlnTime = setStartTime();

#ifndef PTHREADS
//		fprintf(stderr, "no threads\n");
		opt->pe ? PEAlnCore(0, &mseqs, bwt, rev_bwt, lookup, rev_lookup, hsp, opt)
				:SEAlnCore(0, &mseqs, bwt, rev_bwt, lookup, rev_lookup, hsp, opt);
#else 
		if(opt->nthreads <= 1)
			opt->pe?PEAlnCore(0, &mseqs, bwt, rev_bwt, lookup, rev_lookup, hsp, opt)
				:SEAlnCore(0, &mseqs, bwt, rev_bwt, lookup, rev_lookup, hsp, opt);
		else {
			pthread_t *tid;
			pthread_attr_t attr;
			THREADAUX *threadAux;
			int j;
			pthread_attr_init(&attr);
			pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
			threadAux = (THREADAUX *)calloc(opt->nthreads, sizeof(THREADAUX));
			tid = (pthread_t*)calloc(opt->nthreads, sizeof(pthread_t));
			for (j = 0; j < opt->nthreads; ++j) {
				threadAux[j].tid = j;
				threadAux[j].bwt = bwt; threadAux[j].rev_bwt = rev_bwt;
				threadAux[j].lookup = lookup; threadAux[j].rev_lookup = rev_lookup;
				threadAux[j].hsp = hsp;
				threadAux[j].mseqs = &mseqs; 
				threadAux[j].o = opt;
				pthread_create(&tid[j], &attr, Workers, threadAux + j);
			}
			pthread_attr_destroy(&attr);
			for (j = 0; j < opt->nthreads; ++j) pthread_join(tid[j], 0);
			free(threadAux); free(tid);
		}
#endif

		fprintf(stderr, "%d ok %7.2f sec\n", nseq, getElapsedTime(startAlnTime));
		DumpAln(&mseqs, &o_aux, ofp, &nAln, &nSE);
		FreeMultiSeq(&mseqs);
	}
	if (opt->pe) 
		fprintf(stderr, "Total Pairs: %d PE\n"
				"Paired:      %d (%5.2f%%) PE\n"
				"Singled:     %d (%5.2f%%) SE\n", nseq/2, nAln/2, (float)nAln/nseq*100, nSE, (float)nSE/(nseq)*100);
	else 
		fprintf(stderr, "Total Reads: %d\n"
				"Alignment:   %d (%5.2f%%)\n", nseq, nAln, (float)nAln/nseq*100);
	free(mseqs.seqList);
	fclose(ifp->ifpA); fclose(ofp->ofpAln);
	if(opt->pe){fclose(ifp->ifpB); fclose(ofp->ofpSe);}
	if(opt->unmapped) fclose(ofp->ofpUn);
	free(ifp); free(ofp);
}                               /*                              */
