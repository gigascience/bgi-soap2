/*
 * =============================================================================
 *
 *       Filename:  PariMatch.c
 *
 *    Description:  
 *
 *       Revision:  none
 *       Compiler:  gcc 4.3.2 or aboAve
 *
 *         Author:  Chang Yu (yc), yuchang@genomics.org.cn
 *        Company:  BGI Shenzhen
 *      CopyRight:  Copyright (c) 2009, BGI Shenzhen
 *
 * =============================================================================
 */


#include "Match.h"

inline int CheckIns(HITITEM *p, HITITEM *q, PEAUX *o) {
	int strain1 = (p->info >> 24)&1;
	int strain2 = (q->info >> 24)&1;
	if(p->chr != q->chr || strain1 == strain2) return FALSE;
	else if(o->FR) {
		if(!strain1 && q->pos-p->pos+o->len >= o->min_ins && q->pos-p->pos+o->len <= o->max_ins) return TRUE;
		else if(strain1 && p->pos-q->pos+o->len >= o->min_ins && p->pos-q->pos+o->len <= o->max_ins) return TRUE;
		else{ 
			return FALSE;
		};
	}
	else if(!o->FR){
		if(strain1 && q->pos-p->pos >= o->min_ins && q->pos-p->pos+o->len <= o->max_ins) return TRUE;
		else if(!strain1 && p->pos-q->pos >= o->min_ins && p->pos-q->pos+o->len <= o->max_ins) return TRUE;
		else{
			return FALSE;
		}
	}
	return TRUE;
}

int HITCMP(const void *a, const void *b){
       if ((*(HITITEM *)a).chr != (*(HITITEM *)b).chr){
		return (*(HITITEM *)a).chr - (*(HITITEM *)b).chr;
	} else 
		return (*(HITITEM *)a).pos - (*(HITITEM *)b).pos;
}

int GenPair(HITTABLE **hitse, PEAUX *po,  HITTABLE **hitpe) {
	if(!hitse[0]->n || !hitse[1]->n) return 0;
	HITITEM *p, *q;
	p = hitse[0]->itemList; q = hitse[1]->itemList;
	const int cutoff = po->cutoff;
	if(hitse[0]->n == 1 && hitse[1]->n == 1 ){
		if(CheckIns(p, q, po)){
			HITCPY(hitpe[0]->itemList+hitpe[0]->n, p);
			HITCPY(hitpe[1]->itemList+hitpe[1]->n, q);
			hitpe[0]->n++; hitpe[1]->n++;
			return 1;
		}else{
			return  0;
		}
	}else{
		if(hitse[0]->n > 1)qsort(hitse[0]->itemList, hitse[0]->n, sizeof(HITITEM), HITCMP);
		if(hitse[1]->n > 1)qsort(hitse[1]->itemList, hitse[1]->n, sizeof(HITITEM), HITCMP);
		int n = hitpe[0]->n;
//		fprintf(stderr, "%d\n", n);
		p=hitse[0]->itemList; q = hitse[1]->itemList;
		while (p!=hitse[0]->itemList+hitse[0]->n) {
			while(p!=hitse[0]->itemList+hitse[0]->n && p->chr<q->chr)p++;
			while(q!=hitse[1]->itemList+hitse[1]->n && p->chr>q->chr)q++;
			if (p==(hitse[0]->itemList+hitse[0]->n) || q==(hitse[1]->itemList+hitse[1]->n)) return n ;
			while (p->chr==q->chr && q!=(hitse[1]->itemList+hitse[1]->n)){
				if (CheckIns(p, q, po)) {
					HITCPY(hitpe[0]->itemList+n, p);
					HITCPY(hitpe[1]->itemList+n, q);
					n++; hitpe[0]->n++; hitpe[1]->n++;
					if(n >= cutoff){
						hitpe[0]->n = hitpe[1]->n = n;
						return n;
					}
				}
				q++;
			}
			p++;
			q = hitse[1]->itemList;
		}
		hitpe[0]->n = hitpe[1]->n = n;
		return n;
	}
}

#if 1	
unsigned short *SWRescue(const ALNSEQ *alnSeq, const BWTOPT *bo, const PEAUX *po, const int rescue, HITTABLE **hitse, HITTABLE **hitpe, int *nc, int *n_rescue){

	HITITEM *hitf = hitse[rescue^1]->itemList;
	int nfound = hitse[rescue^1]->n;

	const unsigned int *pacRef = bo->pacRef;
	const unsigned int dnaLen = bo->dnaLen;
	char *seq;
	ChrBlock *blockList = bo->blockList;

	int minIns, maxIns, len, keyLength, n;
	unsigned int occPos, beg;

	keyLength = alnSeq->len;
	minIns = po->min_ins; maxIns = po->max_ins;
	len = maxIns-minIns+3*keyLength;
	occPos = beg = 0;
	AlnParam ap = aln_param_bwa;

	path_t *path, *p;
	int i, path_len, n_cigar;
	path_len = n_cigar = 0;

	cigar_t * cigar = NULL;
	path = (path_t *)calloc((len+keyLength), sizeof(path_t));


	unsigned char *refSeq = (unsigned char *)calloc(len, sizeof(unsigned char));

	int SWCutoffX, SWCutoffY;
	SWCutoffX = SWCutoffY = bo->min_len < keyLength ? bo->min_len : (keyLength < 17 ? keyLength : 17);

	HITITEM *peItem1, *peItem2;
	peItem1 = hitpe[rescue^1]->itemList;
	peItem2 = hitpe[rescue]->itemList;

	int mm = 10;

	int n_mm, n_gapo, n_gape, gap_beg, ed_dist;
	n_mm = n_gapo = n_gape = gap_beg = 0;
	ed_dist = keyLength;

	unsigned short tmp_cigar[16];

//	fprintf(stderr, "%d\n", nfound);
	for (i = 0; i < nfound; ++i) {

		occPos = (hitf+i)->occ_pos;

		if((((hitf+i)->info>>24)&0x7) > mm) continue;
		mm = (hitf+i)->info>>24&0x7;
		n = (hitf+i)->blockid;
		if(po->FR ^ (hitf+i)->strain) {
			beg = occPos + minIns - keyLength;
			if(beg + len >= (blockList + n)->blockEnd) continue;
			seq = (hitf+i)->strain ? alnSeq->seq : alnSeq->rc;
		} else {
			beg = occPos - maxIns - keyLength;
			if(beg < (blockList + n)->blockStart) continue;
			seq = (hitf + i)->strain ? alnSeq->seq : alnSeq->rc;
		}

		{
			unsigned char *p = refSeq;
			unsigned int j, l;
			for(j=beg, l=0; l<len && j<dnaLen; ++l, ++j)
				*p++ = (((*(pacRef+(j>>4)))>>(((~j)&0xf)<<1))&0x3);
		}

//		fprintf(stderr, "%d	n_cigar %d\n", i, n_cigar);
		if (n_cigar) {free(cigar); n_cigar = 0;}
//		fprintf(stderr, "%d	n_cigar %d\n", i, n_cigar);

		aln_local_core(refSeq, len, (unsigned char *)seq, keyLength, &ap, path, &path_len, 1);
		cigar = aln_path2cigar(path, path_len, &n_cigar);
		int k, x, y;
		x = y = k = 0;
		for (k = 0, x = y = 0; k < n_cigar; ++k) {
			unsigned short c = cigar[k];
			if (c>>14 == FROM_M) x += c&0x3fff, y += c&0x3fff;
			else if (c>>14 == FROM_D) x += c&0x3fff;
			else y += c&0x3fff;
		}

		if (x < SWCutoffX  && y < SWCutoffY) continue;

		{ // update cigar and coordinate;
			SWCutoffX = x; SWCutoffY = y;
			int start, end;
			p = path + path_len - 1;
			beg += (p->i? p->i : 1) - 1;
			start = (p->j? p->j : 1) - 1;
			end = path->j;
			cigar = (unsigned short*)realloc(cigar, 2 * (n_cigar + 2));
			if (start) {
				memmove(cigar + 1, cigar, 2 * (n_cigar));
				cigar[0] = 3<<14 | start;
				++(n_cigar);
			}
			if (end < keyLength) {
				cigar[n_cigar] = 3<<14 | (keyLength - end);
				++(n_cigar);
			}
		}

		n_mm = n_gapo = n_gape = gap_beg = 0;
		int indel = 3;
		{
			p = path + path_len - 1;
			x = p->i? p->i - 1 : 0; y = p->j? p->j - 1 : 0;
			int l=0;
			for (k = 0; k < n_cigar; ++k) {
				unsigned short c = cigar[k];
				if (c>>14 == FROM_M) {
					for (l = 0; l < (c&0x3fff); ++l)
						if (refSeq[x+l] < 4 && seq[y+l] < 4 && refSeq[x+l] != seq[y+l]) ++n_mm;
					x += c&0x3fff, y += c&0x3fff;
				} else if (c>>14 == FROM_D) {
					indel = 3; gap_beg = y; x += c&0x3fff; ++n_gapo; n_gape += (c&0x3fff) - 1;
				} else if (c>>14 == FROM_I){
					indel = 4; gap_beg = y; y += c&0x3fff; ++n_gapo; n_gape += (c&0x3fff) - 1;
				}
			}
			if(n_mm >= bo->max_mm || n_gapo > 1 || n_gape + n_gapo > bo->gap_len) continue;
			if (!n_gapo) indel=0;
		}

		*n_rescue += 1;

		if (n_gape + n_gapo + n_mm < ed_dist)
		{// update pe hit
			hitpe[rescue^1]->n = hitpe[rescue]->n = 1;
			HITCPY(peItem1, hitf+i);
			peItem2->chr = peItem1->chr;
			peItem2->pos = beg-occPos+(hitf+i)->pos;
			peItem2->occ_pos = beg;
			peItem2->strain = 1 ^ peItem1->strain;
			peItem2->n_mm = n_mm;
			peItem2->n_gapo = n_gapo;
			peItem2->n_gape = n_gape;
			peItem2->info = 0;
			peItem2->info |= ((indel<<25) | ((gap_beg&0xff)<<12) | ((n_gape+1)&0xff));
			peItem2->gap_beg = gap_beg;
			peItem2->n_cigar = n_cigar;
			ed_dist = n_gape + n_gapo + n_mm;
			for(k=0; k < n_cigar; ++k)tmp_cigar[k] = cigar[k];
		}

	}

	if (*n_rescue) {
		if (n_cigar < peItem2->n_cigar) cigar = (unsigned short *) calloc (peItem2->n_cigar, sizeof(unsigned short));
		for(i = 0; i < peItem2->n_cigar; ++i) cigar[i] = tmp_cigar[i];
		*nc = peItem2->n_cigar;
	}

	free(path);
	free(refSeq);
	return cigar;

}

#endif

void PEAlnCore(int tid, MULTISEQ *mseqs, BWT *bwt, BWT *rev_bwt, LOOKUPTABLE *lookup, LOOKUPTABLE *rev_lookup, HSP *hsp,const SOAPOPT *opt) {
	int i;
	ALNSEQ *alnSeq[2];
 	HITTABLE *hitse[2], *hitpe[2];
	hitse[0] = (HITTABLE *)malloc(sizeof(HITTABLE));
	hitse[1] = (HITTABLE *)malloc(sizeof(HITTABLE));
	hitpe[0] = (HITTABLE *)malloc(sizeof(HITTABLE));
	hitpe[1] = (HITTABLE *)malloc(sizeof(HITTABLE));
	hitse[0]->itemList = (HITITEM *) malloc (sizeof(HITITEM) * MAX_ALN);
	hitse[1]->itemList = (HITITEM *) malloc (sizeof(HITITEM) * MAX_ALN);
	hitpe[0]->itemList = (HITITEM *) malloc (sizeof(HITITEM) * (MAX_ALN+1));
	hitpe[1]->itemList = (HITITEM *) malloc (sizeof(HITITEM) * (MAX_ALN+1));
	const int multiTotal = mseqs->n;
	PEAUX pe_aux;
	BWTOPT boA, boB;
	int mode, cutoff, ns, seedLen, rr;
	mode = opt->mode; cutoff = opt->cutoff; ns = opt->ns; seedLen = opt->aln_len; rr = opt->rr;
	boB.nblock = boA.nblock = hsp->numOfBlock;boB.blockList = boA.blockList = hsp->blockList;
	boB.cutoff=boA.cutoff = MAX_ALN; boB.gap_len = boA.gap_len = opt->gap_len; boB.gap_fb = boA.gap_fb = opt->gap_fb;
	boB.max_mm = boA.max_mm = opt->max_mm; boB.pacRef = boA.pacRef = hsp->packedDNA;
	boB.dnaLen = boA.dnaLen = hsp->dnaLength;
	boA.min_len = boB.min_len = opt->min_len;
	boA.h = boA.x = boA.y = boB.h = boB.x = boB.y = 0;
        pe_aux.min_ins = opt->min_ins; pe_aux.max_ins = opt->max_ins; 
	pe_aux.FR = opt->FR; pe_aux.len = 0; pe_aux.cutoff = MAX_ALN;
	int x = 0;
	int se, pe, non;
	se=pe=non=0;
	double swBeg, swTime;
	swBeg = swTime = 0;
//	int swRun=0;
	int nRescue = 0;
	for(i=0; i < multiTotal; i += 2){
//		fprintf(stderr, "%d\n", i);
#ifdef  PTHREADS
		if (opt->nthreads > 1) {
			pthread_mutex_lock(&lock);
			ALNSEQ *p = mseqs->seqList+i;
			if (p->tid < 0) {
				int j;
				int pend = multiTotal-i;
				for (j = 0; j < pend && j < NSEQ_PER_THREAD; j+=2){
					(p+j)->tid = (p+j+1)->tid = tid;
				}
			} else if (p->tid != tid) {
				pthread_mutex_unlock(&lock);
				continue;
			}
			pthread_mutex_unlock(&lock);
		}
#endif

		alnSeq[0] = mseqs->seqList+i;
		alnSeq[1] = mseqs->seqList+i+1;
		hitse[0]->n = hitse[1]->n = hitpe[0]->n = hitpe[1]->n = 0;

		if(alnSeq[0]->ns <= ns || alnSeq[1]->ns <= ns){
			int nc = 0;
			nRescue = 0;
			x+=2;
			int ah0, ah1, ah2, bh0, bh1, bh2, ah3, bh3;
			boA.seqLen = boA.alnLen = alnSeq[0]->len;
			boB.seqLen = boB.alnLen = alnSeq[1]->len;
			unsigned int extLen = 0;
			pe_aux.len = alnSeq[0]->len;
			boA.fw = alnSeq[0]->seq; boA.rc = alnSeq[0]->rc;
			boB.fw = alnSeq[1]->seq; boB.rc = alnSeq[1]->rc;

ALIGN:
			ah0 = ah1 = ah2 = ah3 = bh0 = bh1 = bh2 = bh3 = 0;
			boA.h = boA.alnLen>>1;
			boA.x = boA.y = boA.alnLen>=39?boA.alnLen/3:(boA.alnLen>=32 && boA.alnLen<39)?10:7;
			boB.h = boB.alnLen>>1;
			boB.x = boB.y = boB.alnLen>=39?boB.alnLen/3:(boB.alnLen>=32 && boB.alnLen<39)?10:7;
			switch (mode) {
				case 5:
				case 4: cutoff = opt->cutoff;
				case 0:
					ah0  = BWTExactMatching((unsigned char*)alnSeq[0]->seq, &boA, FORWARD, bwt, lookup, hitse[0]);
					ah0 += BWTExactMatching((unsigned char*)alnSeq[0]->rc+extLen, &boA, REVERSE, bwt, lookup, hitse[0]);
					bh0  = BWTExactMatching((unsigned char*)alnSeq[1]->seq, &boB, FORWARD, bwt, lookup, hitse[1]);
					bh0 += BWTExactMatching((unsigned char*)alnSeq[1]->rc+extLen, &boB, REVERSE, bwt, lookup, hitse[1]);
					if (ah0 && bh0) {
							 GenPair(hitse, &pe_aux, hitpe);}
					if (hitpe[0]->n >= cutoff || mode == 0) break;
				case 1:
					ah1  = BWT1ErrorMatching((unsigned char*)alnSeq[0]->seq, &boA,  FORWARD, bwt, rev_bwt, lookup, rev_lookup, hitse[0]);
					ah1 += BWT1ErrorMatching((unsigned char*)alnSeq[0]->rc+extLen, &boA, REVERSE, bwt, rev_bwt, lookup, rev_lookup, hitse[0]);
					bh1  = BWT1ErrorMatching((unsigned char*)alnSeq[1]->seq, &boB, FORWARD, bwt, rev_bwt, lookup, rev_lookup, hitse[1]);
					bh1 += BWT1ErrorMatching((unsigned char*)alnSeq[1]->rc+extLen, &boB, REVERSE, bwt, rev_bwt, lookup, rev_lookup, hitse[1]);
					if (ah1 || bh1) {
							 GenPair(hitse, &pe_aux, hitpe);}
					if (hitpe[0]->n >= cutoff || mode == 1) break;
				case 2:
					ah2  = BWT2ErrorMatching((unsigned char*)alnSeq[0]->seq, &boA, FORWARD, bwt, rev_bwt,  lookup, rev_lookup, hitse[0]);
					ah2 += BWT2ErrorMatching((unsigned char*)alnSeq[0]->rc+extLen, &boA, REVERSE, bwt, rev_bwt,  lookup, rev_lookup, hitse[0]);
					bh2  = BWT2ErrorMatching((unsigned char*)alnSeq[1]->seq, &boB, FORWARD, bwt, rev_bwt, lookup, rev_lookup, hitse[1]);
					bh2 += BWT2ErrorMatching((unsigned char*)alnSeq[1]->rc+extLen, &boB, REVERSE, bwt, rev_bwt,  lookup, rev_lookup, hitse[1]);
					if (ah2 || bh2){ 
							 GenPair(hitse, &pe_aux, hitpe);}
					if (hitpe[0]->n >= cutoff || mode == 4 || mode == 2) break;
			}

			if (seedLen<boB.alnLen) {
				if (!hitse[0]->n && !hitse[1]->n && (seedLen<boA.alnLen && seedLen< boB.alnLen)) {
					boB.alnLen = boA.alnLen = seedLen;
					boA.extLen = alnSeq[0]->len-seedLen;
					boB.extLen = alnSeq[1]->len-seedLen;
					if (alnSeq[0]->len < seedLen || alnSeq[1]->len <seedLen){
						fprintf(stderr, "read_len shorter than seed_len%d. Continue\n", seedLen);
						goto OUTPUT;
					}
					goto ALIGN;
				} else if (!hitpe[0]->n &&  !hitse[0]->n && seedLen<boA.alnLen) {
					if (alnSeq[0]->len < seedLen){
						fprintf(stderr, "read_len shorter than seed_len%d. Continue\n", seedLen);
						goto OUTPUT;
					}
					boA.alnLen = seedLen; boA.extLen = alnSeq[0]->len - seedLen; boA.h = boA.alnLen>>1;
					boA.x = boA.y = boA.alnLen>=39?boA.alnLen/3:(boA.alnLen>=32 && boA.alnLen<39)?10:7;
					ah0  = BWTExactMatching((unsigned char*)alnSeq[0]->seq, &boA, FORWARD, bwt, lookup, hitse[0]);
					ah0 += BWTExactMatching((unsigned char*)alnSeq[0]->rc+boA.extLen, &boA, REVERSE, bwt, lookup, hitse[0]);
					if (ah0 && GenPair(hitse, &pe_aux, hitpe)) goto OUTPUT;
					ah1  = BWT1ErrorMatching((unsigned char*)alnSeq[0]->seq, &boA,  FORWARD, bwt, rev_bwt, lookup, rev_lookup, hitse[0]);
					ah1 += BWT1ErrorMatching((unsigned char*)alnSeq[0]->rc+boA.extLen, &boA, REVERSE, bwt, rev_bwt, lookup, rev_lookup, hitse[0]);
					if (ah1 && GenPair(hitse, &pe_aux,  hitpe))goto OUTPUT;
					ah2  = BWT2ErrorMatching((unsigned char*)alnSeq[0]->seq, &boA, FORWARD, bwt, rev_bwt, lookup, rev_lookup, hitse[0]);
					ah2 += BWT2ErrorMatching((unsigned char*)alnSeq[0]->rc+boA.extLen, &boA, REVERSE, bwt, rev_bwt, lookup, rev_lookup, hitse[0]);
					if (ah2 && GenPair(hitse, &pe_aux, hitpe))goto OUTPUT;
				} else if (!hitpe[1]->n &&  !hitse[1]->n &&  seedLen <boB.alnLen) {
					if (alnSeq[1]->len < seedLen){
						fprintf(stderr, "read_len shorter than seed_len%d. Continue\n", seedLen);
						goto OUTPUT;
					}
					boB.alnLen = seedLen; boB.extLen = alnSeq[1]->len - seedLen; boB.h = boB.alnLen>>1;
					boB.x = boB.y = boB.alnLen>=39?boB.alnLen/3:(boB.alnLen>=32 && boB.alnLen<39)?10:7;
					bh0  = BWTExactMatching((unsigned char*)alnSeq[1]->seq, &boB, FORWARD, bwt, lookup, hitse[1]);
					bh0 += BWTExactMatching((unsigned char*)alnSeq[1]->rc+boB.extLen, &boB, REVERSE, bwt, lookup, hitse[1]);
					if(bh0 && GenPair(hitse, &pe_aux,  hitpe)) goto OUTPUT;
					bh1  = BWT1ErrorMatching((unsigned char*)alnSeq[1]->seq, &boB, FORWARD, bwt, rev_bwt, lookup, rev_lookup, hitse[1]);
					bh1 += BWT1ErrorMatching((unsigned char*)alnSeq[1]->rc+boB.extLen, &boB, REVERSE, bwt, rev_bwt, lookup, rev_lookup, hitse[1]);
					if(bh1 && GenPair(hitse, &pe_aux, hitpe)) goto OUTPUT;
					bh2  = BWT2ErrorMatching((unsigned char*)alnSeq[1]->seq, &boB, FORWARD, bwt, rev_bwt, lookup, rev_lookup, hitse[1]);
					bh2 += BWT2ErrorMatching((unsigned char*)alnSeq[1]->rc+boB.extLen, &boB, REVERSE, bwt, rev_bwt, lookup, rev_lookup, hitse[1]);
					if(bh2 && GenPair(hitse, &pe_aux,  hitpe)) goto OUTPUT;
				}
			}
			
			unsigned short * cigar = NULL;
//			if (hitse[1]->n && !hitse[0]->n && boA.gap_len){                       /* gap goto sw */
			if (hitse[1]->n && !hitse[0]->n && (boA.min_len < alnSeq[0]->len || boA.gap_len)){
//				swBeg = setStartTime();
				cigar = SWRescue(alnSeq[0], &boA, &pe_aux, 0, hitse, hitpe, &nc, &nRescue); 
//				swTime += getElapsedTime(swBeg);
//				swRun++;
				goto OUTPUT;
//			} else if (!hitse[1]->n && hitse[0]->n && boA.gap_len) {            /* gap goto sw */
			} else if (!hitse[1]->n && hitse[0]->n && (boB.min_len < alnSeq[1]->len || boA.gap_len)) {
//				swBeg = setStartTime();
				cigar = SWRescue(alnSeq[1], &boB, &pe_aux, 1, hitse, hitpe, &nc, &nRescue);
				
//				swTime += getElapsedTime(swBeg);
//				swRun++;
				goto OUTPUT;
			}

OUTPUT:
			if (hitpe[0]->n && hitpe[1]->n){
				pe+=2;
//				assert(hitpe[0]->n==hitpe[1]->n);
//				printf("site: %d\n", hitpe[0]->n);
				int site = (hitpe[0]->n == 1)?0:rand()%hitpe[0]->n;
				alnSeq[0]->flag = alnSeq[1]->flag = 0x3;
//				printf("site: %d\n", site);
//				printf("out:%d\t%d\n", hitpe[0]->itemList[site].n_cigar, nc);
				PickupHit(alnSeq[0], rr, &site, hitpe[0], hsp->packedDNA, hsp->dnaLength, cigar);
//				printf("out:%d\n", hitpe[1]->itemList[site].n_cigar);
				PickupHit(alnSeq[1], rr, &site, hitpe[1], hsp->packedDNA, hsp->dnaLength, cigar);
				if (nRescue) alnSeq[0]->nhits = alnSeq[1]->nhits = nRescue;
			} else {
				int site = 0;
				if (hitse[0]->n && hitse[1]->n) {
					se+=2;
					site = hitse[0]->n == 1?0:rand()%hitse[0]->n;
					PickupHit(alnSeq[0], rr, &site, hitse[0], hsp->packedDNA, hsp->dnaLength, cigar);
					site = hitse[1]->n == 1?0:rand()%hitse[1]->n;
					PickupHit(alnSeq[1], rr, &site, hitse[1], hsp->packedDNA, hsp->dnaLength, cigar);
					alnSeq[0]->flag = alnSeq[1]->flag = 1;
				} else if(!hitse[1]->n && hitse[0]->n) {
					se++;
					site = hitse[0]->n == 1?0:rand()%hitse[0]->n;
					PickupHit(alnSeq[0], rr, &site, hitse[0], hsp->packedDNA, hsp->dnaLength, cigar);
					alnSeq[1]->flag |= 0x8;
					non++;
				} else if(!hitse[0]->n && hitse[1]->n) {
					se++;
					site = hitse[1]->n == 1?0:rand()%hitse[1]->n;
					PickupHit(alnSeq[1], rr, &site, hitse[1], hsp->packedDNA, hsp->dnaLength, cigar);
					alnSeq[0]->flag |= 0x8;
					non++;
				} else {
					non+=2;
					alnSeq[0]->flag |= 0x12;
					alnSeq[1]->flag |= 0x12;
					alnSeq[0]->report = 0;
					alnSeq[1]->report = 0;
				}
			}
			if(nc) {
				free(cigar);
				nc = 0;
			}
		} else {
			non+=2;
			alnSeq[0]->flag |= 0x12;
			alnSeq[1]->flag |= 0x12;
			alnSeq[0]->report = 0;
			alnSeq[1]->report = 0;
		}
	}
	free(hitse[0]->itemList);
	free(hitse[1]->itemList);
	free(hitpe[0]->itemList);
	free(hitpe[1]->itemList);
	free(hitse[0]);free(hitse[1]);free(hitpe[0]);free(hitpe[1]);
}
