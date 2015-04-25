/*
 * =============================================================================
 *
 *       Filename:  soapio.c
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
#include "soapio.h"
#include <assert.h>

#define BUF_SIZE 0xF00000
#define BUF_PER_LINE 0x800
       
void FreeMultiSeq(MULTISEQ *mseqs){
	int i; 
	int n = mseqs->n;
	ALNSEQ *p ;
	for(i=0; i<n; ++i){
		p = mseqs->seqList+i;
		free(p->name); free(p->seq); free(p->rc); free(p->qual); free(p->rcqual);
		if(p->report>0){
			free(p->itemList->md);
			free(p->itemList->cigar);
			free(p->itemList);
		}
	}
}

static inline  char * reverse(const  char *seq, int len){
	char *rc = ( char *)malloc(sizeof( char) * (len+1));
	int i;
	for(i=len;i>0;--i) *(rc+len-i) = *(seq+i-1);
	return rc;
}

#define SEQDUP(dest, ori, pe) {         \
	dest->flag = pe;		\
	dest->tid = -1;			\
	dest->id  = id;			\
	dest->len = ori.l;		\
	dest->nhits = dest->report = 0;		\
	dest->ns = ori.ns;		\
	dest->seq = (char *)malloc(sizeof(char)*ori.l);		\
	memcpy(dest->seq, ori.seq, ori.l);		\
	dest->rc = reverse(ori.rc, ori.l);		\
	dest->qual = strdup(ori.qual);		\
	dest->rcqual = reverse(ori.qual, ori.l);		\
	dest->rcqual[ori.l] = '\0';			\
	dest->name = strdup(ori.name);		\
	dest++;			\
}

int GenMultiReads(const HSP *hsp, MULTISEQ *mseqs, const int len, const int pe, unsigned int *start, int *nb){
	const unsigned int *pacRef = hsp->packedDNA;
	char **chrName = hsp->chrName;
	const unsigned int refLen = hsp->dnaLength;
	const ChrBlock  *blockList = hsp->blockList;
	ALNSEQ *alnSeq;
	alnSeq = mseqs->seqList;
	int num, i, j;
	num = i = j = 0;
	unsigned int st = *start;
	int n = *nb ;
	unsigned int ori = (blockList+n)->ori;
	unsigned int blockStart = (blockList+n)->blockStart;
	unsigned int blockEnd = (blockList+n)->blockEnd;
	while (num < MAX_MULTI_READS && st < refLen-len) {
		alnSeq->name = (char *)malloc(sizeof(char)*MAX_NAME_LEN);
		alnSeq->seq = (char *)malloc(sizeof(char)*(len+1));
		alnSeq->rc = (char *)malloc(sizeof(char)*(len+1));
		alnSeq->qual = (char *)malloc(sizeof(char)*(len+1));
		for (i=0; i<len; ++i) {
			*(alnSeq->seq+i) = ((*(pacRef+((st+i)>>4)))>>(30-((st+i)&0xf)*2))&0x3;
			*(alnSeq->rc+len-i-1) = (~((*(alnSeq->seq+i))&0x3))&0x3;
			*(alnSeq->qual+i) = 'h';
		}
		if (st+len-1>blockEnd) {
			n++;
			ori = (blockList+n)->ori;
			blockStart = (blockList+n)->blockStart;
			blockEnd = (blockList+n)->blockEnd;
			continue;
		}
		sprintf(alnSeq->name, ">%s_%d", chrName[(blockList+n)->chrID], st-blockStart+ori+1);
		alnSeq->name[strlen(alnSeq->name)] = '\0';
		alnSeq->report = 0;
		alnSeq->nhits = 0;
		alnSeq->qual[i] = '\0';
		alnSeq->ns = 0;
		alnSeq->tid = -1;
		alnSeq->flag = 0;
		alnSeq->len = len;
		alnSeq->id = st++;
		num++;
		alnSeq++;
	}
	*start = st;
	*nb = n;
	mseqs->n = num;
	return num;
}

int GetMultiSeq (InFileList *ifp, MULTISEQ *mseqs, const int pe, int(*get_read)(FILE * , seq_t * , const int)){
#ifdef DEBUG
//	fprintf(stderr, "Get Multi Seqs\n");
#endif
	ALNSEQ *alnSeq;
	alnSeq= mseqs->seqList;
	int num, len, id;
	FILE * ifpA, * ifpB;
	num = 0;
	ifpA = ifp->ifpA; ifpB = ifp->ifpB;
	id  = ifp->id;
	seq_t tmp;
	tmp.max = tmp.l = 0; 
	tmp.seq = tmp.rc = tmp.qual = NULL;
	while(num < MAX_MULTI_READS){
		tmp.ns = 0;
		if (feof(ifpA)||(pe && feof(ifpB))) break;
		if ((len=get_read(ifpA, &tmp, TRUE))>0){
			SEQDUP(alnSeq, tmp, pe);
			++num;
			if(pe &&(len=get_read(ifpB, &tmp, TRUE))>0){
				SEQDUP(alnSeq, tmp, pe);
				++num ;
			}
			++id;
		}
	}
	mseqs->n = num;
	ifp->id = id;
//	fprintf(stderr, "%d\n", num);
	/*
#ifdef DEBUG
	int j;
	fprintf(stderr, "len :%d\n", len);
	fprintf(stderr, "fw\n");
	for(j = 0; j<len; j++) fprintf(stderr, "%d", *(tmp.seq+j));
	fprintf(stderr, "\n");
	fprintf(stderr, "rc\n");
	for(j = 0; j<len; j++) fprintf(stderr, "%d", *(tmp.rc+j));
	fprintf(stderr, "\n");
	alnSeq--;
	for(j = 0; j < alnSeq->len; j++) fprintf(stderr, "rc %d", *(alnSeq->rc+j));
	fprintf(stderr, "\n");
	fprintf(stderr, "%s\n%s\n", alnSeq->name, alnSeq->qual);
	fprintf(stderr, "soap get multisequences ...\n");
#endif
	exit(0);
	//*/
	free(tmp.seq); free(tmp.qual);free(tmp.rc);
	return num;
}

#define SOAPOUT(file){	\
/*	fprintf (stderr, "generate format\n");		\
	*/			\
	int k = 0;		\
	if(o->id)		\
		ksprintf(str, "%d\t", alnSeq->id);\
	else 		\
		ksprintf(str, "%s\t", alnSeq->name);	\
	int n_cigar = hit->n_cigar;			\
	int beg=0, end=len;			\
	if(hit->cigar[0]>>14 == 3) beg = hit->cigar[0]&0x3ff;		\
	if(hit->cigar[n_cigar-1]>>14 == 3) end = len - (hit->cigar[n_cigar-1]&0x3ff);		\
	if(strain){		\
		for(k=beg; k<end; ++k){		\
			kputc("ACGTN"[(int)*(rc+k)], str);		\
		}				\
		ksprintf(str, "\t");			\
		for(k=beg; k<end; ++k){		\
			kputc(alnSeq->rcqual[k], str);		\
		}				\
		ksprintf(str, "\t");			\
	} else {		\
		for(k=beg; k<end; ++k){		\
			kputc("ACGTN"[(int)*(seq+k)], str);		\
		}				\
		ksprintf(str, "\t");			\
		for(k=beg; k<end; ++k){		\
			kputc(alnSeq->qual[k], str);		\
		}				\
		ksprintf(str, "\t");			\
	}			\
	/*			\
	fprintf(stderr, "%s\n", alnSeq->qual);			\
	fprintf(stderr, "%d\n", alnSeq->nhits);			\
	fprintf(stderr, "%c\n", "ab"[file]); fprintf(stderr, "%d\n", alnSeq->len);fprintf(stderr,"%c\n", "+-"[strain]);fprintf(stderr, "%s\n", chrName[(alnSeq->itemList+j)->chr]);fprintf(stderr, "%u\n", alnSeq->itemList->pos);			\
	*/		\
	ksprintf(str, "%d\t%c\t%d\t%c\t%s\t%d\t", alnSeq->nhits,"ab"[file], end-beg, "+-"[strain], chrName[hit->chr], hit->pos);			\
	if(!n_seedMM)ksprintf(str, "0\t");	\
	else if (n_seedMM == 1) 			\
		ksprintf(str, "1\t%c->%d%c%d\t", "ACGT"[(info_seedMM>>8)&3], info_seedMM&0xff, strain?"ACGT"[(int)rc[info_seedMM&0xff]]:"ACGT"[(int)seq[info_seedMM&0xff]], (strain?alnSeq->qual[len-(info_seedMM&0xff)]:alnSeq->qual[info_seedMM&0xff])-'@');		\
	else if (n_seedMM == 2) 		\
		ksprintf(str, "2\t%c->%d%c%d\t%c->%d%c%d\t", 	\
				"ACGT"[(info_seedMM>>8)&3], info_seedMM&0xff, strain?"ACGT"[(int)rc[info_seedMM&0xff]]:"ACGT"[(int)seq[info_seedMM&0xff]], (strain?alnSeq->qual[len-(info_seedMM&0xff)-1]:alnSeq->qual[info_seedMM&0xff])-'@',				\
				"ACGT"[(info_seedMM>>20)&3], (info_seedMM>>12)&0xff, strain?"ACGT"[(int)rc[(info_seedMM>>12)&0xff]]:"ACGT"[(int)seq[(info_seedMM>>12)&0xff]], (strain?alnSeq->qual[len-1-((info_seedMM>>12)&0xff)]:alnSeq->qual[(info_seedMM>>12)&0xff])-'@');		\
	else if (n_seedMM == 3) {		\
		ksprintf(str, "%d\t%d\t", (100+1+hit->n_gape), (info_seedMM>>12)&0xff);		\
	}		\
	else if (n_seedMM == 4) {		\
		ksprintf(str, "%d\t%d\t", (200+1+hit->n_gape), (info_seedMM>>12)&0xff);		\
	}		\
	if ((alnSeq->itemList+j)->n_cigar){		\
		for (k=0; k<hit->n_cigar;k++)			\
			ksprintf(str, "%d%c", hit->cigar[k]&0x3ff, "MIDS"[(hit->cigar[k]>>14)]);		\
	}else ksprintf(str, "%dM\t", alnSeq->len);		\
	ksprintf(str, "\t%s\n", hit->md);			\
}

#if 0
#define BINARY_SOAP() {			\
	fwrite(&(alnSeq->itemList+j)->id, sizeof(unsigned int), 1, ofp);				\
	fwrite(&(alnSeq->itemList+j)->chr, sizeof(unsigned int), 1, ofp);				\
	fwrite(&(alnSeq->itemList+j)->pos, sizeof(unsigned int), 1, ofp);				\
	fwrite(&(alnSeq->itemList+j)->len, sizeof(unsigned int), 1, ofp);			\
	int k;				\
	if (strain) {			\
		for(k=0; k<(alnSeq->itemList+j)->len; ++k)		\
			fwrite();						\
	} else {		\
		for(k=0; k<(alnSeq->itemList+j)->len; ++k)		\
			fwrite();						\
	}				\
}
#endif

#include <assert.h>

void DumpAln(MULTISEQ *mseqs, OUTAUX *o, OutFileList *ofp,unsigned int *nAln, unsigned int *nSE){
	int n, i;
	n = mseqs->n;
	unsigned int n_aln = *nAln;
	unsigned int n_se  = *nSE;
	char **chrName = o->chrName ;
	ALNSEQ *alnSeq;
	HITITEM *hit;
	FILE *ofpAln, *ofpSe, *ofpUn;
	ofpAln = ofp->ofpAln;
	ofpSe  = ofp->ofpSe;
	ofpUn  = ofp->ofpUn;
	kstring_t *str = (kstring_t *)calloc(1, sizeof(kstring_t));
	for (i=0; i<n; ++i) {
		alnSeq = mseqs->seqList + i;
		int j = 0;
		if (alnSeq->report) {
			unsigned int strain, n_seedMM, n_mm, info_seedMM, len, flag;
			char *seq = alnSeq->seq;
			char *rc = alnSeq->rc;
			flag = alnSeq->flag;
			len = alnSeq->len;
			if((alnSeq->flag>>1&0x1) || !(alnSeq->flag&0x1)) {
				hit = alnSeq->itemList;
				for(j=0; j<alnSeq->report; ++j) {
					str->l = 0;
					strain = hit->strain;
					n_seedMM = hit->info >> 25 & 0x7;
					info_seedMM = hit->info & 0xffffff;
					n_mm = hit->n_mm;
//					int file =(alnSeq->flag)&1?(i&1):0;
					SOAPOUT((alnSeq->flag&1)?(i&1):0);
					fprintf(ofpAln, "%s", str->s);
					++hit;
				}
				++n_aln;
			} else {
				hit = alnSeq->itemList;
				for(j=0; j<alnSeq->report; ++j) {
					str->l = 0;
					strain = hit->strain;
					n_seedMM = hit->info >> 25 & 0x7;
					info_seedMM = hit->info & 0xffffff;
					n_mm = hit->n_mm;
					SOAPOUT((alnSeq->flag&1)?(i&1):0);
					fprintf(ofpSe, "%s", str->s);
					++hit;
				}
				++n_se;
			}
		} else if (o->un) {
			fprintf(ofpUn, ">%s\n", alnSeq->name);
			int j=0;
			for(;j<alnSeq->len;j++)
				fprintf(ofpUn, "%c", "ACGT"[(int)*(alnSeq->seq+j)]);
			fprintf(ofpUn, "\n");
		}
	}
	*nAln = n_aln;
	*nSE  = n_se;
	free(str->s);
	free(str);
}
