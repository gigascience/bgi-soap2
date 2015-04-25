/*
 * =============================================================================
 *
 *       Filename:  SeqIO.h
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

#ifndef  __SEQIO_H__
#define  __SEQIO_H__			/*  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include "HSP.h"
#define MAX_NAME_LEN 256
#define QUERY_LEN    256
#define FASTA        0
#define FASTQ        1

typedef struct _SEQ_T_{
	int max, l, ns;
	char name[MAX_NAME_LEN];
	char *seq, *rc, *qual;
}seq_t;

int CheckFast(int fd);
int fasta(FILE *fp, seq_t *seq, const int CONV);
int fastq(FILE *fp, seq_t *seq, const int CONV);

#endif     /* -----  __SEQIO_H__  ----- */
