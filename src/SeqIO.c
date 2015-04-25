/*
 * =============================================================================
 *
 *       Filename:  SeqIO.c
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
#include "SeqIO.h"
extern unsigned char charMap[256];
extern unsigned char complementMap[256];
extern const char ambiguityCount[16];

int CheckFast (int fd) {
	char c;
	if (read(fd, &c, 1)>0) {
		lseek(fd, -1, SEEK_CUR);
		if (c == '>') return FASTA;
		else if (c == '@') return FASTQ;
		else {
			fprintf(stderr, "File Error: unrecognized file\n");
			close(fd);
			exit(EXIT_FAILURE);
		}
	}
	return FASTA;
}

int fasta(FILE *fp, seq_t *seq, const int CONV){

	int l, max, ns;
	int c;
	char *p, *q;
	
	l = ns = 0;
	max=seq->max;
	while (!feof(fp) && getc(fp)!= (int)'>');
	if (feof(fp)) return -1;

	p = seq->name;
	while ((c= getc(fp)) != ' ' && c != '\r' && c != '\t' && c != '\n' && ++l < MAX_NAME_LEN) *p++ = c;
	/*parse RG ID for SAM
	if(o->SAM){
		if(p[l-1] == '1' && p[l-2] == '/') {r->RG_ID = 1; l-=2;}
		else if(p[l-1] == '2' && p[l-2] == '/'){r->RG_ID = 2; l-=2;}
		else {seq->RG_ID = 1;}
	}
		//*/
	*p = '\0';

	while (c != '\n') c = (char) getc(fp);
	
	if (feof(fp)) {
		fprintf(stderr, "\nFile Error: unexpected feof\n");
		exit(EXIT_FAILURE);
	}

	l = 0;
	p=seq->seq; q=seq->rc;
	while ((c = getc(fp)) != '>' && !feof(fp)) {
		if (c != '\n' && c != '\r') {
			if (l >= max) {
				max += QUERY_LEN;
				seq->seq  = (char *)realloc(seq->seq, sizeof(char)* max);
				seq->rc   = (char *)realloc(seq->rc, sizeof(char)* max);
				seq->qual = (char *)realloc(seq->qual, sizeof(char) * max);
				p = seq->seq + l;
				q = seq->rc + l;
			}
			if (ambiguityCount[charMap[c]] != 1) {
				*p++ = charMap['G'];
				*q++ = complementMap['G'];
				ns++;
			} else {
				*p++ = charMap[c];
				*q++ = complementMap[c];
			}
			seq->qual[l] = 'h';
			l++;
		}
	}

	seq->qual[l] = *p   = *q = '\0';
	seq->l       = l;
	seq->max     = max;
	seq->ns      = ns;
	if (c == '>') ungetc(c,fp);
	return l;
}		/* -----  end of function fasta  ----- */

int fastq (FILE *fp, seq_t *seq, const int CONV) {

#ifdef DEBUG
//	fprintf(stderr, "get read\n");
#endif
	int l, max, ns;
	int c;
	char *p;
	l = ns = 0;
	max = seq->max;
	while (!feof(fp) && getc(fp)!= '@');
	if (feof(fp)) return -1;

	l = 0;
	p=seq->name;
	while ((c = getc(fp)) != '\t' && c != ' ' && c != '\n' && c != '\r' && l++ < MAX_NAME_LEN) *p++ = c;
	*p = '\0';
//	fprintf(stderr, "%s\n", seq->name);

	/* RG ID for SAM
	if(o->SAM){
		if(p[l-1] == '1' && p[l-2] == '/') {r->RG_ID = 1; l-=2;}
		else if(p[l-1] == '1' && p[l-2] == '/'){r->RG_ID = 2; l-=2;}
		else {r->RG_ID = 1;}
	}
		//*/

	while (c != '\n') c = getc(fp);

	if (feof(fp)) {
		fprintf(stderr, "\nFile Error: unexpected feof\n");
		exit(EXIT_FAILURE);
	}

	l = 0;
//	p = seq->seq; q = seq->rc;
	while ((c = getc(fp)) != '+' && !feof(fp)) {
		if (c != '\n' && c != '\r') {
			if (l >= max) {
				max += QUERY_LEN;
				seq->seq  = (char *)realloc(seq->seq, sizeof(char)*max);
				seq->rc   = (char *)realloc(seq->rc, sizeof(char)* max);
				seq->qual = (char *)realloc(seq->qual, sizeof(char)*max);
//				fprintf(stderr, "%d\n", max);
//				p = seq->seq + l;
//				q = seq->rc + l;
			}
//			fprintf(stdout, "%c", c);
			if(ambiguityCount[charMap[c]] == 1){
				seq->seq[l] = charMap[c];
				seq->rc[l++] = complementMap[c];
			}else{
				seq->seq[l] = charMap['G'];
				seq->rc[l++] = complementMap['G'];
				ns++;
			}
		}
	}
//	*p     = '\0'; *q = '\0';
//	fprintf(stderr, "\n");
//	for(j=0; j<l;j++)fprintf(stderr, "%d", seq->seq[j]);
//	fprintf(stderr, "\n");
	seq->l = l;
//	fprintf(stdout, "\n");
	while (!feof(fp) && (c= getc(fp))!= '\n');
	if (feof(fp)) {
		fprintf(stderr, "\nFile Error: unexpected feof\n");
		return 0;
	}

	l = 0;
	p = seq->qual;
	while ((c = (char) getc(fp)) != '\n' && c != '\r' && !feof(fp)) {
		if (l > max) {
			max += QUERY_LEN;
			seq->qual = (char *)realloc(seq->qual, sizeof(char)*max);
			p=seq->qual; p+=l;
		}
		*p++ = c;
		l++;
	}
	*p = '\0';

	if (l != seq->l) {
		fprintf(stderr, "Length Error: incompitable seq and qual length\n");
		fprintf(stderr, "       %s\n", seq->name);
		return 0;
	}
	if (c == '@') ungetc(c,fp);
	seq->max = max;
	seq->ns  = ns;
//	fprintf(stderr, "%d:%d\n", seq->l, l);
	return seq->l;
}

/* test
int main(int argc, char *argv[]){

	int fd = open(argv[1])

	return EXIT_SUCCESS;
}				
//*/
