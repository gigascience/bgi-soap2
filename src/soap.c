/* * =============================================================================
 *
 *       Filename:  soap.c
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <getopt.h>
#include "Match.h"
#include "HSP.h"
#include "BWT.h"
#include "TypeNLimit.h"
#include "extratools.h"
#include "MemManager.h"

#ifndef MAKE_TIME
#define MAKE_TIME "00:00:00"
#endif

#ifndef VID
//#define VID "2.11"
//#define VID "2.15"		// with gap and extend cigar
//#define VID "2.16"		// -r 2 segment fault
//#define VID "2.17"		// mm in rc
//#define VID "2.18"		// 07/05/2009
//#define VID "2.19"			// 13/07/2009
#define VID "2.20"			// 23/07/2009 gap missed in forward strand
#endif

#define MAX_FILENAME_LEN 1024
#define MAX_SUFFIX_LEN 255

const char *PROGRAM = "SOAPaligner/soap2";
const char *AUTHOR  = "BGI shenzhen";
const char *VERSION = VID;			/*release date: 14/01/2009*/
const char *CONTACT = "soap@genomics.org.cn";

char readAFileName[MAX_FILENAME_LEN]                           = "";
char readBFileName[MAX_FILENAME_LEN]                           = "";
char outFileName[MAX_FILENAME_LEN]                             = "";
char outUnpairFileName[MAX_FILENAME_LEN]                       = "";
char outUnmapFileName[MAX_FILENAME_LEN]                        = "";

char database_prefix[MAX_FILENAME_LEN]                         = "";

char AnnotationSuffix[MAX_SUFFIX_LEN]                          = ".ann";
char PackedDNASuffix[MAX_SUFFIX_LEN]                           = ".pac";
char BWTCodeSuffix[MAX_SUFFIX_LEN]                             = ".bwt";
char BWTOccValueSuffix[MAX_SUFFIX_LEN]                         = ".fmv";
char SaValueSuffix[MAX_SUFFIX_LEN]                             = ".sa";

char RevPackedDNASuffix[MAX_SUFFIX_LEN]                        = ".rev.pac";
char RevBWTCodeSuffix[MAX_SUFFIX_LEN]                          = ".rev.bwt";
char RevBWTOccValueSuffix[MAX_SUFFIX_LEN]                      = ".rev.fmv";

char LookupTableSuffix[MAX_SUFFIX_LEN]                         = ".lkt";
char RevLookupTableSuffix[MAX_SUFFIX_LEN]                      = ".rev.lkt";
char HighOccHashTableSuffix[MAX_SUFFIX_LEN]                    = ".hot";

// DatabaseFiles parameters
char AnnotationFileName[MAX_FILENAME_LEN+MAX_SUFFIX_LEN]       = "";
char PackedDNAFileName[MAX_FILENAME_LEN+MAX_SUFFIX_LEN]        = "";
char BWTCodeFileName[MAX_FILENAME_LEN+MAX_SUFFIX_LEN]          = "";
char BWTOccValueFileName[MAX_FILENAME_LEN+MAX_SUFFIX_LEN]      = "";
char SaValueFileName[MAX_FILENAME_LEN+MAX_SUFFIX_LEN]          = "";

//For Reversed BWT
char RevPackedDNAFileName[MAX_FILENAME_LEN+MAX_SUFFIX_LEN]     = "";
char RevBWTCodeFileName[MAX_FILENAME_LEN+MAX_SUFFIX_LEN]       = "";
char RevBWTOccValueFileName[MAX_FILENAME_LEN+MAX_SUFFIX_LEN]   = "";

//For Extra Data Structures
char LookupTableFileName[MAX_FILENAME_LEN+MAX_SUFFIX_LEN]      = "";
char RevLookupTableFileName[MAX_FILENAME_LEN+MAX_SUFFIX_LEN]   = "";
char HighOccHashTableFileName[MAX_FILENAME_LEN+MAX_SUFFIX_LEN] = "";
unsigned int LookUpTableSize  = 13;
unsigned int RevLookUpTableSize = 13;

// Memory parameters
	/*
	int PoolSize = 20971520;				// 2M  - fixed; not configurable through ini
	int WorkingMemorySize = 67108864;	// 64M - good for 8M hit; configurable through ini
	int AlignmentMemorySize = 4194304;	// 4M
//*/

int PoolSize            = 20971520;
int WorkingMemorySize   = 1073741824;
int AlignmentMemorySize = 536870912;

unsigned char charMap[256];
unsigned char complementMap[256];

static SOAPOPT *opt;

void Usage(void) {
	fprintf(stdout, "\nProgram: %s\n", PROGRAM);
	fprintf(stdout, "Compile Date: "
			MAKE_TIME"\n");
	fprintf(stdout, "Author:  %s\n", AUTHOR);
	fprintf(stdout, "Version: %s\n", VERSION);
	fprintf(stdout, "Contact: %s\n", CONTACT);
	fprintf(stdout, "\nUsage:\tsoap [options]\n");
	fprintf(stdout,	"\t-a  <str>   query a file, *.fq, *.fa\n");
	fprintf(stdout,	"\t-b  <str>   query b file\n");
	fprintf(stdout,	"\t-D  <str>   reference sequences indexing table, *.index format\n");
	fprintf(stdout,	"\t-o  <str>   output alignment file(txt)\n");
	fprintf(stdout, "\t-M  <int>   match mode for each read or the seed part of read, which shouldn't contain more than 2 mismaches, [4]\n"
			"\t            0: exact match only\n"
			"\t            1: 1 mismatch match only\n"
			"\t            2: 2 mismatch match only\n"
			"\t            4: find the best hits\n");
	fprintf(stdout,	"\t-u  <str>   output unmapped reads file\n");
	fprintf(stdout, "\t-t          output reads id instead reads name, [none]\n");
	fprintf(stdout, "\t-l  <int>   align the initial n bps as a seed [%d] means whole length of read\n", opt->aln_len);
	fprintf(stdout,	"\t-n  <int>   filter low-quality reads containing >n Ns before alignment, [%d]\n", opt->ns);
	fprintf(stdout,	"\t-r  [0,1,2] how to report repeat hits, 0=none; 1=random one; 2=all, [%d]\n", opt->rr);
	fprintf(stdout,	"\t-m  <int>   minimal insert size allowed, [%d]\n", opt->min_ins); //minimal insert size
	fprintf(stdout,	"\t-x  <int>   maximal insert size allowed, [%d]\n", opt->max_ins); //max_insert_size
	fprintf(stdout,	"\t-2  <str>   output file of unpaired alignment hits\n");
	fprintf(stdout,	"\t-v  <int>   maximum number of mismatches allowed on a read. [%d] bp\n", opt->max_mm);
	fprintf(stdout,	"\t-s  <int>   minimal alignment length (for soft clip) [%d] bp\n", opt->min_len);
//	fprintf(stdout,	"\t-U          only find uniq mapped reads with n mismatches for single-end, [%d]\n", opt->uniq);
//	fprintf(stdout,	"\t-A          report all mismatches reads in SOAP Format, default [none] report number \n");
	fprintf(stdout,	"\t-g  <int>   one continuous gap size allowed on a read. [%d] bp\n", opt->gap_len );//max_gap, allowed_gap
	fprintf(stdout,	"\t-R          for long insert size of pair end reads RF. [none](means FR pair)\n");
	fprintf(stdout,	"\t-e  <int>   will not allow gap exist inside n-bp edge of a read, default=5\n");     //gap_edge
//fprintf(stdout,	"\t\t-z  <char>  initial quality, default=@ [Illumina is using '@', Sanger Institute is using '!']\n");//zero_quality
//	fprintf(stdout,	"\t\t-c  [0,1,2] do alignment on which reference chain? 0:both; 1:forward only; 2:reverse only. default=%d");//chains
#ifdef PTHREADS
	fprintf(stdout,	"\t-p  <int>   number of processors to use, [%d]\n", opt->nthreads);   //number of processors
#endif
	fprintf(stdout,	"\n\t-h          this help\n\n");
	exit(1);
	//*/
}		/* -----  end of function Usage  ----- */

SOAPOPT *OptIni(){
	SOAPOPT *o;
	o = (SOAPOPT *) malloc (1 * sizeof(SOAPOPT));
	o->fast     = FASTQ;
	o->aln_len  = 256;
	o->ns       = 5;
	o->max_mm   = 5;
	o->gap_len  = 5;
	o->gap_fb   = 5;
	o->nthreads = 1;
	o->min_ins  = 400;
	o->max_ins  = 600;
	o->unmapped = 0;
	o->rr       = 1;
	o->gap_len  = 0;
	o->pe       = FALSE;
	o->cutoff   = 1;
	o->mode     = 4;
	o->id       = FALSE;
	o->FR       = TRUE;
	o->allErr   = FALSE;
	o->min_len  = 255;
	o->uniq     = 0;
	return o;
}

void ParseOpt(int argc, char *argv[]){
	char c;
	while((c = getopt(argc, argv, "a:b:D:o:2:u:m:x:M:AK:l:v:U:g:w:i:e:q:c:Rz:r:B:s:p:tn:h"))!=-1){
		switch(c){
			//basic IO
			case 'a':
				snprintf(readAFileName, MAX_FILENAME_LEN, "%s", optarg);
				break;
			case 'D':
				snprintf(database_prefix, MAX_FILENAME_LEN, "%s", optarg);
				break;
			case 'o':
				snprintf(outFileName, MAX_FILENAME_LEN, "%s", optarg);
				break;
			case 'b':
				opt->pe = TRUE;
				snprintf(readBFileName, MAX_FILENAME_LEN, "%s", optarg);
				break;
			case '2':
				snprintf(outUnpairFileName, MAX_FILENAME_LEN, "%s", optarg);
				break;
			case 'm':
				opt->min_ins = atoi(optarg);
				break;
			case 'x':
				opt->max_ins = atoi(optarg);
				break;

//advance options
			case 'u':
				opt->unmapped = TRUE;
				snprintf(outUnmapFileName, MAX_FILENAME_LEN, "%s", optarg);
				break;
			case 'l':
				opt->aln_len = atoi(optarg);
				break;
			case 'M':
				{
					opt->mode = atoi(optarg);
					
					if(opt->mode == 4) opt->cutoff = 1;
					else if(opt->mode == 5) opt->cutoff = MAX_ALN;
					break;
				}
			case 'K':
				opt->cutoff = min(atoi(optarg), MAX_ALN);
				break;
			case 'v':
				opt->max_mm = min(atoi(optarg), MAX_MISMATCH);
				break;
			case 'g':
				opt->gap_len = min(atoi(optarg), MAX_GAP_LEN);
				break;
			case 'e':
				opt->gap_fb = atoi(optarg);
				break;
			case 'R':
				opt->FR = 0;
				break;
			case 'z':
				opt->zero_qual = atoi(optarg);
				break;
			case 'r':
				opt->rr = atoi(optarg);
				break;
			case 't':
				opt->id = 1;
				break;
			case 'n':
				opt->ns = atoi(optarg);
				break;
			case 'B':
				opt->bisulfite = atoi(optarg);
				break;
			case 'U':
				opt->uniq = atoi(optarg);
				break;
			case 's':
				opt->min_len = atoi(optarg);
				break;
#ifdef PTHREADS
#define MAX_PTHREADS 20
			case 'p':
				opt->nthreads = min(atoi(optarg), MAX_PTHREADS);
				break;
#endif
			case 'c':
				opt->chain = atoi(optarg);
				break;
			//unrecognizable input
			case 'h':
			case '?':
				Usage();
		}
	}
}

void FileNameIni(){
	snprintf(AnnotationFileName, MAX_FILENAME_LEN+MAX_SUFFIX_LEN, "%s%s", database_prefix, AnnotationSuffix);
	snprintf(PackedDNAFileName, MAX_FILENAME_LEN+MAX_SUFFIX_LEN, "%s%s", database_prefix, PackedDNASuffix);
	snprintf(BWTCodeFileName, MAX_FILENAME_LEN+MAX_SUFFIX_LEN, "%s%s", database_prefix, BWTCodeSuffix);
	snprintf(BWTOccValueFileName, MAX_FILENAME_LEN+MAX_SUFFIX_LEN, "%s%s", database_prefix, BWTOccValueSuffix);
	snprintf(SaValueFileName, MAX_FILENAME_LEN+MAX_SUFFIX_LEN, "%s%s", database_prefix, SaValueSuffix);
	snprintf(RevPackedDNAFileName, MAX_FILENAME_LEN+MAX_SUFFIX_LEN, "%s%s", database_prefix, RevPackedDNASuffix);
	snprintf(RevBWTCodeFileName, MAX_FILENAME_LEN+MAX_SUFFIX_LEN, "%s%s", database_prefix, RevBWTCodeSuffix);
	snprintf(RevBWTOccValueFileName, MAX_FILENAME_LEN+MAX_SUFFIX_LEN, "%s%s",database_prefix , RevBWTOccValueSuffix);
	snprintf(LookupTableFileName, MAX_FILENAME_LEN+MAX_SUFFIX_LEN, "%s%s", database_prefix, LookupTableSuffix);
	snprintf(RevLookupTableFileName, MAX_FILENAME_LEN+MAX_SUFFIX_LEN, "%s%s", database_prefix, RevLookupTableSuffix);
	snprintf(HighOccHashTableFileName, MAX_FILENAME_LEN+MAX_SUFFIX_LEN, "%s%s", database_prefix, HighOccHashTableSuffix);
}
#define MODE S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH
void FileTest(FILEDS *fds){
	if((fds->ifdA=open(readAFileName, O_RDONLY))==-1){
		fprintf(stderr, "Query File Error: Can't read  %s\n", readAFileName);
		exit(EXIT_FAILURE);
	}
	if((fds->ofdAln=creat(outFileName, MODE))==-1){
		fprintf(stderr, "Output File Error: Can't write %s\n",outFileName);
		exit(EXIT_FAILURE);
	}
	if(opt->pe){
		if((fds->ifdB=open(readBFileName, O_RDONLY))==-1){
			fprintf(stderr, "Query File Error: Can't read  %s\n", readBFileName);
			exit(EXIT_FAILURE);
		}
		if( (fds->ofdSe=creat(outUnpairFileName, MODE))==-1){
			fprintf(stderr, "Output File Error: Can't write %s\n", outUnpairFileName);
			exit(EXIT_FAILURE);
		}
	}
	fprintf(stderr, "Query File a: %s\n", readAFileName);
	if (opt->pe)fprintf(stderr, "Query File b: %s\n", readBFileName);
	fprintf(stderr, "Output File: %s\n", outFileName);
	if (opt->pe)fprintf(stderr, "             %s\n", outUnpairFileName);
	if(opt->unmapped){
		if ((fds->ofdUn=creat(outUnmapFileName, MODE))==-1){
			fprintf(stderr, "Output File Error: Can't write %s\n", outUnmapFileName);
			exit(EXIT_FAILURE);
		} else {
			fprintf(stderr, "             %s\n", outUnmapFileName);
		}
	}
}

int main(int argc, char *argv[]){

	opt = OptIni();
	if (argc < 3) {
		Usage();
	}
	
	HSP *hsp;

	BWT *bwt;
	BWT *rev_bwt;

	LOOKUPTABLE lookup;
	LOOKUPTABLE rev_lookup;

	HASHTABLE hashtable;

	MMPool *mmPool;

	MMMasterInitialize(3, 0, FALSE, NULL);
	mmPool = MMPoolCreate(PoolSize);

	HSPFillCharMap(charMap);
	HSPFillComplementMap(complementMap);
	asciiTime("Begin Program SOAPaligner/soap2");
	double startTime = setStartTime();
	double elapsedTime = 0;
	double loadTime=0;
	
	ParseOpt(argc, argv);
	FileNameIni();
	fprintf(stderr, "Reference: %s\n", database_prefix);
	FILEDS fds;
	FileTest(&fds);
	fprintf(stderr, "Load Index Table ...\n");
//*
	hsp = HSPLoad(mmPool, PackedDNAFileName, AnnotationFileName);
	bwt = BWTLoad(mmPool, BWTCodeFileName, BWTOccValueFileName, SaValueFileName, NULL, NULL, NULL);
	rev_bwt = BWTLoad(mmPool, RevBWTCodeFileName, RevBWTOccValueFileName, NULL, NULL, NULL, NULL);
	LoadLookupTable(&lookup,LookupTableFileName,LookUpTableSize);
	LoadLookupTable(&rev_lookup,RevLookupTableFileName,RevLookUpTableSize);
	LoadHashTable(&hashtable,HighOccHashTableFileName);
	RegisterDecoder(bwt,&hashtable);
	loadTime = getElapsedTime(startTime);
	//*/
	fprintf(stderr, "Load Index Table OK\n");
	fprintf(stderr, "Begin Alignment ...\n");
	MatchProcess(&fds, bwt, rev_bwt, &lookup, &rev_lookup, hsp, opt);
	//*/
	elapsedTime = getElapsedTime(startTime);
	fprintf(stderr, "Total Elapsed Time:       %9.2f\n"
			"      - Load Index Table: %9.2f\n"
			"      - Alignment:        %9.2f\n", elapsedTime, loadTime, (elapsedTime-loadTime));

	FreeLookupTable(&lookup);
	FreeLookupTable(&rev_lookup);
	FreeHashTable(&hashtable);
	HSPFree(mmPool, hsp);
	BWTFree(mmPool, bwt);
	BWTFree(mmPool, rev_bwt);
	MMPoolFree(mmPool);
	close(fds.ifdA);close(fds.ofdAln);
	if(opt->pe){close(fds.ifdB); close(fds.ofdSe);}
	if(opt->unmapped)close(fds.ofdUn);
	free(opt);

	asciiTime("SOAPaligner/soap2 End");
	fprintf(stderr, "\n");
	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
