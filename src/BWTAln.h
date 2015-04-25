/*
 * =============================================================================
 *
 *       Filename:  BWTAln.h
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


#ifndef  _BWTALN_H__INC
#define  _BWTALN_H__INC

#include "BWT.h"
#include "extratools.h"
#include "HSP.h"

unsigned int REVBWTForwardSearch(const unsigned char *convertedKey, const unsigned int keyLength, const BWT *rev_bwt,unsigned int *resultSaIndexLeft, unsigned int *resultSaIndexRight,unsigned int *rev_resultSaIndexLeft, unsigned int *rev_resultSaIndexRight);
unsigned int REVBWTContForwardSearch(const unsigned char *convertedKey, const unsigned int start, const unsigned int len,const BWT *rev_bwt,unsigned int *saL, unsigned int *saR,unsigned int *rev_saL, unsigned int *rev_saR);
unsigned int BWTContBackwardSearch(const unsigned char *convertedKey, const unsigned int start, const unsigned int len, const BWT *bwt, unsigned int *saL, unsigned int *saR);
unsigned int BWTBackward1Error(const unsigned char *querypattern, const BWTOPT *bo, BWT *bwt, unsigned int start, unsigned int len, unsigned int pl, unsigned int pr, unsigned int info, HITTABLE *hits);
unsigned int REVBWTForward1Error(const unsigned char *queryPattern,const BWTOPT *bo,  BWT *bwt, BWT * rev_bwt, unsigned int start,unsigned int len, unsigned int pl,unsigned int pr, unsigned int rev_pl,unsigned int rev_pr, unsigned int info, HITTABLE *hits);
int BWTExactMatching(const unsigned char *convertedKey, const BWTOPT *bo, const int chain, BWT *bwt, LOOKUPTABLE *lookup, HITTABLE *hits);
int BWT1ErrorMatching(const unsigned char *convertedKey, const BWTOPT *bo, const int chain, BWT *bwt, BWT *rev_bwt, LOOKUPTABLE *lookup, LOOKUPTABLE *rev_lookup, HITTABLE *hits);
int BWT2ErrorMatching(const unsigned char *convertedKey, const BWTOPT *bo, const int chain, BWT *bwt, BWT *rev_bwt, LOOKUPTABLE *lookup, LOOKUPTABLE *rev_lookup, HITTABLE *hits);
 int BWTGapMatching(const unsigned char *convertedKey, const BWTOPT *bo, const int chain, BWT *bwt, BWT *rev_bwt, LOOKUPTABLE *lookup, LOOKUPTABLE *rev_lookup, HITTABLE *hits);

#endif   /* ----- #ifndef _BWTALN_H__INC  ----- */



