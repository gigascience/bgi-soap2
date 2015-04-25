/*

   HSP.h		BWTBlastn functions

   This module contains miscellaneous BWTBlastn functions.

   Copyright (C) 2004, Wong Chi Kwong.

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

*/

#ifndef __HSP_H__
#define __HSP_H__

#include "TypeNLimit.h"
#include "MemManager.h"
#include "TextConverter.h"

#define ALPHABET_SIZE				4
#define BIT_PER_CHAR				2
#define CHAR_PER_128				64
#define CHAR_PER_WORD				16
#define CHAR_PER_BYTE				4

#define MAX_ALIGNMENT_LENGTH	131072
#define SHORTEST	70

typedef struct _ChrBlock{
	int chrID;
	unsigned int blockStart;
	unsigned int blockEnd;
	unsigned int ori;
}ChrBlock;

typedef struct _NewAnnotation{
	char chrName[MAX_SEQ_NAME_LENGTH];
	int nameLen;
	unsigned int chrStart;
	unsigned int chrEnd;
	int blockNum;
	ChrBlock *blockInChr;
}NewAnnotation;

typedef struct Annotation {
	int gi;
	char text[MAX_SEQ_NAME_LENGTH+1];
} Annotation;

typedef struct HSP {
	unsigned int* packedDNA;
	int chrNum;
	char **chrName;
	int numOfBlock;
	ChrBlock *blockList;
	unsigned int dnaLength;
}HSP;

#define MAX_SEQ_NAME_LENGTH				256

#define MAX_HISTO_SIZE					256

#define INVALID_CHAR_INDEX				15

#define ALIGN_MATCH					0
#define ALIGN_MISMATCH_AMBIGUITY	1
#define ALIGN_INSERT				2
#define ALIGN_DELETE				3

#define ALIGN_PER_WORD				16
#define ALIGN_BIT					2

#define AUX_TEXT_PER_WORD			8
#define AUX_TEXT_BIT				4

static const char lowercaseDnaCharIndex = 14;	// Seems that BLAST treat masked characters as 'N' (still have 1/4 chance of matching)
static const char nonMatchDnaCharIndex  = 15;
static const char dnaChar[16]			= {'A', 'C', 'G', 'T', 'M', 'R', 'S', 'V', 'W', 'Y', 'H', 'K', 'D', 'B', 'N', 'L'};
static const char dnaComplement[16]		= {'T', 'G', 'C', 'A', 'K', 'Y', 'S', 'B', 'W', 'R', 'D', 'M', 'H', 'V', 'N', 'L'};
static const char ambiguityCount[16]    = { 1 ,  1 ,  1 ,  1 ,  2 ,  2 ,  2 ,  3 ,  2 ,  2 ,  3 ,  2 ,  3 ,  3 ,  4 ,  0 };
static const char ambiguityMatch[16][4] = {{0, 0, 0, 0},
	{1, 0, 0, 0},
	{2, 0, 0, 0},
	{3, 0, 0, 0},
	{0, 1, 0, 0},
	{0, 2, 0, 0},
	{1, 2, 0, 0},
	{0, 1, 2, 0},
	{0, 3, 0, 0},
	{1, 3, 0, 0},
	{0, 1, 3, 0},
	{2, 3, 0, 0},
	{0, 2, 3, 0},
	{1, 2, 3, 0},
	{0, 1, 2, 3},
	{0, 0, 0, 0}
};

// Map must be allocated with char[256]
void HSPFillCharMap(unsigned char *charMap);
void HSPFillComplementMap(unsigned char *complementMap);

HSP *HSPLoad(MMPool *mmPool, const char *PackedDNAFileName, const char *AnnotationFileName);
HSP *HSPConvertFromText(MMPool *mmPool, const unsigned char *text, const unsigned int textLength,
						const unsigned int FASTARandomSeed, const int maskLowerCase,
						const int gi, const char *seqName);
void HSPFree(MMPool *mmPool, HSP *hsp);

unsigned int HSPParseFASTAToPacked(const char* FASTAFileName, const char* annotationFileName, const char* packedDNAFileName, const char* ambiguityFileName,
					  const unsigned int FASTARandomSeed, const int maskLowerCase);
unsigned int HSPPackedToFASTA(const char* FASTAFileName, const char* annotationFileName, const char* packedDNAFileName, const char* ambiguityFileName);


#endif
