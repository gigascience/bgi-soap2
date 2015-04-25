/*

   BWT.h	BWT-Index

   This module contains an implementation of BWT-index for alphabet size = 4.
   The functions provided include:
    Load functions for loading BWT to memory;
    Core functions for accessing core Inverse Psi values;
	Search functions for searching patterns from text;
	Text retrieval functions for retrieving text from BWT.

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

#ifndef __BWT_H__
#define __BWT_H__
#include "HSP.h"
#include "TypeNLimit.h"
#include "MemManager.h"
#include "TextConverter.h"


#define BITS_PER_OCC_VALUE			16
#define OCC_VALUE_PER_WORD			2
#define OCC_INTERVAL				256
#define WORD_BETWEEN_OCC			16
#define OCC_INTERVAL_MAJOR			65536

#define SORT_ALL					0
#define SORT_16_BIT					1
#define SORT_NONE					2

#define BUCKET_BIT					16
#define NUM_BUCKET					65536

#define MAX_APPROX_MATCH_ERROR	7
#define MAX_ARPROX_MATCH_LENGTH	32

#define BWTDP_MAX_SUBSTRING_LENGTH	512

typedef struct _BWTOPT_TYPE_{
	int cutoff;
	int alnLen, seqLen;
	int min_len;
	int h, x, y;
	int max_mm, gap_len, gap_fb;
	int nblock;
	ChrBlock *blockList;
	unsigned int *pacRef;
	unsigned int dnaLen;
	unsigned int extLen;
	char *fw, *rc;
}BWTOPT;

typedef struct SaIndexRange {
	unsigned int startSaIndex;
	unsigned int endSaIndex;
} SaIndexRange;


typedef struct BWT {
	unsigned int textLength;			// length of the text
	unsigned int saInterval;			// interval between two SA values stored explicitly
	unsigned int inverseSaInterval;		// interval between two inverse SA stored explicitly
	unsigned int inverseSa0;			// SA-1[0]
	unsigned int *cumulativeFreq;		// cumulative frequency
	unsigned int *bwtCode;				// BWT code
	unsigned int *occValue;				// Occurrence values stored explicitly
	unsigned int *occValueMajor;		// Occurrence values stored explicitly
	unsigned int *saValue;				// SA values stored explicitly
	unsigned int *inverseSa;			// Inverse SA stored explicitly
	SaIndexRange *saIndexRange;			// SA index range
	int saIndexRangeNumOfChar;			// Number of characters indexed in SA index range
	unsigned int *saValueOnBoundary;	// Pre-calculated frequently referred data
	unsigned int *decodeTable;			// For decoding BWT by table lookup
	unsigned int decodeTableGenerated;	// == TRUE if decode table is generated on load and will be freed
	unsigned int bwtSizeInWord;			// Temporary variable to hold the memory allocated
	unsigned int occSizeInWord;			// Temporary variable to hold the memory allocated
	unsigned int occMajorSizeInWord;	// Temporary variable to hold the memory allocated
	unsigned int saValueSize;			// Temporary variable to hold the memory allocated
	unsigned int inverseSaSize;			// Temporary variable to hold the memory allocated
	unsigned int saIndexRangeSize;		// Temporary variable to hold the memory allocated
} BWT;

#define MAX_DIAGONAL_LEVEL 4				// Number of sub-pattern to keep for detecting diagonal hit

// Error information is stored as:
// 1. bitVector
//	  After hamming distance match
// 2. count
//    After edit distance match
// 3. score
//    After the hits are processed with scoring functions

typedef struct SaIndexGroupNew {	// SA index range and information of a particular error arrangement of a matched sub-pattern
	unsigned int startSaIndex;			// starting SA index
	unsigned int numOfMatch;				// number of match
	unsigned int posQuery;				// position in query; used for detecting diagonal hits
	unsigned int info;					// extra hit information; to be copied to hitList.info
} SaIndexGroupNew;

typedef struct SaIndexGroupOld {	// SA index range and information of a particular error arrangement of a matched sub-pattern
	unsigned int startSaIndex;			// starting SA index
	unsigned int numOfMatch;				// number of match
	unsigned int info;					// extra hit information; to be copied to hitList.info
} SaIndexGroupOld;

typedef struct SaIndexGroup {	// SA index range and information of a particular error arrangement of a matched sub-pattern
	unsigned int startSaIndex;			// starting SA index
	unsigned int numOfMatch;			// number of match
	unsigned int info;					// extra hit information
} SaIndexGroup;

typedef struct SaIndexGroupWithErrorBitVector {	// SA index range and information of a particular error arrangement of a matched sub-pattern
	unsigned int startSaIndex;			// starting SA index
	unsigned int numOfMatch;			// number of match
	unsigned int errorBitVector;			// error bit vector
} SaIndexGroupWithErrorBitVector;

typedef struct SaIndexGroupWithLengthError {	// SA index range and information of a particular error arrangement of a matched sub-pattern
	unsigned int startSaIndex;			// starting SA index
	unsigned int numOfMatch;			// number of match
	unsigned posQuery : 16;		// position in query
	unsigned length   : 8;		// length of hit
	unsigned error    : 8;		// error in hit
} SaIndexGroupWithLengthError;

typedef struct SaIndexGroupProcessed {	// Alternative usage of SaIndexGroup - once processed, error bit vector is replaced by index to text position
	unsigned int startSaIndex;			// starting SA index
	unsigned int numOfMatch;			// number of match
	unsigned int textPositionIndex;		// storing the pointer to text position
} SaIndexGroupProcessed;

typedef struct DupSaIndexGroup {	// Alternative usage of SaIndexGroup - the group duplicates another group
	unsigned int lastDupSaIndexGroupIndex;	// index to last duplicated group
	unsigned int saIndexGroupIndex;			// index to the first SA into group among the duplicates
	unsigned int textPositionIndex;			// storing the pointer to text position
} DupSaIndexGroup;

typedef struct SaIndexGroupHash {	// Hash table for checking duplicate SA index group
	unsigned int startSaIndex;
	unsigned int saIndexGroupIndex;
} SaIndexGroupHash;

typedef struct BWTSaRetrievalStatistics {
	unsigned int bwtSaRetrieved;
	unsigned int saDiagonalLinked;
	unsigned int saDiagonalFiltered;
	unsigned int saDuplicated;
} BWTSaRetrievalStatistics;

typedef struct BWTDPStatistics {
	int maxDepth;
	int maxDPCell;
	int maxDPMemoryInWord;
	int totalMaxDepth;
	int totalMaxDPCell;
	int totalMaxDPMemoryInWord;
	LONG acceptedPathDepth;
	LONG acceptedPath;
	LONG rejectedPathDepth;
	LONG rejectedPath;
	LONG* __restrict totalNode;
	LONG* __restrict rejectedNode;
	LONG* __restrict totalDPCell;
} BWTDPStatistics;

typedef struct SaIndexList {
	unsigned int saIndex;
	unsigned int textPositionIndex;
} SaIndexList;

typedef struct HitCombination {
	int numOfCombination;
	int maxError;
	int keyLength;
	int skipTableWidth;
	int *errorPos;
	int *skip;
	int *skipErrorIndex;
} HitCombination;

typedef struct DPText {
	int charBeingProcessed;
	int dpCellIndex;
	int numOfDpCellSegment;
	unsigned int dummy1;	// Must not be removed; so that saIndexLeft and saIndexRight are aligned to 16 byte boundary
	unsigned int saIndexLeft[ALPHABET_SIZE];
	unsigned int saIndexRight[ALPHABET_SIZE];
} DPText;

typedef struct DPScanDepth {
	unsigned P				:	31;
	unsigned withAmbiguity	:	1;
} DPScanDepth;


// Load / unload functions
BWT *BWTCreate(MMPool *mmPool, const unsigned int textLength, unsigned int *decodeTable);
BWT *BWTLoad(MMPool *mmPool, const char *bwtCodeFileName, const char *occValueFileName, 
			 const char *saValueFileName, const char *inverseSaFileName, const char *saIndexRangeFileName,
			 unsigned int *decodeTable);
void BWTFree(MMPool *mmPool, BWT *bwt);
//void BWTPrintMemoryUsage(const BWT *bwt, FILE *output, const unsigned int packedDNASize);

// Precalculate frequenctly accessed data
void BWTGenerateSaValueOnBoundary(MMPool *mmPool, BWT *bwt);

// Core functions
// The following must be customized for differenet compression schemes ***
unsigned int BWTDecode(const BWT *bwt, const unsigned int index1, const unsigned int index2, const unsigned int character);
void BWTDecodeAll(const BWT *bwt, const unsigned int index1, const unsigned int index2, unsigned int* __restrict occValue);
unsigned int BWTOccValue(const BWT *bwt, unsigned int index, const unsigned int character);
void BWTOccValueTwoIndex(const BWT *bwt, unsigned int index1, unsigned int index2, const unsigned int character, unsigned int* __restrict occValue);
void BWTAllOccValue(const BWT *bwt, unsigned int index, unsigned int* __restrict occValue);
void BWTAllOccValueTwoIndex(const BWT *bwt, unsigned int index1, unsigned int index2, unsigned int* __restrict occValue1, unsigned int* __restrict occValue2);
unsigned int BWTOccValueOnSpot(const BWT *bwt, unsigned int index, unsigned int* __restrict character);
unsigned int BWTSearchOccValue(const BWT *bwt, const unsigned int character, const unsigned int searchOccValue);


// Utility functions for no compression only
unsigned int BWTResidentSizeInWord(const unsigned int numChar);
unsigned int BWTFileSizeInWord(const unsigned int numChar);
void BWTClearTrailingBwtCode(BWT *bwt);

// These are generic to different compression schemes (and generic to no compression as well)
unsigned int BWTPsiMinusValue(const BWT *bwt, const unsigned int index);
unsigned int BWTPsiPlusValue(const BWT *bwt, const unsigned int index);
unsigned int BWTSaValue(const BWT *bwt, unsigned int index);
unsigned int BWTInverseSa(const BWT *bwt, unsigned int saValue);
unsigned int BWTOccIntervalMajor(const unsigned int occInterval);
unsigned int BWTOccValueMinorSizeInWord(const unsigned int numChar);
unsigned int BWTOccValueMajorSizeInWord(const unsigned int numChar);

// Search functions
// packedText should be allocated with at least 1 Word buffer initialized to zero

// Text retrieval functions
// Position in text will be placed at the first word of hitListSizeInWord

// startSaIndex + resultInfo must be sorted in increasing order; there must be no overlapping groups except that one group can completely enclose another

// QSort comparison functions
int SaIndexGroupStartSaIndexOrder(const void *saIndexGroup, const int index1, const int index2);
int SaIndexGroupStartSaIndexLengthErrorOrder(const void *saIndexGroup, const int index1, const int index2);
int HitListPosTextErrorLengthOrder(const void *hitList, const int index1, const int index2);
int HitListPosText16BitOrder(const void *hitList, const int index1, const int index2);
int HitListPosTextOrder(const void *hitList, const int index1, const int index2);
int GappedHitListScorePosTextOrder(const void *gappedHitList, const int index1, const int index2);
int GappedHitListDbSeqIndexScorePosTextOrder(const void *gappedHitList, const int index1, const int index2);


#endif
