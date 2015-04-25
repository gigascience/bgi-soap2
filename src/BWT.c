/*

   BWT.c	BWT-Index

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
   GNU General Public License for more details.L

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <emmintrin.h>
#include <mmintrin.h>
#include "BWT.h"
#include "MiscUtilities.h"
#include "DNACount.h"
#include "TextConverter.h"
#include "MemManager.h"
#include "r250.h"
#include "HSP.h"

// static functions
static INLINE unsigned int BWTOccValueExplicit(const BWT *bwt, const unsigned int occIndexExplicit, const unsigned int character);
static INLINE void BWTAllOccValueExplicit(const BWT *bwt, const unsigned int occIndexExplicit, unsigned int* __restrict occValueExplicit);
static INLINE unsigned int BWTSaIndexToChar(const BWT *bwt, const unsigned int saIndex);
static INLINE unsigned int BWTGetWordPackedText(const unsigned int *packedText, const unsigned int index, const unsigned int shift, const unsigned int numOfBit);

static INLINE void BWTPrefetchOccValueExplicit(const BWT *bwt, const unsigned int occIndexExplicit);
static INLINE void BWTPrefetchBWT(const BWT *bwt, const unsigned int index);


int SaIndexGroupDPHitOrder1(const void *saIndexGroup, const int index1, const int index2);
int SaIndexGroupDPHitOrder2(const void *saIndexGroup, const int index1, const int index2);


static INLINE unsigned int BWTSaIndexToChar(const BWT *bwt, const unsigned int saIndex) {

	return (saIndex > bwt->cumulativeFreq[1]) + (saIndex > bwt->cumulativeFreq[2])
										   + (saIndex > bwt->cumulativeFreq[3]);

}

BWT *BWTCreate(MMPool *mmPool, const unsigned int textLength, unsigned int *decodeTable) {

	BWT *bwt;

	bwt = MMPoolDispatch(mmPool, sizeof(BWT));

	bwt->textLength = 0;
	bwt->inverseSa = 0;

	bwt->cumulativeFreq = MMPoolDispatch(mmPool, (ALPHABET_SIZE + 1) * sizeof(unsigned int));
	initializeVAL(bwt->cumulativeFreq, ALPHABET_SIZE + 1, 0);

	bwt->bwtSizeInWord = 0;
	bwt->saValueOnBoundary = NULL;

	// Generate decode tables
	if (decodeTable == NULL) {
		bwt->decodeTable = MMPoolDispatch(mmPool, DNA_OCC_CNT_TABLE_SIZE_IN_WORD * sizeof(unsigned int));
		GenerateDNAOccCountTable(bwt->decodeTable);
	} else {
		bwt->decodeTable = decodeTable;
	}

	bwt->occMajorSizeInWord = BWTOccValueMajorSizeInWord(textLength);
	bwt->occValueMajor = MMPoolDispatch(mmPool, bwt->occMajorSizeInWord * sizeof(unsigned int));

	bwt->occSizeInWord = 0;
	bwt->occValue = NULL;

	bwt->saInterval = ALL_ONE_MASK;
	bwt->saValueSize = 0;
	bwt->saValue = NULL;

	bwt->inverseSaInterval = ALL_ONE_MASK;
	bwt->inverseSaSize = 0;
	bwt->inverseSa = NULL;

	return bwt;

}

BWT *BWTLoad(MMPool *mmPool, const char *bwtCodeFileName, const char *occValueFileName, 
			 const char *saValueFileName, const char *inverseSaFileName, const char *saIndexRangeFileName,
			 unsigned int *decodeTable) {

	unsigned int i;
	FILE *bwtCodeFile, *occValueFile, *saValueFile = NULL, *inverseSaFile = NULL, *saIndexRangeFile = NULL;
	BWT *bwt;
	unsigned int tmp;
	unsigned int bwtCodeLengthInFile;
	unsigned int numOfSaIndexRange;

	bwtCodeFile = (FILE*)fopen64(bwtCodeFileName, "rb");
	if (bwtCodeFile == NULL) {
		fprintf(stderr, "BWTLoad() : cannot open bwtCodeFile!\n");
		exit(1);
	}

	occValueFile = (FILE*)fopen64(occValueFileName, "rb");
	if (occValueFile == NULL) {
		fprintf(stderr, "BWTLoad() : cannot open occValueFile!\n");
		exit(1);
	}

	if (saValueFileName != NULL && saValueFileName[0] != '\0' && saValueFileName[0] != '-') {
		saValueFile = (FILE*)fopen64(saValueFileName, "rb");
		if (saValueFile == NULL) {
			fprintf(stderr, "BWTLoad() : cannot open saValueFile!\n");
			exit(1);
		}
	}

	if (inverseSaFileName != NULL && inverseSaFileName[0] != '\0' && inverseSaFileName[0] != '-') {
		inverseSaFile = (FILE*)fopen64(inverseSaFileName, "rb");
		if (inverseSaFile == NULL) {
			fprintf(stderr, "BWTLoad() : cannot open inverseSaFile!\n");
			exit(1);
		}
	}

	if (saIndexRangeFileName != NULL && saIndexRangeFileName[0] != '\0' && saIndexRangeFileName[0] != '-') {
		saIndexRangeFile = (FILE*)fopen64(saIndexRangeFileName, "rb");
		if (saIndexRangeFile == NULL) {
			fprintf(stderr, "BWTLoad() : cannot open saIndexRangeFile!\n");
			exit(1);
		}
	}

	bwt = MMPoolDispatch(mmPool, sizeof(BWT));

	fread(&bwt->inverseSa0, sizeof(unsigned int), 1, bwtCodeFile);

	bwt->cumulativeFreq = MMPoolDispatch(mmPool, (ALPHABET_SIZE + 1) * sizeof(unsigned int));
	bwt->cumulativeFreq[0] = 0;
	fread(bwt->cumulativeFreq + 1, sizeof(unsigned int), ALPHABET_SIZE, bwtCodeFile);
	bwt->textLength = bwt->cumulativeFreq[ALPHABET_SIZE];

	fread(&tmp, sizeof(unsigned int), 1, occValueFile);
	if (tmp != bwt->inverseSa0) {
		fprintf(stderr, "BWTLoad(): OccValue inverseSa0 not match!\n");
		exit(1);
	}
	for (i=1; i<=ALPHABET_SIZE; i++) {
		fread(&tmp, sizeof(unsigned int), 1, occValueFile);
		if (tmp != bwt->cumulativeFreq[i]) {
			fprintf(stderr, "BWTLoad(): OccValue cumulativeFreq not match!\n");
			exit(1);
		}
	}

	bwt->bwtSizeInWord = BWTResidentSizeInWord(bwt->textLength) + WORD_BETWEEN_OCC / 2;	// + 8 words so that the 128 bits before and after an explicit occ are in the same aligned 64 byte
	bwtCodeLengthInFile = BWTFileSizeInWord(bwt->textLength);
	bwt->bwtCode = MMUnitAllocate(bwt->bwtSizeInWord * sizeof(unsigned int));
	fread(bwt->bwtCode, sizeof(unsigned int), bwtCodeLengthInFile, bwtCodeFile);
	fclose(bwtCodeFile);
	BWTClearTrailingBwtCode(bwt);

	bwt->occSizeInWord = BWTOccValueMinorSizeInWord(bwt->textLength) ;
	bwt->occMajorSizeInWord = BWTOccValueMajorSizeInWord(bwt->textLength);
	bwt->occValue = MMUnitAllocate(bwt->occSizeInWord * sizeof(unsigned int));
	fread(bwt->occValue, sizeof(unsigned int), bwt->occSizeInWord, occValueFile);
	bwt->occValueMajor = MMUnitAllocate(bwt->occMajorSizeInWord * sizeof(unsigned int));
	fread(bwt->occValueMajor, sizeof(unsigned int), bwt->occMajorSizeInWord, occValueFile);
	fclose(occValueFile);

	if (decodeTable == NULL) {
		bwt->decodeTable = MMUnitAllocate(DNA_OCC_CNT_TABLE_SIZE_IN_WORD * sizeof(unsigned int));
		GenerateDNAOccCountTable(bwt->decodeTable);
		bwt->decodeTableGenerated = TRUE;
	} else {
		bwt->decodeTable = decodeTable;
		bwt->decodeTableGenerated = FALSE;
	}

	bwt->saValueOnBoundary = NULL;
	if (saValueFile == NULL) {
		bwt->saInterval = ALL_ONE_MASK;
		bwt->saValueSize = 0;
		bwt->saValue = NULL;
	} else {
		fread(&tmp, sizeof(unsigned int), 1, saValueFile);
		if (tmp != bwt->inverseSa0) {
			fprintf(stderr, "BWTLoad(): SaValue inverseSa0 not match!\n");
			exit(1);
		}
		for (i=1; i<=ALPHABET_SIZE; i++) {
			fread(&tmp, sizeof(unsigned int), 1, saValueFile);
			if (tmp != bwt->cumulativeFreq[i]) {
				fprintf(stderr, "BWTLoad(): SaValue cumulativeFreq not match!\n");
				exit(1);
			}
		}
		fread(&bwt->saInterval, sizeof(unsigned int), 1, saValueFile);
		bwt->saValueSize = (bwt->textLength + bwt->saInterval) / bwt->saInterval * sizeof(unsigned int);
		bwt->saValue = MMUnitAllocate(bwt->saValueSize);
		fread(bwt->saValue, 1, bwt->saValueSize, saValueFile);
		bwt->saValue[0] = (unsigned int)-1;	// Special handling for bwt
		fclose(saValueFile);

		BWTGenerateSaValueOnBoundary(mmPool, bwt);
	}

	if (inverseSaFile == NULL) {
		bwt->inverseSaInterval = ALL_ONE_MASK;
		bwt->inverseSaSize = 0;
		bwt->inverseSa = NULL;
	} else {
		fread(&tmp, sizeof(unsigned int), 1, inverseSaFile);
		if (tmp != bwt->inverseSa0) {
			fprintf(stderr, "BWTLoad(): InverseSaValue inverseSa0 not match!\n");
			exit(1);
		}
		for (i=1; i<=ALPHABET_SIZE; i++) {
			fread(&tmp, sizeof(unsigned int), 1, inverseSaFile);
			if (tmp != bwt->cumulativeFreq[i]) {
				fprintf(stderr, "BWTLoad(): InverseSaValue cumulativeFreq not match!\n");
				exit(1);
			}
		}
		fread(&bwt->inverseSaInterval, sizeof(unsigned int), 1, inverseSaFile);
		bwt->inverseSaSize = (bwt->textLength + bwt->inverseSaInterval) / bwt->inverseSaInterval * sizeof(unsigned int);
		bwt->inverseSa = MMUnitAllocate(bwt->inverseSaSize);
		fread(bwt->inverseSa, 1, bwt->inverseSaSize, inverseSaFile);
		fclose(inverseSaFile);
	}

	// Load Sa index range
	if (saIndexRangeFile == NULL) {
		bwt->saIndexRange = NULL;
		bwt->saIndexRangeNumOfChar = 0;
		bwt->saIndexRangeSize = 0;
	} else {
		fread(&tmp, sizeof(unsigned int), 1, saIndexRangeFile);
		if (tmp != bwt->inverseSa0) {
			fprintf(stderr, "BWTLoad(): SaIndex inverseSa0 not match!\n");
			exit(1);
		}
		for (i=1; i<=ALPHABET_SIZE; i++) {
			fread(&tmp, sizeof(unsigned int), 1, saIndexRangeFile);
			if (tmp != bwt->cumulativeFreq[i]) {
				fprintf(stderr, "BWTLoad(): SaIndex cumulativeFreq not match!\n");
				exit(1);
			}
		}
		fread(&bwt->saIndexRangeNumOfChar, sizeof(unsigned int), 1, saIndexRangeFile);
		numOfSaIndexRange = 1 << (bwt->saIndexRangeNumOfChar * 2);	// 4^saIndexRangeNumOfChar
		bwt->saIndexRange = MMUnitAllocate(numOfSaIndexRange * sizeof(SaIndexRange));
		fread(bwt->saIndexRange, sizeof(SaIndexRange), numOfSaIndexRange, saIndexRangeFile);
		bwt->saIndexRangeSize = numOfSaIndexRange * sizeof(SaIndexRange);
		fclose(saIndexRangeFile);
	}

	return bwt;

}


void BWTFree(MMPool *mmPool, BWT *bwt) {

	MMPoolReturn(mmPool, bwt->cumulativeFreq, ALPHABET_SIZE * sizeof(unsigned int));
	MMUnitFree(bwt->bwtCode, bwt->bwtSizeInWord * sizeof(unsigned int));

	if (bwt->occValue != NULL) {
		MMUnitFree(bwt->occValue, bwt->occSizeInWord * sizeof(unsigned int));
	}
	if (bwt->occValueMajor != NULL) {
		MMUnitFree(bwt->occValueMajor, bwt->occMajorSizeInWord * sizeof(unsigned int));
	}

	if (bwt->saValue != NULL) {
		MMUnitFree(bwt->saValue, bwt->saValueSize);
	}
	if (bwt->inverseSa != NULL) {
		MMUnitFree(bwt->inverseSa, bwt->inverseSaSize);
	}

	if (bwt->decodeTableGenerated == TRUE) {
		MMUnitFree(bwt->decodeTable, DNA_OCC_CNT_TABLE_SIZE_IN_WORD * sizeof(unsigned int));
	}

	if (bwt->saIndexRange != NULL) {
		MMUnitFree(bwt->saIndexRange, bwt->saIndexRangeSize);
	}

	if (bwt->saValueOnBoundary != NULL) {
		MMPoolReturn(mmPool, bwt->saValueOnBoundary, sizeof(unsigned int) * 2 * ALPHABET_SIZE);
	}

	MMPoolReturn(mmPool, bwt, sizeof(BWT));

}
/*
void BWTPrintMemoryUsage(const BWT *bwt, FILE *output, const unsigned int packedDNASize) {

	unsigned int totalMemorySize;

	fprintf(output, "BWT code size    : %u\n", bwt->bwtSizeInWord * sizeof(unsigned int));
	fprintf(output, "Occ value size   : %u\n", (bwt->occSizeInWord + bwt->occMajorSizeInWord) * sizeof(unsigned int));
	if (bwt->saValueSize > 0) {
		fprintf(output, "SA value size    : %u\n", bwt->saValueSize);
	}
	if (bwt->inverseSaSize > 0) {
		fprintf(output, "Inverse SA size  : %u\n", bwt->inverseSaSize);
	}
	if (bwt->saIndexRange > 0) {
		fprintf(output, "SA index rangee  : %u\n", bwt->saIndexRangeSize);
	}
	if (packedDNASize > 0) {
		fprintf(output, "Packed DNA size  : %u\n", packedDNASize);
	}
	
	totalMemorySize = (bwt->bwtSizeInWord + bwt->occSizeInWord + bwt->occMajorSizeInWord) * sizeof(unsigned int)
					   + bwt->saValueSize + bwt->inverseSaSize + bwt->saIndexRangeSize + packedDNASize;
	fprintf(output, "Total memory     : %u\n", totalMemorySize);
	fprintf(output, "Bit per char     : %.2f\n", 
			(float)totalMemorySize / ((float)bwt->textLength / BITS_IN_BYTE));

}
//*/
void BWTGenerateSaValueOnBoundary(MMPool *mmPool, BWT *bwt) {

	unsigned int i;

	if (bwt->saValueOnBoundary == NULL) {
		bwt->saValueOnBoundary = MMPoolDispatch(mmPool, sizeof(unsigned int) * 2 * ALPHABET_SIZE);
	}

	for (i=0; i<ALPHABET_SIZE; i++) {
		bwt->saValueOnBoundary[i * 2 + 1] = BWTSaValue(bwt, bwt->cumulativeFreq[i + 1]);
		if (bwt->cumulativeFreq[i] < bwt->textLength) {
			bwt->saValueOnBoundary[i * 2] = BWTSaValue(bwt, bwt->cumulativeFreq[i] + 1);
		} else {
			bwt->saValueOnBoundary[i * 2] = bwt->saValueOnBoundary[i * 2 + 1];
		}
	}

}

// Ordering of index1 and index2 is not important; this module will handle the ordering
// index1 and index2 can be on the same aligned 128 bit region or can be on adjacant aligned 128 bit region
// If index1 and index2 are in the same aligned 128 bit region, one of them must be on the boundary
// These requirements are to reduce the no. of branches in the program flow

unsigned int BWTDecode(const BWT *bwt, const unsigned int index1, const unsigned int index2, const unsigned int character) {

	unsigned int numChar1, numChar2, minIndex, maxIndex, minIndex128, maxIndex128;
	unsigned int r;

	const static unsigned int ALIGN_16 partitionOne1[4]  = { 47, 31, 15, 0 };
	const static unsigned int ALIGN_16 partitionOne2[4]  = { 0, 15, 31, 47 };
	const static unsigned int ALIGN_16 partitionZero1[4]  = { 63, 47, 31, 15 };
	const static unsigned int ALIGN_16 partitionZero2[4]  = { 15, 31, 47, 63 };

	// SSE registers
	__m128i r1e, r2e;
	__m128i mcl;
	__m128i m0, m1;
	__m128i r1a, r1b, r1c;
	__m128i r2a, r2b, r2c;

	// Sort index1 and index2
	r = (index1 - index2) & -(index1 < index2);
	minIndex = index2 + r;
	maxIndex = index1 - r;

	// Locate 128 bit boundary
	minIndex128 = lastAlignedBoundary(minIndex, CHAR_PER_128);
	maxIndex128 = lastAlignedBoundary(maxIndex - (maxIndex - minIndex > CHAR_PER_128), CHAR_PER_128);

	// Determine no.of characters to count
	numChar1 = maxIndex128 - minIndex;
	numChar2 = maxIndex - maxIndex128;

	// Load encoding into register here in the hope of hiding some memory latency
	r1e = _mm_load_si128((__m128i *)(bwt->bwtCode + minIndex128 / CHAR_PER_WORD));	// Load encoding into register
	r2e = _mm_load_si128((__m128i *)(bwt->bwtCode + maxIndex128 / CHAR_PER_WORD));	// Load encoding into register

	// Set character extraction masks 
	m0 = _mm_set1_epi32(0xFFFFFFFF + (character & 1));	// Character selection mask for even bits
	m1 = _mm_set1_epi32(0xFFFFFFFF + (character >> 1));	// Character selection mask for odd bits
	mcl = _mm_set1_epi32(0x55555555);					// Set bit-clearing mask to 0x55555555....(alternate 1-bit)

	// Set counting mask for 2 x 128 bits

	r1a = _mm_set1_epi32(numChar1);		// Load number of characters into register
	r2a = _mm_set1_epi32(numChar2);		// Load number of characters into register

	r1b = _mm_load_si128((__m128i*)partitionOne1);	// Load partition into register
	r2b = _mm_load_si128((__m128i*)partitionOne2);	// Load partition into register

	r1c = _mm_load_si128((__m128i*)partitionZero1);	// Load partition into register
	r2c = _mm_load_si128((__m128i*)partitionZero2);	// Load partition into register

	r1b = _mm_cmpgt_epi32(r1a, r1b);				// Compare to generate 4x32 bit mask; the word with counting boundary is all ones
	r2b = _mm_cmpgt_epi32(r2a, r2b);				// Compare to generate 4x32 bit mask; the word with counting boundary is all ones
				
	r1c = _mm_cmpgt_epi32(r1a, r1c);				// Compare to generate 4x32 bit mask; the word with counting boundary is all zeros
	r2c = _mm_cmpgt_epi32(r2a, r2c);				// Compare to generate 4x32 bit mask; the word with counting boundary is all zeros

	r1b = _mm_srli_epi32(r1b, (16 - numChar1 % 16) * 2);	// Shift bits so that all word comform to the requirement of counting the word with counting boundary 
	r2b = _mm_slli_epi32(r2b, (16 - numChar2 % 16) * 2);	// Shift bits so that all word comform to the requirement of counting the word with counting boundary

	r1c = _mm_or_si128(r1b, r1c);	// Combine two masks
	r2c = _mm_or_si128(r2b, r2c);	// Combine two masks

	r1c = _mm_and_si128(r1c, mcl);	// Combine with bit-clearing mask (now = 0x55555555....)
	r2c = _mm_and_si128(r2c, mcl);	// Combine with bit-clearing mask (now = 0x55555555....)

	// Start counting; encoding has been loaded into register earlier

	r1b = _mm_srli_epi32(r1e, 1);	// Shift encoding to right by 1 bit
	r2b = _mm_srli_epi32(r2e, 1);	// Shift encoding to right by 1 bit

	r1a = _mm_xor_si128(r1e, m0);	// Check even-bits with mask
	r2a = _mm_xor_si128(r2e, m0);	// Check even-bits with mask

	r1b = _mm_xor_si128(r1b, m1);	// Check odd-bits with mask
	r2b = _mm_xor_si128(r2b, m1);	// Check odd-bits with mask

	r1a = _mm_and_si128(r1a, r1b);	// Combine even and odd bits
	r2a = _mm_and_si128(r2a, r2b);	// Combine even and odd bits

	r1a = _mm_and_si128(r1a, r1c);	// Combine with counting mask, which has been combined with bit-clearing mask of 0x55555555.... 
	r2a = _mm_and_si128(r2a, r2c);	// Combine with counting mask, which has been combined with bit-clearing mask of 0x55555555.... 

	// Combine 2 x 128 bits and continue counting

	r1a = _mm_add_epi32(r1a, r2a);		// Combine 2 x 128 bits by adding them together

	mcl = _mm_set1_epi32(0x33333333);	// Set bit-clearing mask to 0x33333333....(alternate 2-bits)

	r1b = _mm_srli_epi32(r1a, 2);		// Shift intermediate result to right by 2 bit
	r1a = _mm_and_si128(r1a, mcl);		// Clear alternate 2-bits of intermediate result by combining with bit-clearing mask (now = 0x33333333....)
	r1b = _mm_and_si128(r1b, mcl);		// Clear alternate 2-bits of shifted intermediate result by combining with bit-clearing mask (now = 0x33333333....)
	r1a = _mm_add_epi32(r1a, r1b);		// Combine shifted and non-shifted intermediate results by adding them together

	mcl = _mm_set1_epi32(0x0F0F0F0F);	// Set bit-clearing mask to 0x0F0F0F0F....(alternate 4-bits)
	m0 = _mm_setzero_si128();			// Set an all-zero mask

	r1b = _mm_srli_epi32(r1a, 4);		// Shift intermediate result to right by 2 bit
	r1a = _mm_add_epi32(r1a, r1b);		// Combine shifted and non-shifted intermediate results by adding them together
	r1a = _mm_and_si128(r1a, mcl);		// Clear alternate 4-bits of intermediate result by combining with bit-clearing mask (now = 0xOFOFOFOF....)

	r1a = _mm_sad_epu8(r1a, m0);		// Treating the 128 bit as 16 x 8 bit; summing up the 1st 8 x 8 bit into 1st 64-bit and 2nd 8 x 8 bit into 2nd 64-bit

	return _mm_extract_epi16(r1a, 0) + _mm_extract_epi16(r1a, 4);	// Extract and return result from register

}

// Ordering of index1 and index2 is not important; this module will handle the ordering
// index1 and index2 can be on the same aligned 128 bit region or can be on adjacant aligned 128 bit region
// If index1 and index2 are in the same aligned 128 bit region, one of them must be on the boundary
// These requirements are to reduce the no. of branches in the program flow

void BWTDecodeAll(const BWT *bwt, const unsigned int index1, const unsigned int index2, unsigned int* __restrict occValue) {

	unsigned int numChar1, numChar2, minIndex, maxIndex, minIndex128, maxIndex128;
	unsigned int r;

	const static unsigned int ALIGN_16 partitionOne1[4]  = { 47, 31, 15, 0 };
	const static unsigned int ALIGN_16 partitionOne2[4]  = { 0, 15, 31, 47 };
	const static unsigned int ALIGN_16 partitionZero1[4]  = { 63, 47, 31, 15 };
	const static unsigned int ALIGN_16 partitionZero2[4]  = { 15, 31, 47, 63 };

	// SSE registers
	__m128i r1e, r2e;
	__m128i mcl;
	__m128i rc, rg, rt;
	__m128i ra1, ra2;
	__m128i rc1, rc2;
	__m128i rg1, rg2;
	__m128i rt1, rt2;


	// Sort index1 and index2
	r = (index1 - index2) & -(index1 < index2);
	minIndex = index2 + r;
	maxIndex = index1 - r;

	// Locate 128 bit boundary
	minIndex128 = lastAlignedBoundary(minIndex, CHAR_PER_128);
	maxIndex128 = lastAlignedBoundary(maxIndex - (maxIndex - minIndex > CHAR_PER_128), CHAR_PER_128);

	// Determine no.of characters to count
	numChar1 = maxIndex128 - minIndex;
	numChar2 = maxIndex - maxIndex128;

	// Load encoding into register here in the hope of hiding some memory latency
	r1e = _mm_load_si128((__m128i *)(bwt->bwtCode + minIndex128 / CHAR_PER_WORD));	// Load encoding into register
	r2e = _mm_load_si128((__m128i *)(bwt->bwtCode + maxIndex128 / CHAR_PER_WORD));	// Load encoding into register

	// Set character extraction masks 
	mcl = _mm_set1_epi32(0x55555555);						// Set bit-clearing mask to 0x55555555....(alternate 1-bit)

	// Set counting mask for 2 x 128 bits

	ra1 = _mm_set1_epi32(numChar1);		// Load number of characters into register
	ra2 = _mm_set1_epi32(numChar2);		// Load number of characters into register

	rc1 = _mm_load_si128((__m128i*)partitionOne1);	// Load partition into register
	rc2 = _mm_load_si128((__m128i*)partitionOne2);	// Load partition into register

	rg1 = _mm_load_si128((__m128i*)partitionZero1);	// Load partition into register
	rg2 = _mm_load_si128((__m128i*)partitionZero2);	// Load partition into register

	rc1 = _mm_cmpgt_epi32(ra1, rc1);				// Compare to generate 4x32 bit mask; the word with counting boundary is all ones
	rc2 = _mm_cmpgt_epi32(ra2, rc2);				// Compare to generate 4x32 bit mask; the word with counting boundary is all ones

	rg1 = _mm_cmpgt_epi32(ra1, rg1);				// Compare to generate 4x32 bit mask; the word with counting boundary is all zeros
	rg2 = _mm_cmpgt_epi32(ra2, rg2);				// Compare to generate 4x32 bit mask; the word with counting boundary is all zeros

	rc1 = _mm_srli_epi32(rc1, (16 - numChar1 % 16) * 2);	// Shift bits so that all word comform to the requirement of counting the word with counting boundary 
	rc2 = _mm_slli_epi32(rc2, (16 - numChar2 % 16) * 2);	// Shift bits so that all word comform to the requirement of counting the word with counting boundary

	ra1 = _mm_or_si128(rc1, rg1);	// Combine two masks
	ra2 = _mm_or_si128(rc2, rg2);	// Combine two masks

	// Start counting; encoding has been loaded into register earlier
	r1e = _mm_and_si128(r1e, ra1);	// Combine encoding with counting mask
	r2e = _mm_and_si128(r2e, ra2);	// Combine encoding with counting mask

	// ra1, ra2, rc1, rc2, rg1, rg2, rt1, rt2 all retired

	// Shift and combine with character selection mask

	ra1 = _mm_srli_epi32(r1e, 1);	// Shift encoding to right by 1 bit
	ra2 = _mm_srli_epi32(r2e, 1);	// Shift encoding to right by 1 bit

	rt1 = _mm_and_si128(r1e, mcl);	// Check even-bits = '1'
	rt2 = _mm_and_si128(r2e, mcl);	// Check even-bits = '1'

	rg1 = _mm_and_si128(ra1, mcl);	// Check odd-bits = '1'
	rg2 = _mm_and_si128(ra2, mcl);	// Check odd-bits = '1'

	rc1 = _mm_andnot_si128(r1e, mcl);	// Check even-bits = '0'
	rc2 = _mm_andnot_si128(r2e, mcl);	// Check even-bits = '0'

	ra1 = _mm_andnot_si128(ra1, mcl);	// Check odd-bits = '0'
	ra2 = _mm_andnot_si128(ra2, mcl);	// Check odd-bits = '0'

	// r1e, r2e retired

	// Count for 'c' 'g' 't'

	r1e = _mm_and_si128(ra1, rt1);		// Combine even and odd bits
	r2e = _mm_and_si128(ra2, rt2);		// Combine even and odd bits
	ra1 = _mm_and_si128(rg1, rc1);		// Combine even and odd bits
	ra2 = _mm_and_si128(rg2, rc2);		// Combine even and odd bits
	rc1 = _mm_and_si128(rg1, rt1);		// Combine even and odd bits
	rc2 = _mm_and_si128(rg2, rt2);		// Combine even and odd bits

	rc = _mm_add_epi32(r1e, r2e);		// Combine 2 x 128 bits by adding them together
	rg = _mm_add_epi32(ra1, ra2);		// Combine 2 x 128 bits by adding them together
	rt = _mm_add_epi32(rc1, rc2);		// Combine 2 x 128 bits by adding them together

	// All except rc, rg, rt retired

	// Continue counting rc, rg, rt

	mcl = _mm_set1_epi32(0x33333333);	// Set bit-clearing mask to 0x33333333....(alternate 2-bits)

	rc1 = _mm_srli_epi32(rc, 2);		// Shift intermediate result to right by 2 bit
	rg1 = _mm_srli_epi32(rg, 2);		// Shift intermediate result to right by 2 bit
	rt1 = _mm_srli_epi32(rt, 2);		// Shift intermediate result to right by 2 bit

	rc2 = _mm_and_si128(rc, mcl);		// Clear alternate 2-bits of intermediate result by combining with bit-clearing mask (now = 0x33333333....)
	rg2 = _mm_and_si128(rg, mcl);		// Clear alternate 2-bits of intermediate result by combining with bit-clearing mask (now = 0x33333333....)
	rt2 = _mm_and_si128(rt, mcl);		// Clear alternate 2-bits of intermediate result by combining with bit-clearing mask (now = 0x33333333....)

	rc1 = _mm_and_si128(rc1, mcl);		// Clear alternate 2-bits of shifted intermediate result by combining with bit-clearing mask (now = 0x33333333....)
	rg1 = _mm_and_si128(rg1, mcl);		// Clear alternate 2-bits of shifted intermediate result by combining with bit-clearing mask (now = 0x33333333....)
	rt1 = _mm_and_si128(rt1, mcl);		// Clear alternate 2-bits of shifted intermediate result by combining with bit-clearing mask (now = 0x33333333....)

	rc = _mm_add_epi32(rc1, rc2);		// Combine shifted and non-shifted intermediate results by adding them together
	rg = _mm_add_epi32(rg1, rg2);		// Combine shifted and non-shifted intermediate results by adding them together
	rt = _mm_add_epi32(rt1, rt2);		// Combine shifted and non-shifted intermediate results by adding them together

	mcl = _mm_set1_epi32(0x0F0F0F0F);	// Set bit-clearing mask to 0x0F0F0F0F....(alternate 4-bits)
	r1e = _mm_setzero_si128();			// Set an all-zero mask

	rc1 = _mm_srli_epi32(rc, 4);		// Shift intermediate result to right by 2 bit
	rg1 = _mm_srli_epi32(rg, 4);		// Shift intermediate result to right by 2 bit
	rt1 = _mm_srli_epi32(rt, 4);		// Shift intermediate result to right by 2 bit

	rc2 = _mm_add_epi32(rc, rc1);		// Combine shifted and non-shifted intermediate results by adding them together
	rg2 = _mm_add_epi32(rg, rg1);		// Combine shifted and non-shifted intermediate results by adding them together
	rt2 = _mm_add_epi32(rt, rt1);		// Combine shifted and non-shifted intermediate results by adding them together

	rc = _mm_and_si128(rc2, mcl);		// Clear alternate 4-bits of intermediate result by combining with bit-clearing mask (now = 0xOFOFOFOF....)
	rg = _mm_and_si128(rg2, mcl);		// Clear alternate 4-bits of intermediate result by combining with bit-clearing mask (now = 0xOFOFOFOF....)
	rt = _mm_and_si128(rt2, mcl);		// Clear alternate 4-bits of intermediate result by combining with bit-clearing mask (now = 0xOFOFOFOF....)

	rc = _mm_sad_epu8(rc, r1e);			// Treating the 128 bit as 16 x 8 bit; summing up the 1st 8 x 8 bit into 1st 64-bit and 2nd 8 x 8 bit into 2nd 64-bit
	rg = _mm_sad_epu8(rg, r1e);			// Treating the 128 bit as 16 x 8 bit; summing up the 1st 8 x 8 bit into 1st 64-bit and 2nd 8 x 8 bit into 2nd 64-bit
	rt = _mm_sad_epu8(rt, r1e);			// Treating the 128 bit as 16 x 8 bit; summing up the 1st 8 x 8 bit into 1st 64-bit and 2nd 8 x 8 bit into 2nd 64-bit

	occValue[1] = _mm_extract_epi16(rc, 0) + _mm_extract_epi16(rc, 4);	// Extract result from register and store into variable
	occValue[2] = _mm_extract_epi16(rg, 0) + _mm_extract_epi16(rg, 4);	// Extract result from register and store into variable
	occValue[3] = _mm_extract_epi16(rt, 0) + _mm_extract_epi16(rt, 4);	// Extract result from register and store into variable
	occValue[0] = maxIndex - minIndex - occValue[1] - occValue[2] - occValue[3];

}


unsigned int BWTOccValue(const BWT *bwt, unsigned int index, const unsigned int character) {

	unsigned int occValue, decodeValue;
	unsigned int occExplicitIndex, occIndex;
	unsigned int r;

	// $ is supposed to be positioned at inverseSa0 but it is not encoded
	// therefore index is subtracted by 1 for adjustment
	index -= (index > bwt->inverseSa0);

#ifdef DEBUG
	if (index > bwt->textLength) {
		fprintf(stderr, "BWTOccValue() : index > textLength!\n");
		exit(1);
	}
#endif

	occExplicitIndex = (index + OCC_INTERVAL / 2 - 1) / OCC_INTERVAL;	// Bidirectional encoding
	occIndex = occExplicitIndex * OCC_INTERVAL;

	//_mm_prefetch((char*)(memory + address[j+1]), _MM_HINT_NTA);


	occValue = BWTOccValueExplicit(bwt, occExplicitIndex, character);
#ifdef DEBUG
	if (occValue > occIndex) {
		fprintf(stderr, "BWTOccValue() : occValueExplicit > occIndex!\n");
		exit(1);
	}
#endif

	if (occIndex != index) {
		decodeValue = BWTDecode(bwt, occIndex, index, character);
		r = -(occIndex > index);
		return occValue + (decodeValue & ~r) - (decodeValue & r);
	} else {
		return occValue;
	}

}

void BWTOccValueTwoIndex(const BWT *bwt, unsigned int index1, unsigned int index2, const unsigned int character, unsigned int* __restrict occValue) {

	unsigned int decodeValue, tempExplicit1, tempExplicit2, tempOccValue1, tempOccValue2;
	unsigned int occExplicitIndex1, occIndex1;
	unsigned int occExplicitIndex2, occIndex2;
	unsigned int r;

	// $ is supposed to be positioned at inverseSa0 but it is not encoded
	// therefore index is subtracted by 1 for adjustment
	index1 -= (index1 > bwt->inverseSa0);
	index2 -= (index2 > bwt->inverseSa0);

#ifdef DEBUG
	if (index1 > bwt->textLength) {
		fprintf(stderr, "BWTOccValueTwoIndex() : index1 > textLength!\n");
		exit(1);
	}
	if (index2 > bwt->textLength) {
		fprintf(stderr, "BWTOccValueTwoIndex() : index2 > textLength!\n");
		exit(1);
	}
#endif

	// Pre-fetch memory to be accessed
	BWTPrefetchBWT(bwt, index1);
	BWTPrefetchBWT(bwt, index2);

	occExplicitIndex1 = (index1 + OCC_INTERVAL / 2 - 1) / OCC_INTERVAL;	// Bidirectional encoding
	occIndex1 = occExplicitIndex1 * OCC_INTERVAL;
	occExplicitIndex2 = (index2 + OCC_INTERVAL / 2 - 1) / OCC_INTERVAL;	// Bidirectional encoding
	occIndex2 = occExplicitIndex2 * OCC_INTERVAL;

	// Pre-fetch memory to be accessed
	BWTPrefetchOccValueExplicit(bwt, occExplicitIndex1);
	BWTPrefetchOccValueExplicit(bwt, occExplicitIndex2);


	if (occIndex1 != index1) {
		decodeValue = BWTDecode(bwt, occIndex1, index1, character);
		r = -(occIndex1 > index1);
		tempOccValue1 = (decodeValue & ~r) - (decodeValue & r);
	} else {
		tempOccValue1 = 0;
	}

	if (occIndex2 != index2) {
		decodeValue = BWTDecode(bwt, occIndex2, index2, character);
		r = -(occIndex2 > index2);
		tempOccValue2 = (decodeValue & ~r) - (decodeValue & r);
	} else {
		tempOccValue2 = 0;
	}

	tempExplicit1 = BWTOccValueExplicit(bwt, occExplicitIndex1, character);
	tempExplicit2 = BWTOccValueExplicit(bwt, occExplicitIndex2, character);
#ifdef DEBUG
	if (tempExplicit1 > occIndex1) {
		fprintf(stderr, "BWTOccValueTwoIndex() : occValueExplicit1 > occIndex1!\n");
		exit(1);
	}
	if (tempExplicit2 > occIndex2) {
		fprintf(stderr, "BWTOccValueTwoIndex() : occValueExplicit2 > occIndex2!\n");
		exit(1);
	}
#endif

	occValue[0] = tempOccValue1 + tempExplicit1;
	occValue[1] = tempOccValue2 + tempExplicit2;

}


void BWTAllOccValue(const BWT *bwt, unsigned int index, unsigned int* __restrict occValue) {

	unsigned int occExplicitIndex, occIndex;
	unsigned int ALIGN_16 tempOccValue[ALPHABET_SIZE];
	unsigned int r;

	// SSE registers
	__m128i rtov, rov, rc, t1, t2;

	// $ is supposed to be positioned at inverseSa0 but it is not encoded
	// therefore index is subtracted by 1 for adjustment
	index -= (index > bwt->inverseSa0);

#ifdef DEBUG
	if (index > bwt->textLength) {
		fprintf(stderr, "BWTOccValue() : index > textLength!\n");
		exit(1);
	}
#endif

	occExplicitIndex = (index + OCC_INTERVAL / 2 - 1) / OCC_INTERVAL;	// Bidirectional encoding
	occIndex = occExplicitIndex * OCC_INTERVAL;

	BWTAllOccValueExplicit(bwt, occExplicitIndex, occValue);

	if (occIndex != index) {

		BWTDecodeAll(bwt, occIndex, index, tempOccValue);

		// The following code add tempOccvalue to occValue if index > occIndex and subtract tempOccValue from occValue if occIndex > index
		r = -(occIndex > index);
		rc = _mm_set1_epi32(r);				// Set rc = r r r r
		rtov = _mm_load_si128((__m128i*)tempOccValue);
		rov = _mm_load_si128((__m128i*)occValue);
		t1 = _mm_andnot_si128(rc, rtov);
		t2 = _mm_and_si128(rc, rtov);
		rov = _mm_add_epi32(rov, t1);
		rov = _mm_sub_epi32(rov, t2);
		_mm_store_si128((__m128i*)occValue, rov);

	} else {
		return;
	}

}

void BWTAllOccValueTwoIndex(const BWT *bwt, unsigned int index1, unsigned int index2, unsigned int* __restrict occValue1, unsigned int* __restrict occValue2) {

	unsigned int occExplicitIndex1, occIndex1;
	unsigned int occExplicitIndex2, occIndex2;
	unsigned int ALIGN_16 tempOccValue1[ALPHABET_SIZE];
	unsigned int ALIGN_16 tempOccValue2[ALPHABET_SIZE];
	unsigned int r;

	// SSE registers
	__m128i rtov, rc, t1, t2, o1, o2;

	// $ is supposed to be positioned at inverseSa0 but it is not encoded
	// therefore index is subtracted by 1 for adjustment
	index1 -= (index1 > bwt->inverseSa0);
	index2 -= (index2 > bwt->inverseSa0);

#ifdef DEBUG
	if (index1 > index2) {
		fprintf(stderr, "BWTAllOccValueTwoIndex() : index1 > index2!\n");
		exit(1);
	}
	if (index2 > bwt->textLength) {
		fprintf(stderr, "BWTAllOccValueTwoIndex() : index2 > textLength!\n");
		exit(1);
	}
#endif

	// Pre-fetch memory to be accessed
	BWTPrefetchBWT(bwt, index1);
	BWTPrefetchBWT(bwt, index2);

	occExplicitIndex1 = (index1 + OCC_INTERVAL / 2 - 1) / OCC_INTERVAL;	// Bidirectional encoding
	occIndex1 = occExplicitIndex1 * OCC_INTERVAL;
	occExplicitIndex2 = (index2 + OCC_INTERVAL / 2 - 1) / OCC_INTERVAL;	// Bidirectional encoding
	occIndex2 = occExplicitIndex2 * OCC_INTERVAL;

	// Pre-fetch memory to be accessed
	BWTPrefetchOccValueExplicit(bwt, occExplicitIndex1);
	BWTPrefetchOccValueExplicit(bwt, occExplicitIndex2);

	if (occIndex1 != index1) {

		BWTDecodeAll(bwt, occIndex1, index1, tempOccValue1);

		// The following code add tempOccvalue to occValue if index > occIndex and subtract tempOccValue from occValue if occIndex > index
		r = -(occIndex1 > index1);
		rtov = _mm_load_si128((__m128i*)tempOccValue1);
		rc = _mm_set1_epi32(r);				// Set rc = r r r r
		t1 = _mm_andnot_si128(rc, rtov);
		t2 = _mm_and_si128(rc, rtov);
		o1 = _mm_sub_epi32(t1, t2);
	} else {
		o1 = _mm_setzero_si128();
	}
/*
	if (occIndex1 != index1) {
		if (occIndex1 < index1) {
			ForwardDNAAllOccCount(bwt->bwtCode + occIndex1 / CHAR_PER_WORD, index1 - occIndex1, tempOccValue, bwt->decodeTable);
			occValue1[0] += tempOccValue[0];
			occValue1[1] += tempOccValue[1];
			occValue1[2] += tempOccValue[2];
			occValue1[3] += tempOccValue[3];
		} else {
			BackwardDNAAllOccCount(bwt->bwtCode + occIndex1 / CHAR_PER_WORD, occIndex1 - index1, tempOccValue, bwt->decodeTable);
			occValue1[0] -= tempOccValue[0];
			occValue1[1] -= tempOccValue[1];
			occValue1[2] -= tempOccValue[2];
			occValue1[3] -= tempOccValue[3];
		}
	}
*/
	if (occIndex2 != index2) {

		BWTDecodeAll(bwt, occIndex2, index2, tempOccValue2);

		// The following code add tempOccvalue to occValue if index > occIndex and subtract tempOccValue from occValue if occIndex > index
		r = -(occIndex1 > index2);
		rc = _mm_set1_epi32(r);				// Set rc = r r r r
		rtov = _mm_load_si128((__m128i*)tempOccValue2);
		t1 = _mm_andnot_si128(rc, rtov);
		t2 = _mm_and_si128(rc, rtov);
		o2 = _mm_sub_epi32(t1, t2);

	} else {
		o2 = _mm_setzero_si128();
	}

	BWTAllOccValueExplicit(bwt, occExplicitIndex1, occValue1);
	BWTAllOccValueExplicit(bwt, occExplicitIndex2, occValue2);

	t1 = _mm_load_si128((__m128i*)occValue1);
	t2 = _mm_load_si128((__m128i*)occValue2);

	t1 = _mm_add_epi32(t1, o1);
	t2 = _mm_add_epi32(t2, o2);

	_mm_store_si128((__m128i*)occValue1, t1);
	_mm_store_si128((__m128i*)occValue2, t2);



/*
	if (occIndex2 != index2) {
		if (occIndex2 < index2) {
			ForwardDNAAllOccCount(bwt->bwtCode + occIndex2 / CHAR_PER_WORD, index2 - occIndex2, tempOccValue, bwt->decodeTable);
			occValue2[0] += tempOccValue[0];
			occValue2[1] += tempOccValue[1];
			occValue2[2] += tempOccValue[2];
			occValue2[3] += tempOccValue[3];
		} else {
			BackwardDNAAllOccCount(bwt->bwtCode + occIndex2 / CHAR_PER_WORD, occIndex2 - index2, tempOccValue, bwt->decodeTable);
			occValue2[0] -= tempOccValue[0];
			occValue2[1] -= tempOccValue[1];
			occValue2[2] -= tempOccValue[2];
			occValue2[3] -= tempOccValue[3];
		}
	}
*/

}

unsigned int BWTOccValueOnSpot(const BWT *bwt, unsigned int index, unsigned int* __restrict character) {

	unsigned int occExplicitIndex, occIndex;
	unsigned int occValue, decodeValue;
	unsigned int r;

	// The bwt character before index will be returned and the count will be up to that bwt character
	#ifdef DEBUG
	if (index == bwt->inverseSa0 + 1) {
		fprintf(stderr, "BWTOccValueOnSpot(): index = inverseSa0 + 1!\n");
		exit(1);
	}
	if (index > bwt->textLength + 1) {
		fprintf(stderr, "BWTOccValueOnSpot() : index > textLength!\n");
		exit(1);
	}
	if (index == 0) {
		fprintf(stderr, "BWTOccValueOnSpot() : index = 0!\n");
		exit(1);
	}
	#endif

	// $ is supposed to be positioned at inverseSa0 but it is not encoded
	// therefore index is incremented for adjustment
	index -= (index > bwt->inverseSa0);

	// Bidirectional encoding
	occExplicitIndex = (index + OCC_INTERVAL / 2 - 1) / OCC_INTERVAL;
	occIndex = occExplicitIndex * OCC_INTERVAL;

	*character = bwt->bwtCode[(index - 1) / CHAR_PER_WORD] << (((index - 1) % CHAR_PER_WORD) * BIT_PER_CHAR) >> (BITS_IN_WORD - BIT_PER_CHAR);
	occValue = BWTOccValueExplicit(bwt, occExplicitIndex, *character);

	if (occIndex != index) {
		decodeValue = BWTDecode(bwt, occIndex, index, *character);
		r = -(occIndex > index);
		return occValue + (decodeValue & ~r) - (decodeValue & r);
	} else {
		return occValue;
	}

}

unsigned int BWTSearchOccValue(const BWT *bwt, const unsigned int character, const unsigned int searchOccValue) {

	unsigned int occValue;
	unsigned int i,j;
	unsigned int c;
	unsigned int bwtPos;
	unsigned int occExplicitIndexLeft, occExplicitIndexRight, occExplicitIndexMiddle;

	#ifdef DEBUG
	if (searchOccValue == 0 || searchOccValue > bwt->textLength) {
		fprintf(stderr, "BWTSearchOccValue() : searchOccValue out of bound!\n");
		exit(1);
	}
	#endif

	// Search Occurrence value

	occExplicitIndexLeft = 0;
	occExplicitIndexRight = (bwt->textLength + OCC_INTERVAL - 1) / OCC_INTERVAL;

	while (occExplicitIndexLeft + 1 < occExplicitIndexRight) {
		occExplicitIndexMiddle = average(occExplicitIndexLeft, occExplicitIndexRight);
		if (searchOccValue > BWTOccValueExplicit(bwt, occExplicitIndexMiddle, character)) {
			occExplicitIndexLeft = occExplicitIndexMiddle;
		} else {
			occExplicitIndexRight = occExplicitIndexMiddle;
		}
	}

	// Not tuned for DNA
	occValue = BWTOccValueExplicit(bwt, occExplicitIndexLeft, character);
	bwtPos = occExplicitIndexLeft * OCC_INTERVAL / CHAR_PER_WORD;

	for (i=0; i < OCC_INTERVAL / CHAR_PER_WORD; i++) {
		c = bwt->bwtCode[bwtPos + i];
		for (j=0; j < CHAR_PER_WORD && occValue < searchOccValue; j++) {
			if (c >> (BITS_IN_WORD - BIT_PER_CHAR) == character) {
				occValue++;
				if (occValue >= searchOccValue) {
					return occExplicitIndexLeft * OCC_INTERVAL + i * CHAR_PER_WORD + j;
				}
			}
			c <<= BIT_PER_CHAR;
		}
	}

	fprintf(stderr, "BWTSearchOccValue() : unexpected error!\n");
	exit(1);

}

static INLINE unsigned int BWTOccValueExplicit(const BWT *bwt, const unsigned int occIndexExplicit, const unsigned int character) {

	unsigned int occIndexMajor;
	unsigned int compareMask, shift, mask;

	occIndexMajor = occIndexExplicit * OCC_INTERVAL / OCC_INTERVAL_MAJOR;

	compareMask = (-(occIndexExplicit % OCC_VALUE_PER_WORD == 0));
	shift = 16 & compareMask;
	mask = 0x0000FFFF | compareMask;

	return bwt->occValueMajor[occIndexMajor * ALPHABET_SIZE + character] +
			((bwt->occValue[occIndexExplicit / OCC_VALUE_PER_WORD * ALPHABET_SIZE + character] >> shift) & mask);

}

static INLINE void BWTAllOccValueExplicit(const BWT *bwt, const unsigned int occIndexExplicit, unsigned int* __restrict occValueExplicit) {

	unsigned int occIndexMajor;
	unsigned int compareMask, shift, mask;

	__m128i v1, v2, m;

	occIndexMajor = occIndexExplicit * OCC_INTERVAL / OCC_INTERVAL_MAJOR;

	compareMask = (-(occIndexExplicit % OCC_VALUE_PER_WORD == 0));
	shift = 16 & compareMask;
	mask = 0x0000FFFF | compareMask;

	v2 = _mm_load_si128((__m128i *)(bwt->occValue + occIndexExplicit / OCC_VALUE_PER_WORD * ALPHABET_SIZE));
	v1 = _mm_load_si128((__m128i *)(bwt->occValueMajor + occIndexMajor * ALPHABET_SIZE));

	m = _mm_set1_epi32(mask);

	v2 = _mm_srli_epi32(v2, shift);
	v2 = _mm_and_si128(v2, m);

	v1 = _mm_add_epi32(v1, v2);

	_mm_store_si128((__m128i*)occValueExplicit, v1);

}

static INLINE void BWTPrefetchOccValueExplicit(const BWT *bwt, const unsigned int occIndexExplicit) {

	unsigned int occIndexMajor;

	occIndexMajor = occIndexExplicit * OCC_INTERVAL / OCC_INTERVAL_MAJOR;

	_mm_prefetch((char*)(bwt->occValueMajor + occIndexMajor * ALPHABET_SIZE), _MM_HINT_T0);
	_mm_prefetch((char*)(bwt->occValue + occIndexExplicit / OCC_VALUE_PER_WORD * ALPHABET_SIZE), _MM_HINT_NTA);

}

static INLINE void BWTPrefetchBWT(const BWT *bwt, const unsigned int index) {

	_mm_prefetch((char*)(bwt->bwtCode + index / CHAR_PER_WORD), _MM_HINT_NTA);

}


unsigned int BWTResidentSizeInWord(const unsigned int numChar) {

	unsigned int numCharRoundUpToOccInterval;

	// The $ in BWT at the position of inverseSa0 is not encoded
	numCharRoundUpToOccInterval = (numChar + OCC_INTERVAL - 1) / OCC_INTERVAL * OCC_INTERVAL;

	return (numCharRoundUpToOccInterval + CHAR_PER_WORD - 1) / CHAR_PER_WORD;

}

unsigned int BWTFileSizeInWord(const unsigned int numChar) {

	// The $ in BWT at the position of inverseSa0 is not encoded
	return (numChar + CHAR_PER_WORD - 1) / CHAR_PER_WORD;

}

unsigned int BWTOccValueMinorSizeInWord(const unsigned int numChar) {

	unsigned int numOfOccValue;

	numOfOccValue = (numChar + OCC_INTERVAL - 1) / OCC_INTERVAL + 1;		// Value at both end for bi-directional encoding
	return (numOfOccValue + OCC_VALUE_PER_WORD - 1) / OCC_VALUE_PER_WORD * ALPHABET_SIZE;

}

unsigned int BWTOccValueMajorSizeInWord(const unsigned int numChar) {

	unsigned int numOfOccValue;
	unsigned int numOfOccIntervalPerMajor;

	numOfOccValue = (numChar + OCC_INTERVAL - 1) / OCC_INTERVAL + 1;				// Value at both end for bi-directional encoding
	numOfOccIntervalPerMajor = OCC_INTERVAL_MAJOR / OCC_INTERVAL;

	return (numOfOccValue + numOfOccIntervalPerMajor - 1) / numOfOccIntervalPerMajor * ALPHABET_SIZE;

}

void BWTClearTrailingBwtCode(BWT *bwt) {

	unsigned int bwtResidentSizeInWord;
	unsigned int wordIndex, offset;
	unsigned int i;

	bwtResidentSizeInWord = BWTResidentSizeInWord(bwt->textLength);

	wordIndex = bwt->textLength / CHAR_PER_WORD;
	offset = (bwt->textLength - wordIndex * CHAR_PER_WORD) * BIT_PER_CHAR;
	if (offset > 0) {
		bwt->bwtCode[wordIndex] = truncateRight(bwt->bwtCode[wordIndex], BITS_IN_WORD - offset);
	} else {
		if (wordIndex < bwtResidentSizeInWord) {
			bwt->bwtCode[wordIndex] = 0;
		}
	}

	for (i=wordIndex+1; i<bwtResidentSizeInWord; i++) {
		bwt->bwtCode[i] = 0;
	}

}

unsigned int BWTPsiMinusValue(const BWT *bwt, const unsigned int index) {

	unsigned int c;
	unsigned int occValue;

	#ifdef DEBUG
	if (index > bwt->textLength) {
		fprintf(stderr, "BWTPsiMinusValue() : index out of range!\n");
		exit(1);
	}
	#endif

	if (index != bwt->inverseSa0) {

		occValue = BWTOccValueOnSpot(bwt, index + 1, &c);
		occValue += bwt->cumulativeFreq[c];

		return occValue;

	} else {
		return 0;
	}

}

unsigned int BWTPsiPlusValue(const BWT *bwt, const unsigned int index) {

	unsigned int c;
	unsigned int psiPlusValue;

	#ifdef DEBUG
	if (index > bwt->textLength) {
		fprintf(stderr, "BWTPsiPlusValue() : index out of range!\n");
		exit(1);
	}
	#endif

	if (index == 0) {
		return bwt->inverseSa0;
	}

	// Find the BWT of PSI+
	c = (index > bwt->cumulativeFreq[1]) + (index > bwt->cumulativeFreq[2])
										 + (index > bwt->cumulativeFreq[3]);

	psiPlusValue = BWTSearchOccValue(bwt, c, index - bwt->cumulativeFreq[c]);
	if (psiPlusValue >= bwt->inverseSa0) {
		psiPlusValue++;
	}
	return psiPlusValue;

}

unsigned int BWTSaValue(const BWT *bwt, unsigned int saIndex) {

	unsigned int saValueSkipped = 0;

	#ifdef DEBUG
	if (saIndex > bwt->textLength) {
		fprintf(stderr, "BWTSaValue() : Index out of range!\n");
		exit(1);
	}
	if (bwt->saValue == NULL) {
		fprintf(stderr, "BWTSaValue() : Explicit SA value is not loaded!\n");
		exit(1);
	}
	#endif

	while (saIndex % bwt->saInterval != 0) {
		saValueSkipped++;
		saIndex = BWTPsiMinusValue(bwt, saIndex);
	}

	#ifdef DEBUG
	if (bwt->saValue[saIndex/bwt->saInterval] + saValueSkipped > bwt->textLength) {
		fprintf(stderr, "BWTSaValue() : saValue out of range!\n");
		exit(1);
	}
	#endif

	// SA[0] stores -1 although it should be textLength
	// PsiMinusValue returns 0 on inverseSa0
	return bwt->saValue[saIndex/bwt->saInterval] + saValueSkipped;

}

unsigned int BWTInverseSa(const BWT *bwt, unsigned int saValue) {

	unsigned int i;
	unsigned int saIndex;
	unsigned int inverseSaExplicitIndex;
	unsigned int saValueToSkip;

	#ifdef DEBUG
	if (saValue > bwt->textLength) {
		fprintf(stderr, "BWTInverseSa() : Index out of range!\n");
		exit(1);
	}
	if (bwt->inverseSa == NULL) {
		fprintf(stderr, "BWTInverseSa() : Explicit inverse SA is not loaded!\n");
		exit(1);
	}
	#endif

	inverseSaExplicitIndex = (saValue + bwt->inverseSaInterval - 1) / bwt->inverseSaInterval;
	if (inverseSaExplicitIndex * bwt->inverseSaInterval > bwt->textLength) {
		saIndex = 0;
		saValueToSkip = bwt->textLength - saValue;
	} else {
		saIndex = bwt->inverseSa[inverseSaExplicitIndex];
		saValueToSkip = inverseSaExplicitIndex * bwt->inverseSaInterval - saValue;
	}

	for (i=0; i<saValueToSkip; i++) {
		saIndex = BWTPsiMinusValue(bwt, saIndex);
	}

	return saIndex;

}

static INLINE unsigned int BWTGetWordPackedText(const unsigned int *packedText, const unsigned int index, const unsigned int shift, const unsigned int numOfBit) {

	unsigned int text;
	const static unsigned int mask[32] = { 0x00000000, 0x80000000, 0xC0000000, 0xE0000000,
								  0xF0000000, 0xF8000000, 0xFC000000, 0xFE000000,
								  0xFF000000, 0xFF800000, 0xFFC00000, 0xFFE00000,
								  0xFFF00000, 0xFFF80000, 0xFFFC0000, 0xFFFE0000,
								  0xFFFF0000, 0xFFFF8000, 0xFFFFC000, 0xFFFFE000,
								  0xFFFFF000, 0xFFFFF800, 0xFFFFFC00, 0xFFFFFE00,
								  0xFFFFFF00, 0xFFFFFF80, 0xFFFFFFC0, 0xFFFFFFE0,
								  0xFFFFFFF0, 0xFFFFFFF8, 0xFFFFFFFC, 0xFFFFFFFE };

	if (shift > 0) {
		// packedText should be allocated with at least 1 Word buffer initialized to zero
		text = (packedText[index] << shift) | (packedText[index + 1] >> (BITS_IN_WORD - shift));
	} else {
		text = packedText[index];
	}

	if (numOfBit < BITS_IN_WORD) {
		// Fill unused bit with zero
		text &= mask[numOfBit];
	}

	return text;
}

int BWTForwardSearch(const unsigned int *packedKey, const unsigned int keyLength, const BWT *bwt, const unsigned int *packedText) {

	unsigned int startSaIndex, endSaIndex, saIndexMiddle;
	unsigned int saExplicitIndexLeft, saExplicitIndexRight, saExplicitIndexMiddle;
	unsigned int saValue;

	unsigned int firstChar;
	unsigned int index, shift;
	unsigned int packedKeyLength, keyLengthInBit;
	unsigned int llcp, rlcp, mlcp, maxlcp;
	unsigned int p = 0;	// to avoid compiler warning only

	if (keyLength % CHAR_PER_WORD == 0) {
		packedKeyLength = keyLength / CHAR_PER_WORD;
		keyLengthInBit = packedKeyLength * BITS_IN_WORD;
	} else {
		packedKeyLength = keyLength / CHAR_PER_WORD + 1;
		keyLengthInBit = (keyLength / CHAR_PER_WORD) * BITS_IN_WORD + 
						 (keyLength % CHAR_PER_WORD) * BIT_PER_CHAR;
	}

	// Get the SA index initial range by retrieving cumulative frequency
	firstChar = packedKey[0] >> (BITS_IN_WORD - BIT_PER_CHAR);

	startSaIndex = bwt->cumulativeFreq[firstChar] + 1;
	endSaIndex = bwt->cumulativeFreq[firstChar + 1];

	if (startSaIndex > endSaIndex) {
		// The first character of search pattern does not exists in text
		return 0;
	}

	// Find lcp for left boundary
	saValue = bwt->saValueOnBoundary[firstChar * 2];		// Pre-calculated

	// restriction for positions near the end of text
	maxlcp = min(packedKeyLength, (bwt->textLength - saValue + CHAR_PER_WORD - 1) / CHAR_PER_WORD);

	shift = BIT_PER_CHAR * (saValue % CHAR_PER_WORD);
	index = saValue / CHAR_PER_WORD;

	llcp = 0;
	while (llcp < maxlcp && packedKey[llcp] == 
					BWTGetWordPackedText(packedText, index + llcp, shift, keyLengthInBit - llcp * BITS_IN_WORD)) {
		llcp++;
	}
	if ((saValue + keyLength > bwt->textLength) && llcp == maxlcp) {
		llcp--;
	}
	if (llcp == packedKeyLength) {
		return 1;
	}

	// Find lcp for right boundary
	saValue = bwt->saValueOnBoundary[firstChar * 2 + 1];	// Pre-calculated

	// restriction for positions near the end of text
	maxlcp = min(packedKeyLength, (bwt->textLength - saValue + CHAR_PER_WORD - 1) / CHAR_PER_WORD);

	shift = BIT_PER_CHAR * (saValue % CHAR_PER_WORD);
	index = saValue / CHAR_PER_WORD;

	rlcp = 0;
	while (rlcp < maxlcp && packedKey[rlcp] == 
					BWTGetWordPackedText(packedText, index + rlcp, shift, keyLengthInBit - rlcp * BITS_IN_WORD)) {
		rlcp++;
	}
	if ((saValue + keyLength > bwt->textLength) && rlcp == maxlcp) {
		rlcp--;
	}
	if (rlcp == packedKeyLength) {
		return 1;
	}

	// Locate in SA index explicitly stored
	saExplicitIndexLeft = startSaIndex / bwt->saInterval;
	saExplicitIndexRight = (endSaIndex - 1) / bwt->saInterval + 1;

	// loop until two adjacent SA explicit index is found
	while (saExplicitIndexLeft + 1 < saExplicitIndexRight) {

		saExplicitIndexMiddle = average(saExplicitIndexLeft, saExplicitIndexRight);

		saValue = bwt->saValue[saExplicitIndexMiddle];
		shift = BIT_PER_CHAR * (saValue % CHAR_PER_WORD);
		index = saValue / CHAR_PER_WORD;

		// Try to increase mlcp
		mlcp = min(llcp, rlcp);		// mlcp = the characters (in unit of 16 for DNA) matched so far
		// restriction for positions near the end of text
		maxlcp = min(packedKeyLength, (bwt->textLength - saValue + CHAR_PER_WORD - 1) / CHAR_PER_WORD);

		while (mlcp < maxlcp) {
			p = BWTGetWordPackedText(packedText, index + mlcp, shift, keyLengthInBit - mlcp * BITS_IN_WORD);
			if (packedKey[mlcp] != p) {
				break;
			}
			mlcp++;
		}
		if ((saValue + keyLength <= bwt->textLength) || mlcp != maxlcp) {
			if (mlcp == packedKeyLength) {
				return 1;
			}
			if (packedKey[mlcp] > p) {
				llcp = mlcp;
				saExplicitIndexLeft = saExplicitIndexMiddle;
			} else {
				rlcp = mlcp;
				saExplicitIndexRight = saExplicitIndexMiddle;
			}
		} else {
			if (packedKey[mlcp-1] >= p) {
				llcp = mlcp - 1;
				saExplicitIndexLeft = saExplicitIndexMiddle;
			} else {
				rlcp = mlcp - 1;
				saExplicitIndexRight = saExplicitIndexMiddle;
			}
			
		}

	}

	// Two adjacent SA explicit index is found, convert back to SA index
	if (saExplicitIndexLeft == startSaIndex / bwt->saInterval) {
		startSaIndex = bwt->cumulativeFreq[firstChar] + 1;
	} else {
		startSaIndex = saExplicitIndexLeft * bwt->saInterval;
	}
	if (saExplicitIndexRight == (endSaIndex - 1) / bwt->saInterval + 1) {
		endSaIndex = bwt->cumulativeFreq[firstChar + 1];
	} else {
		endSaIndex = saExplicitIndexRight * bwt->saInterval;
	}

	// binary search by decoding bwt

	while (startSaIndex < endSaIndex) {

		saIndexMiddle = average(startSaIndex, endSaIndex);

		saValue = BWTSaValue(bwt, saIndexMiddle);
		shift = BIT_PER_CHAR * (saValue % CHAR_PER_WORD);
		index = saValue / CHAR_PER_WORD;

		// Try to increase mlcp
		mlcp = min(llcp, rlcp);		// mlcp = the characters (in unit of 16 for DNA) matched so far
		// restriction for positions near the end of text
		maxlcp = min(packedKeyLength, (bwt->textLength - saValue + CHAR_PER_WORD - 1) / CHAR_PER_WORD);

		while (mlcp < maxlcp) {
			p = BWTGetWordPackedText(packedText, index + mlcp, shift, keyLengthInBit - mlcp * BITS_IN_WORD);
			if (packedKey[mlcp] != p) {
				break;
			}
			mlcp++;
		}
		if ((saValue + keyLength <= bwt->textLength) || mlcp != maxlcp) {
			if (mlcp == packedKeyLength) {
				return 1;
			}
			if (packedKey[mlcp] > p) {
				llcp = mlcp;
				startSaIndex = saIndexMiddle + 1;
			} else {
				rlcp = mlcp;
				endSaIndex = saIndexMiddle;
			}
		} else {
			if (packedKey[mlcp-1] >= p) {
				llcp = mlcp - 1;
				startSaIndex = saIndexMiddle + 1;
			} else {
				rlcp = mlcp - 1;
				endSaIndex = saIndexMiddle;
			}
			
		}

	}

	// no match found
	return 0;

}
