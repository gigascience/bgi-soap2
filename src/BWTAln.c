#include "BWTAln.h"

unsigned int REVBWTForwardSearch(const unsigned char *convertedkey, const unsigned int keylength, const BWT *rev_bwt, unsigned int *resultsaindexleft, unsigned int *resultsaindexright, unsigned int *rev_resultsaindexleft, unsigned int *rev_resultsaindexright) {

	unsigned int sacount=0;
	unsigned int rev_startsaindex, rev_endsaindex;
	unsigned int startsaindex, endsaindex;
	unsigned int pos = 1;
	int i;
	unsigned int c = convertedkey[0];
	unsigned int occcount_start[4];
	unsigned int occcount_end[4];
	unsigned int occcount[4];

	rev_startsaindex = rev_bwt->cumulativeFreq[convertedkey[0]]+1;
	rev_endsaindex = rev_bwt->cumulativeFreq[convertedkey[0]+1];
	startsaindex = rev_bwt->cumulativeFreq[convertedkey[0]]+1;
	endsaindex = rev_bwt->cumulativeFreq[convertedkey[0]+1];

	while (pos < keylength && startsaindex <= endsaindex) {
		c = convertedkey[pos];

		BWTAllOccValue(rev_bwt,rev_startsaindex,occcount_start);
		BWTAllOccValue(rev_bwt,rev_endsaindex + 1,occcount_end);

		rev_startsaindex = rev_bwt->cumulativeFreq[c] + occcount_start[c] + 1;
		rev_endsaindex = rev_bwt->cumulativeFreq[c] + occcount_end[c];

		occcount[3]=0;
		for (i=2;i>=0;i--) {
			occcount[i]=occcount[i+1]+occcount_end[i+1]-occcount_start[i+1];
		}

		endsaindex = endsaindex - occcount[c];
		startsaindex = endsaindex - (rev_endsaindex-rev_startsaindex);
		pos++;
	}

	*resultsaindexleft = startsaindex;
	*resultsaindexright = endsaindex;
	*rev_resultsaindexleft = rev_startsaindex;
	*rev_resultsaindexright = rev_endsaindex;

	sacount+=endsaindex-startsaindex+1;
	// number of occurrence = endsaindex - startsaindex + 1
	return sacount;

}


unsigned int REVBWTContForwardSearch(const unsigned char *convertedkey, const unsigned int start, const unsigned int len, const BWT *rev_bwt, unsigned int *sal, unsigned int *sar, unsigned int *rev_sal, unsigned int *rev_sar) {

	unsigned int sacount=0;
	unsigned int pos = start;
	unsigned char c;
	unsigned int occcount_start[4];
	unsigned int occcount_end[4];
	unsigned int occcount[4];
	int k;
	while (pos < start+len  && *sal <= *sar) {
		c = convertedkey[pos];

		BWTAllOccValue(rev_bwt,*rev_sal,occcount_start);
		BWTAllOccValue(rev_bwt,*rev_sar + 1,occcount_end);

		*rev_sal = rev_bwt->cumulativeFreq[c] + occcount_start[c] + 1;
		*rev_sar = rev_bwt->cumulativeFreq[c] + occcount_end[c];

		occcount[3]=0;
		for (k=2;k>=0;k--) {
			occcount[k]=occcount[k+1]+occcount_end[k+1]-occcount_start[k+1];
		}

		*sar = *sar - occcount[c];
		*sal = *sar - (*rev_sar-*rev_sal);

		pos++;
	}
	sacount+=*sar-*sal+1;
	return sacount;

}
unsigned int BWTContBackwardSearch(const unsigned char *convertedkey, const unsigned int start, const unsigned int len, const BWT *bwt, unsigned int *sal, unsigned int *sar) {

	unsigned int sacount=0;
	unsigned int pos = len;
	unsigned char c;

	if (*sal > *sar) {
		return 0;
	}

	while (pos > 0 && *sal <= *sar) {
		c = convertedkey[pos-1];
		*sal = bwt->cumulativeFreq[c] + BWTOccValue(bwt, *sal, c) + 1;
		*sar = bwt->cumulativeFreq[c] + BWTOccValue(bwt, *sar + 1, c);
		pos--;
	}
	sacount+=*sar-*sal+1;
	return sacount;
}

unsigned int BWTBackward1Error(const unsigned char *querypattern, const BWTOPT *bo, BWT *bwt, unsigned int start, unsigned int len, unsigned int pl, unsigned int pr, unsigned int info, HITTABLE *hits) {
		unsigned int mk_l=1,mk_r=0;
		unsigned int occcount_pstart[4];
		unsigned int occcount_pend[4];
		unsigned int sacount=0;
		int hitcount = 0;

		unsigned char c;
		unsigned char ec;
		int i;
		//printf("bwtbackward1error %u %u\n",pl,pr);
		//                                  v--start       v----start+len
		// querypattern = xxxxxxxxxxxxxxxxx[xxxxxxxxxxxxxx]xxxxxx
		//                                      <------- search direction
		//                               querypattern[start+len-1], querypattern[start+len-2]...querypattern[start]
		//                           for i=0 to len-1,
		//                                append querypatter[start+len-1-i]!

		for (i=0;(i<len && pl<=pr);i++) {
			//call once only proceduressssss - great
			BWTAllOccValue(bwt,pl,occcount_pstart);
			BWTAllOccValue(bwt,pr + 1,occcount_pend);

			//backward manner
			for (ec=0;ec<4;ec++) {
				if (querypattern[start+len-1-i]==ec)
					continue;

				info &= 0xffff000;
				info |= ((ec&3)<<8|((start+len-1-i)&0xff))&0xfff;
				mk_l=pl;
				mk_r=pr;

				mk_l = bwt->cumulativeFreq[ec] + occcount_pstart[ec] + 1;
				mk_r = bwt->cumulativeFreq[ec] + occcount_pend[ec];

				if (BWTContBackwardSearch(querypattern,start,len-i-1,bwt,&mk_l,&mk_r)) {
					hitcount += OCCProcess(mk_l,mk_r, bo, info, hits);
					sacount+=mk_r-mk_l+1;
				}
			}
			c = querypattern[start+len-1-i];
			pl = bwt->cumulativeFreq[c] + occcount_pstart[c] + 1;
			pr = bwt->cumulativeFreq[c] + occcount_pend[c];

		}
		return hitcount;
}


unsigned int REVBWTForward1Error(const unsigned char *querypattern, const BWTOPT *bo, BWT * bwt,BWT * rev_bwt, unsigned int start,unsigned int len,unsigned int pl,unsigned int pr,unsigned int rev_pl,unsigned int rev_pr, unsigned int info, HITTABLE *hits) {

		unsigned int mk_l=1,mk_r=0,rev_mk_l,rev_mk_r;
		unsigned int occcount_pstart[4];
		unsigned int occcount_pend[4];
		unsigned int occcountp[4];
		unsigned int sacount=0;
		int hitcount = 0;

		unsigned char c;
		unsigned char ec;
		unsigned int i;
		int k;
		const int coord = (info>>24&1)?(bo->seqLen-bo->alnLen):0;

		for (i=0;(i<len && pl<=pr);i++) {
			//call once only proceduressssss - great
			BWTAllOccValue(rev_bwt,rev_pl,occcount_pstart);
			BWTAllOccValue(rev_bwt,rev_pr + 1,occcount_pend);

			occcountp[3]=0;
			for (k=2;k>=0;k--) {
				occcountp[k]=occcountp[k+1]+occcount_pend[k+1]-occcount_pstart[k+1];
			}


			//forward manner
			for (ec=0;ec<4;ec++) {
				if (querypattern[start+i]==ec)
					continue;
				info &= 0xffff000;
				info |= ((ec&3)<<8|((start+i+coord)&0xff))&0xfff;
				mk_l=pl;
				mk_r=pr;
				rev_mk_l=rev_pl;
				rev_mk_r=rev_pr;

				unsigned int pos = i+1;

				rev_mk_l = rev_bwt->cumulativeFreq[ec] + occcount_pstart[ec] + 1;
				rev_mk_r = rev_bwt->cumulativeFreq[ec] + occcount_pend[ec];

				mk_r = mk_r - occcountp[ec];
				mk_l = mk_r - (rev_mk_r-rev_mk_l);

				if (REVBWTContForwardSearch(querypattern,start+pos,len-i-1,rev_bwt,&mk_l,&mk_r,&rev_mk_l,&rev_mk_r)) {
//					printf("%d\t", start+i);
//					printf("%d\n", ec);
					hitcount+= OCCProcess(mk_l,mk_r, bo, info, hits);
					//return mk_l, mk_r
					sacount+=mk_r-mk_l+1;
				}
			}
			c = querypattern[start+i];

			rev_pl = rev_bwt->cumulativeFreq[c] + occcount_pstart[c] + 1;
			rev_pr = rev_bwt->cumulativeFreq[c] + occcount_pend[c];

			pr = pr - occcountp[c];
			pl = pr - (rev_pr-rev_pl);
		}
		return hitcount;
}

int BWTExactMatching(const unsigned char *convertedKey, const BWTOPT *bo, int chain, BWT *bwt, LOOKUPTABLE *lookup, HITTABLE *hits){
	if(convertedKey == NULL) return 0;
	const unsigned int keyLength = bo->alnLen;
	LOOKUPTABLE lookupTable;
	lookupTable.tableSize = lookup->tableSize;
	lookupTable.table = lookup->table;
	unsigned int l, r;
	unsigned int i;
	int hitcount = 0;
//	fprintf(stdout, "BWTExactMatching\n");
	unsigned int info = (chain&1) << 24 ;
	unsigned long long packedPattern = 0;
//	printf("tablesize: %d, keyLength: %d\n", lookupTable.tableSize, keyLength);
	for (i = 0; i <lookupTable.tableSize ; i++) {
		packedPattern<<=2;
		packedPattern |= (convertedKey[keyLength-lookupTable.tableSize+i] & 3);
	}
	l = packedPattern ? lookupTable.table[packedPattern-1]+1 : 1;
	r = lookupTable.table[packedPattern];

	for (i = keyLength-lookupTable.tableSize; i > 0 && l <= r; --i) {
		unsigned char c = convertedKey[i-1];
		l = bwt->cumulativeFreq[c] + BWTOccValue(bwt, l, c) + 1;
		r = bwt->cumulativeFreq[c] + BWTOccValue(bwt, r + 1, c);
	}
	
	if (l<=r && hits->n < bo->cutoff) {
//		fprintf(stderr, "occ find\n");
		hitcount += OCCProcess(l, r, bo, info, hits);
		return hitcount;
	}
	return 0;
}

int BWT1ErrorMatching(const unsigned char * convertedKey, const BWTOPT *bo, const int chain, BWT *bwt, BWT *rev_bwt, LOOKUPTABLE *lookup, LOOKUPTABLE *rev_lookup, HITTABLE *hits) {
	if(convertedKey == NULL) return 0;
	LOOKUPTABLE lookupTable, rev_lookupTable;
	lookupTable.tableSize     = lookup->tableSize;
	rev_lookupTable.tableSize = rev_lookup->tableSize;
	lookupTable.table         = lookup->table;
	rev_lookupTable.table     = rev_lookup->table;
//	unsigned int cutoff       = bo->cutoff;
	unsigned int keyLength    = bo->alnLen;
	unsigned int forwardDepth = bo->h;
	unsigned int info         = (1<<25)|((chain&1)<<24);
	unsigned int l, r;
	unsigned int rev_l, rev_r;
	unsigned int i;
	int hitcount = 0;
	unsigned backwardDepth = keyLength - forwardDepth;
//	fprintf(stdout, "BWT1misMatching\n");
	//1. Backward Case
	//==============================================
	// look-up the last characters (backward)
	unsigned long long packedPattern = 0;
	for (i = 0; i <lookupTable.tableSize ; i++) {
		packedPattern<<=2;
		packedPattern |= (convertedKey[keyLength-lookupTable.tableSize+i] & 3);
	}
	l = packedPattern ? lookupTable.table[packedPattern-1]+1 : 1;
	r = lookupTable.table[packedPattern];

	// backward search with BWT until the forward depth section
	for (i = lookupTable.tableSize; i < backwardDepth && l <= r; ++i) {
		unsigned char c = convertedKey[keyLength-i-1];
		l = bwt->cumulativeFreq[c] + BWTOccValue(bwt, l, c) + 1;
		r = bwt->cumulativeFreq[c] + BWTOccValue(bwt, r + 1, c);
	}
	// error in the forward depth section of the query pattern
	hitcount += BWTBackward1Error(convertedKey,
			bo,
			bwt,
			0, forwardDepth,
			l,r,
			info,
			hits);
	//    fprintf(stdout, "saCount1: %u\n", saCount);
//	if(hits->n >= bo->cutoff) return hitcount;


	//2.Forward Case
	//==============================================
	unsigned int occCount_start[4];
	unsigned int occCount_end[4];
	unsigned int occCount[4];
	unsigned long long l_packedPattern = 0;
	unsigned long long r_packedPattern = 0;
	unsigned long long rev_packedPattern = 0;
	// look-up the first characters
	for (i = 0; i <rev_lookupTable.tableSize ; i++) {
		l_packedPattern<<=2;
		l_packedPattern |= (convertedKey[i]  & 3 );
	}
	r_packedPattern = l_packedPattern;
	//If the look-up tables are of different size
	l_packedPattern <<= (lookupTable.tableSize-rev_lookupTable.tableSize)*2;
	r_packedPattern <<= (lookupTable.tableSize-rev_lookupTable.tableSize)*2;
	r_packedPattern |= (1<<(lookupTable.tableSize-rev_lookupTable.tableSize)*2) - 1;

	for (i = 0; i <rev_lookupTable.tableSize ; i++) {
		rev_packedPattern<<=2;
		rev_packedPattern |= (convertedKey[rev_lookupTable.tableSize-i-1] & 3);
	}

	l = l_packedPattern ? lookupTable.table[l_packedPattern-1]+1 : 1;
	r = lookupTable.table[r_packedPattern];
	rev_l = rev_packedPattern ? rev_lookupTable.table[rev_packedPattern-1]+1 : 1;
	rev_r = rev_lookupTable.table[rev_packedPattern];

	// forward search with BWT until the end of forward depth section
	for (i = rev_lookupTable.tableSize; i < forwardDepth && l <= r; ++i) {
		unsigned char c = convertedKey[i];
		BWTAllOccValue(rev_bwt,rev_l,occCount_start);
		BWTAllOccValue(rev_bwt,rev_r + 1,occCount_end);

		rev_l = rev_bwt->cumulativeFreq[c] + occCount_start[c] + 1;
		rev_r = rev_bwt->cumulativeFreq[c] + occCount_end[c];

		occCount[3]=0;
		int k;
		for (k=2;k>=0;k--) {
			occCount[k]=occCount[k+1]+occCount_end[k+1]-occCount_start[k+1];
		}

		r = r - occCount[c];
		l = r - (rev_r-rev_l);
	}

	// error in the forward depth section of the query pattern
	hitcount+=REVBWTForward1Error(convertedKey,
			bo,
			bwt, rev_bwt,
			forwardDepth, backwardDepth,
			l, r,
			rev_l, rev_r,
			info,
			hits);
	//    fprintf(stdout, "saCount: %u\n", saCount);
	return hitcount;
}

int BWT2ErrorMatching(const unsigned char *convertedKey, const BWTOPT *bo, const int chain, BWT * bwt, BWT * rev_bwt, LOOKUPTABLE *lookup, LOOKUPTABLE *rev_lookup, HITTABLE *hits) {
	if(convertedKey == NULL) return 0;
	LOOKUPTABLE lookupTable, rev_lookupTable;
	lookupTable.tableSize     = lookup->tableSize;
	rev_lookupTable.tableSize = rev_lookup->tableSize;
	lookupTable.table         = lookup->table;
	rev_lookupTable.table     = rev_lookup->table;
	unsigned int keyLength    = bo->alnLen;
//	fprintf(stderr, "keyLength %u\n", keyLength);
	unsigned int sizeX        = bo->x;
	unsigned int sizeY        = bo->y;
	unsigned int cutoff       = bo->cutoff;
	unsigned int info         = (2<<25)|(chain<<24);
	unsigned int l, r;
	unsigned int rev_l, rev_r;
	unsigned int i;

	unsigned char ec;
	unsigned int sizeZ = keyLength - sizeX - sizeY;
	unsigned int occCount_pstart[4];
	unsigned int occCount_pend[4];
	unsigned int occCountp[4];
	unsigned int occCount_start[4];
	unsigned int occCount_end[4];
	unsigned int occCount[4];

	unsigned long long packedPattern = 0;
	unsigned long long l_packedPattern = 0;
	unsigned long long r_packedPattern = 0;
	unsigned long long rev_packedPattern = 0;
	unsigned long long rev_l_packedPattern = 0;
	unsigned long long rev_r_packedPattern = 0;

	unsigned long long mask;
	unsigned long long ALLONE = (1<<(lookupTable.tableSize*2))-1;
	unsigned char c;
	int hitcount = 0;
	const int coord = (info>>24&1)?(bo->seqLen-bo->alnLen):0;
//	fprintf(stdout, "BWT2misMatching\n");
	//Separate into 4 cases according to the documentation.
	//==============================================


	//Case A    Backward Search
	//1.    cellZ
	//2.    2-mismatch cellX+Y
	//==============================================
	// look-up the last characters (backward)

	for (i = 0; i <lookupTable.tableSize ; i++) {
		packedPattern<<=2;
		packedPattern |= (convertedKey[keyLength-lookupTable.tableSize+i] & 3);
	}

	l = packedPattern ? lookupTable.table[packedPattern-1]+1 : 1;
	r = lookupTable.table[packedPattern];

	// backward search with BWT in cellZ
	for (i = lookupTable.tableSize; i < sizeZ && l <= r; ++i) {
		c = convertedKey[keyLength-i-1];
		l = bwt->cumulativeFreq[c] + BWTOccValue(bwt, l, c) + 1;
		r = bwt->cumulativeFreq[c] + BWTOccValue(bwt, r + 1, c);
	}
	// 2 errors in cellX+Y
	for (i=sizeX+sizeY-1;(i>0 && l<=r);i--) {
		BWTAllOccValue(bwt,l,occCount_pstart);
		BWTAllOccValue(bwt,r + 1,occCount_pend);

		//Backward Manner
		for (ec=0;ec<4;ec++) {
			if (convertedKey[i]==ec) continue;
			unsigned int mk_l=l;
			unsigned int mk_r=r;
			info &= 0x7000000;
			info |= ((((ec&0x3)<<8)|((i+coord)&0xff))&0x3ff)<<12;
			mk_l = bwt->cumulativeFreq[ec] + occCount_pstart[ec] + 1;
			mk_r = bwt->cumulativeFreq[ec] + occCount_pend[ec];

			//return mk_l, mk_r
			if (mk_l <= mk_r) {
				//r_count+=mk_r-mk_l+1;
				hitcount+=BWTBackward1Error(convertedKey,
						bo, bwt,
						0, i,
						mk_l,mk_r,
						info,
						hits);
				if(hits->n >= cutoff) return hitcount;
			}
		}
		c = convertedKey[i];
		l = bwt->cumulativeFreq[c] + occCount_pstart[c] + 1;
		r = bwt->cumulativeFreq[c] + occCount_pend[c];
	}
//	if(hits->n >= cutoff) return hitcount;

//	printf("case A %d\n", saCount);

	//Case B    Forward Search
	//1.    cellX+Y
	//2.    2-mismatch cellZ
	//==============================================

	packedPattern = 0;
	l_packedPattern = 0;
	r_packedPattern = 0;
	rev_packedPattern = 0;
	rev_l_packedPattern = 0;
	rev_r_packedPattern = 0;
	// look-up the first characters
	for (i = 0; i <rev_lookupTable.tableSize ; i++) {
		l_packedPattern<<=2;
		l_packedPattern |= (convertedKey[i] & 3);
	}
//	fprintf(stdout, "l_packedPattern: %u\n", l_packedPattern);
	r_packedPattern = l_packedPattern;
	// if the look-up tables are of different size
	l_packedPattern<<= (lookupTable.tableSize-rev_lookupTable.tableSize)*2;
	r_packedPattern<<= (lookupTable.tableSize-rev_lookupTable.tableSize)*2;
	r_packedPattern |= (1<<(lookupTable.tableSize-rev_lookupTable.tableSize)*2) - 1;

	for (i = 0; i <rev_lookupTable.tableSize ; i++) {
		rev_packedPattern<<= 2;
		rev_packedPattern |= (convertedKey[rev_lookupTable.tableSize-i-1] & 3);
	}
//	fprintf(stdout, "rev_packedPattern: %u\n", rev_packedPattern);
	l     = l_packedPattern ? lookupTable.table[l_packedPattern-1]+1 : 1;
	r     = lookupTable.table[r_packedPattern];
	rev_l = rev_packedPattern ? rev_lookupTable.table[rev_packedPattern-1]+1 : 1;
	rev_r = rev_lookupTable.table[rev_packedPattern];

	// forward search with BWT until the end of forward depth section
	for (i = rev_lookupTable.tableSize; i < sizeX+sizeY && l <= r; ++i) {
		c = convertedKey[i] & 0x3 ;
		BWTAllOccValue(rev_bwt,rev_l,occCount_start);
		BWTAllOccValue(rev_bwt,rev_r + 1,occCount_end);

		rev_l = rev_bwt->cumulativeFreq[c] + occCount_start[c] + 1;
		rev_r = rev_bwt->cumulativeFreq[c] + occCount_end[c];

		occCount[3]=0;
		int k;
		for (k=2;k>=0;k--) {
			occCount[k]=occCount[k+1]+occCount_end[k+1]-occCount_start[k+1];
		}

		r = r - occCount[c];
		l = r - (rev_r-rev_l);
	}

	//2 error in cellZ
//	fprintf(stdout, "find errr\n");
//	fprintf(stdout, "SizeX+Y: %d\n"
//			"keylen: %d\n"
//			"l: %d, r: %d\n", sizeX+sizeY, keyLength, l, r);
	for (i=sizeX+sizeY;(i<keyLength && l<=r);i++) {
		//Call Once Only Proceduressssss - Great
		BWTAllOccValue(rev_bwt,rev_l,occCount_pstart);
		BWTAllOccValue(rev_bwt,rev_r + 1,occCount_pend);

		int k;
		occCountp[3]=0;
		for (k=2;k>=0;k--) {
			occCountp[k]=occCountp[k+1]+occCount_pend[k+1]-occCount_pstart[k+1];
		}

		//Forward Manner
		for (ec=0;ec<4;ec++) {
			if (convertedKey[i]==ec) continue;
//			fprintf(stdout, "%d\n", i);
			info &= 0x7000000;
			info |= ((((ec&0x3)<<8)|((i+coord)&0xff))&0x3ff)<<12;
			unsigned int mk_l=l;
			unsigned int mk_r=r;
			unsigned int rev_mk_l=rev_l;
			unsigned int rev_mk_r=rev_r;

			rev_mk_l = rev_bwt->cumulativeFreq[ec] + occCount_pstart[ec] + 1;
			rev_mk_r = rev_bwt->cumulativeFreq[ec] + occCount_pend[ec];

			mk_r = mk_r - occCountp[ec];
			mk_l = mk_r - (rev_mk_r-rev_mk_l);

			//return mk_l, mk_r
			if (mk_l <= mk_r) {
				//r_count+=mk_r-mk_l+1;
				//2-nd Error Matching in sub cellZ Range
//				fprintf(stdout, "%d--%d--%d\n,",convertedKey[i], i, ec);
				hitcount+=REVBWTForward1Error(convertedKey,
						bo, 
						bwt, rev_bwt,
						i+1,
						keyLength-i-1,
						mk_l,mk_r,rev_mk_l,rev_mk_r,
						info,
						hits);
				if(hits->n >= cutoff) return hitcount;
				//saCount+=forward1Error(p,bwt,rev_bwt, i+1,cellX+cellY+cellZ-i-1,mk_l,mk_r,rev_mk_l,rev_mk_r);
//				fprintf(stdout, "\nsaCoutn:%d\n",saCount);
			}
		}
		c = convertedKey[i];

		rev_l = rev_bwt->cumulativeFreq[c] + occCount_pstart[c] + 1;
		rev_r = rev_bwt->cumulativeFreq[c] + occCount_pend[c];

		r = r - occCountp[c];
		l = r - (rev_r-rev_l);
	}
//	if(hits->n >= cutoff)return hitcount;
	//*/

//	printf("case B %d\n", saCount);
	//

	//Case C
	//1.    cellX (forward)
	//2.    1-mismatch cellY (forward)
	//3.    1-mismatch cellZ (forward)
	//==============================================

	packedPattern = 0;
	l_packedPattern = 0;
	r_packedPattern = 0;
	rev_packedPattern = 0;
	rev_l_packedPattern = 0;
	rev_r_packedPattern = 0;

	for (i = 0; i <lookupTable.tableSize ; i++) {
		packedPattern<<=2;
		packedPattern |= (convertedKey[i] & 3);
	}
	for (i = 0; i <rev_lookupTable.tableSize ; i++) {
		rev_packedPattern<<=2;
		rev_packedPattern |= (convertedKey[rev_lookupTable.tableSize-i-1] & 3 );
	}
	//*
	//For error happen in lookup range....
	//for (i = sizeX; i< lookupTable.tableSize ;i++) {
	for (i = 0; i< lookupTable.tableSize ;i++) {
		unsigned char ec;
		for (ec=0;ec<4;ec++) {
			if (ec == (convertedKey[i] & 0x3)) continue;
			info &= 0x7000000;
			info |= ((((ec&0x3)<<8)|((i+coord)&0xff))&0x3ff)<<12;
			unsigned long long err_packedPattern = packedPattern;
			unsigned long long err_rev_packedPattern = rev_packedPattern;

			unsigned int bitPos = (lookupTable.tableSize-i-1)*2;
			mask = ALLONE - (3 << bitPos);
			mask |= ec << bitPos;
			err_packedPattern |= (3 << bitPos);
			err_packedPattern &= mask;

			bitPos = i*2;
			mask = ALLONE - (3 << bitPos);
			mask |= ec << bitPos;
			err_rev_packedPattern |= (3 << bitPos);
			err_rev_packedPattern &= mask;

			l = err_packedPattern ? lookupTable.table[err_packedPattern-1]+1 : 1;
			r = lookupTable.table[err_packedPattern];
			rev_l = err_rev_packedPattern ? rev_lookupTable.table[err_rev_packedPattern-1]+1 : 1;
			rev_r = rev_lookupTable.table[err_rev_packedPattern];
			//LookupSafe(lookupTable,bwt,err_packedPattern,err_packedPattern,&l,&r);
			//LookupSafe(rev_lookupTable,rev_bwt,err_rev_packedPattern,err_rev_packedPattern,&rev_l,&rev_r);

			unsigned int pos = lookupTable.tableSize;
			while (pos < sizeX+sizeY  && l <= r) {
				c = convertedKey[pos] & 0x3;

				BWTAllOccValue(rev_bwt,rev_l,occCount_start);
				BWTAllOccValue(rev_bwt,rev_r + 1,occCount_end);

				rev_l = rev_bwt->cumulativeFreq[c] + occCount_start[c] + 1;
				rev_r = rev_bwt->cumulativeFreq[c] + occCount_end[c];

				int k;
				occCount[3]=0;
				for (k=2;k>=0;k--) {
					occCount[k]=occCount[k+1]+occCount_end[k+1]-occCount_start[k+1];
				}

				r = r - occCount[c];
				l = r - (rev_r-rev_l);

				pos++;
			}


			if (l <= r) {
				//2-nd Error Matching in cellZ Range
				hitcount+=REVBWTForward1Error(convertedKey, 
								bo,
								bwt, rev_bwt,
								sizeX+sizeY, sizeZ,
								l,r,
								rev_l,rev_r,
								info,
								hits);
				if(hits->n >= cutoff) return hitcount;
			}
		}
	}
	//*/
	if(hits->n >= cutoff)return hitcount;


	l = packedPattern ? lookupTable.table[packedPattern-1]+1 : 1;
	r = lookupTable.table[packedPattern];
	rev_l = rev_packedPattern ? rev_lookupTable.table[rev_packedPattern-1]+1 : 1;
	rev_r = rev_lookupTable.table[rev_packedPattern];


	//For error happen outside lookup range..
	for (i=lookupTable.tableSize;(i<sizeX+sizeY && l<=r);i++) {
		//Call Once Only Proceduressssss - Great
		BWTAllOccValue(rev_bwt,rev_l,occCount_pstart);
		BWTAllOccValue(rev_bwt,rev_r + 1,occCount_pend);

		int k;
		occCountp[3]=0;
		for (k=2;k>=0;k--) {
			occCountp[k]=occCountp[k+1]+occCount_pend[k+1]-occCount_pstart[k+1];
		}
		//Forward Manner
		for (ec=0;ec<4;ec++) {
			if ((convertedKey[i]&0x3) ==ec) continue;
			info &= 0x7000000;
			info |= ((((ec&0x3)<<8)|((i+coord)&0xff))&0x3ff)<<12;
			unsigned int mk_l=l;
			unsigned int mk_r=r;
			unsigned int rev_mk_l=rev_l;
			unsigned int rev_mk_r=rev_r;

			unsigned int pos = i+1;

			rev_mk_l = rev_bwt->cumulativeFreq[ec] + occCount_pstart[ec] + 1;
			rev_mk_r = rev_bwt->cumulativeFreq[ec] + occCount_pend[ec];

			mk_r = mk_r - occCountp[ec];
			mk_l = mk_r - (rev_mk_r-rev_mk_l);

			while (pos < sizeX+sizeY  && rev_mk_l <= rev_mk_r) {
				c = convertedKey[pos] & 0x3;

				BWTAllOccValue(rev_bwt,rev_mk_l,occCount_start);
				BWTAllOccValue(rev_bwt,rev_mk_r + 1,occCount_end);

				rev_mk_l = rev_bwt->cumulativeFreq[c] + occCount_start[c] + 1;
				rev_mk_r = rev_bwt->cumulativeFreq[c] + occCount_end[c];

				int k;
				occCount[3]=0;
				for (k=2;k>=0;k--) {
					occCount[k]=occCount[k+1]+occCount_end[k+1]-occCount_start[k+1];
				}

				mk_r = mk_r - occCount[c];
				mk_l = mk_r - (rev_mk_r-rev_mk_l);

				pos++;
			}
			//return mk_l, mk_r
			if (mk_l <= mk_r) {
				//2-nd Error Matching in cellZ Range
				hitcount+=REVBWTForward1Error(convertedKey, 
						bo,
						bwt, rev_bwt,
						sizeX+sizeY,
						sizeZ,
						mk_l,mk_r,
						rev_mk_l,rev_mk_r,
						info,
						hits);
//				if(hits->n >= cutoff) return hitcount;
			}
		}
		c = convertedKey[i];

		rev_l = rev_bwt->cumulativeFreq[c] + occCount_pstart[c] + 1;
		rev_r = rev_bwt->cumulativeFreq[c] + occCount_pend[c];

		r = r - occCountp[c];
		l = r - (rev_r-rev_l);
	}
	if(hits->n >= cutoff)return hitcount;
	//*/

//	printf("case C %d\n", saCount);

/*
	//Case D
	//1.    cellY (forward)
	//2.    1-mismatch cellZ (forward)
	//3.    1-mismatch cellX (backward)
	//==============================================
	packedPattern = 0;
	l_packedPattern = 0;
	r_packedPattern = 0;
	rev_packedPattern = 0;
	rev_l_packedPattern = 0;
	rev_r_packedPattern = 0;

	for (i = 0; i <lookupTable.tableSize ; i++) {
		packedPattern<<=2;
		packedPattern |= (convertedKey[sizeX+i] & 3);
	}
	for (i = 0; i <rev_lookupTable.tableSize ; i++) {
		rev_packedPattern<<=2;
		rev_packedPattern |= (convertedKey[sizeX+rev_lookupTable.tableSize-i-1] & 3);
	}
	//*
	//For error happen in lookup range....
	for (i = sizeY; i< lookupTable.tableSize ;i++) {
		unsigned char ec;
		for (ec=0;ec<4;ec++) {
			if (ec == convertedKey[sizeX+i]) continue;
			info &= 0x7000000;
			info |= ((((ec&0x3)<<8)|((sizeX+i)&0xff))&0x3ff)<<12;
			unsigned long long err_packedPattern = packedPattern;
			unsigned long long err_rev_packedPattern = rev_packedPattern;

			unsigned int bitPos = (lookupTable.tableSize-i-1)*2;
			mask = ALLONE - (3 << bitPos);
			mask |= ec << bitPos;
			err_packedPattern |= (3 << bitPos);
			err_packedPattern &= mask;

			bitPos = i*2;
			mask = ALLONE - (3 << bitPos);
			mask |= ec << bitPos;
			err_rev_packedPattern |= (3 << bitPos);
			err_rev_packedPattern &= mask;

			l = err_packedPattern ? lookupTable.table[err_packedPattern-1]+1 : 1;
			r = lookupTable.table[err_packedPattern];
			rev_l = err_rev_packedPattern ? rev_lookupTable.table[err_rev_packedPattern-1]+1 : 1;
			rev_r = rev_lookupTable.table[err_rev_packedPattern];


			unsigned int pos = sizeX+lookupTable.tableSize;
			while (pos < keyLength  && l <= r) {
				c = convertedKey[pos] & 0x3 ;

				BWTAllOccValue(rev_bwt,rev_l,occCount_start);
				BWTAllOccValue(rev_bwt,rev_r + 1,occCount_end);

				rev_l = rev_bwt->cumulativeFreq[c] + occCount_start[c] + 1;
				rev_r = rev_bwt->cumulativeFreq[c] + occCount_end[c];

				int k;
				occCount[3]=0;
				for (k=2;k>=0;k--) {
					occCount[k]=occCount[k+1]+occCount_end[k+1]-occCount_start[k+1];
				}

				r = r - occCount[c];
				l = r - (rev_r-rev_l);

				pos++;
			}


			if (l <= r) {
				//2-nd Error Matching in cellX Range
				hitcount+=BWTBackward1Error(convertedKey,bo, 
										bwt,
										0, sizeX,
										l,r,
										info,
										hits);
			}
			if (hits->n >= cutoff) return hitcount;
		}
	}
	///

	l = packedPattern ? lookupTable.table[packedPattern-1]+1 : 1;
	r = lookupTable.table[packedPattern];
	rev_l = rev_packedPattern ? rev_lookupTable.table[rev_packedPattern-1]+1 : 1;
	rev_r = rev_lookupTable.table[rev_packedPattern];


    //For error happen outside lookup range..
	for (i=sizeX+lookupTable.tableSize;(i<keyLength && l<=r);i++) {
		BWTAllOccValue(rev_bwt,rev_l,occCount_pstart);
		BWTAllOccValue(rev_bwt,rev_r + 1,occCount_pend);

		int k;
		occCountp[3]=0;
		for (k=2;k>=0;k--) {
			occCountp[k]=occCountp[k+1]+occCount_pend[k+1]-occCount_pstart[k+1];
		}
		//Forward Manner
		for (ec=0;ec<4;ec++) {
			if (convertedKey[i]==ec) continue;
			info &= 0x7000000;
			info |= ((((ec&0x3)<<8)|(i&0xff))&0x3ff)<<12;
			unsigned int mk_l=l;
			unsigned int mk_r=r;
			unsigned int rev_mk_l=rev_l;
			unsigned int rev_mk_r=rev_r;

			unsigned int pos = i+1;

			rev_mk_l = rev_bwt->cumulativeFreq[ec] + occCount_pstart[ec] + 1;
			rev_mk_r = rev_bwt->cumulativeFreq[ec] + occCount_pend[ec];

			mk_r = mk_r - occCountp[ec];
			mk_l = mk_r - (rev_mk_r-rev_mk_l);

			while (pos < keyLength  && rev_mk_l <= rev_mk_r) {
				c = convertedKey[pos];

				BWTAllOccValue(rev_bwt,rev_mk_l,occCount_start);
				BWTAllOccValue(rev_bwt,rev_mk_r + 1,occCount_end);

				rev_mk_l = rev_bwt->cumulativeFreq[c] + occCount_start[c] + 1;
				rev_mk_r = rev_bwt->cumulativeFreq[c] + occCount_end[c];

				int k;
				occCount[3]=0;
				for (k=2;k>=0;k--) {
					occCount[k]=occCount[k+1]+occCount_end[k+1]-occCount_start[k+1];
				}

				mk_r = mk_r - occCount[c];
				mk_l = mk_r - (rev_mk_r-rev_mk_l);

				pos++;
			}
			//return mk_l, mk_r
			if (mk_l <= mk_r) {
				//2-nd Error Matching in cellX Range
				hitcount+=BWTBackward1Error(convertedKey, bo,
						bwt, 0,
						sizeX,
						mk_l,mk_r,
						info,
						hits);
			}
			if (hits->n >= cutoff)return hitcount;
		}
		c = convertedKey[i];

		rev_l = rev_bwt->cumulativeFreq[c] + occCount_pstart[c] + 1;
		rev_r = rev_bwt->cumulativeFreq[c] + occCount_pend[c];

		r = r - occCountp[c];
		l = r - (rev_r-rev_l);
	}
	///
//	printf("case D %d\n", saCount);
	//*/
	return hitcount;

}

static inline int POSCMP(const void *a, const void *b){
	return *(unsigned int *)a - *(unsigned int *)b;
}
