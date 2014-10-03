#include <stdio.h>
#include <string.h>
#include <stdlib.h>

//#include "common.h"
#include "alignment.h"

#define ARRAY(a, r, c) a[(r)*row_len + (c)]

//enum{EDNAFULL, NW_MATRIX, SJ_MATRIX};
double gsm[3][4][4] = {
/* EDNAFULL matrix from ftp://ftp.ncbi.nih.gov/blast/matrices/NUC.4.4 */
	{
		{5.0, -4.0, -4.0, -4.0},
		{-4.0, 5.0, -4.0, -4.0},
		{-4.0, -4.0, 5.0, -4.0},
		{-4.0, -4.0, -4.0, 5.0}
	},

/* SJ_MATRIX If you really want to avoid mismatches, you may use this matrix. */
	{
		{1.0, -10000.0, -10000.0, -10000.0},
		{-10000.0, 1.0, -10000.0, -10000.0},
		{-10000.0, -10000.0, 1.0, -10000.0},
		{-10000.0, -10000.0, -10000.0, 1.0}
	},

/* NW_MATRIX From http://en.wikipedia.org/wiki/Needleman-Wunsch_algorithm */
	{
		{10.0, -1.0, -3.0, -4.0},
		{-1.0, 7.0, -5.0, -3.0},
		{-3.0, -5.0, 9.0, 0.0},
		{-4.0, -3.0, 0.0, 8.0}
	}
};
#if 0
/* EDNAFULL matrix from ftp://ftp.ncbi.nih.gov/blast/matrices/NUC.4.4 */
double sm1[4][4] = { /* similarity matrix */
	{5.0, -4.0, -4.0, -4.0},
	{-4.0, 5.0, -4.0, -4.0},
	{-4.0, -4.0, 5.0, -4.0},
	{-4.0, -4.0, -4.0, 5.0}
};

/* From http://en.wikipedia.org/wiki/Needleman-Wunsch_algorithm */
double sm2[4][4] = { /* similarity matrix */
	{10.0, -1.0, -3.0, -4.0},
	{-1.0, 7.0, -5.0, -3.0},
	{-3.0, -5.0, 9.0, 0.0},
	{-4.0, -3.0, 0.0, 8.0}
};

/* If you really want to avoid mismatches, you may use this matrix. */
double sm3[4][4] = { /* similarity matrix */
	{1.0, -10000.0, -10000.0, -10000.0},
	{-10000.0, 1.0, -10000.0, -10000.0},
	{-10000.0, -10000.0, 1.0, -10000.0},
	{-10000.0, -10000.0, -10000.0, 1.0}
};
#endif

double gop = -10.0; /* gap opening penalty */
double gec = -0.5;
#define INF	1000000000

#define PRINT_ARRAY(a) {\
	int ii, jj;\
	for(ii = 0; ii < (int)strlen(sequence)+1; ii++){\
		for(jj = 0; jj < (int)strlen(reference)+1; jj++)\
			printF("%4.1f\t", ARRAY(a, ii, jj));\
		printF("\n");\
	}\
}

/*  --------------------------- reference (t in the book)
	|
	|
	|
	|
	sequence (s in the book)

Reference: Setubal, Meidanis, "Introduction to Computational Molecular Biology,
		   Chap. 3.2.1 "GLOBAL COMPARISON - THE BASIC ALGORITHM" (pp. 49 - 55) */
alignment NW_alignment(const char *sequence, const char *reference, alignment align, int matrix, int print)
{
	double *a = NULL, val1 = 0, val2 = 0, val3 = 0, sm[4][4];
	double s = 0, s_diag = 0, s_left = 0, s_up = 0;
	int n = strlen(reference) + 1, m = strlen(sequence) + 1, row_len = strlen(reference) + 1;
	int row = 0, column = 0, ri = n+m-1, si = n+m-1, i = 0;
	unsigned int ref_len = 0, seq_len = 0;
	char *ref = NULL, *seq = NULL;

	gop = -5; /* from http://en.wikipedia.org/wiki/Needleman-Wunsch_algorithm */
	memset(&align, 0x00, sizeof(alignment));
	memcpy(sm, gsm[matrix], sizeof(double)*4*4);

	MALLOC(a, double, n*m);
	MALLOC(ref, char, n+m);
	MALLOC(seq, char, n+m);

	for(row = 0; row < m; row++) ARRAY(a, row, 0) = (double)row*gop;
	for(column = 0; column < n; column++) ARRAY(a, 0, column) = (double)column*gop;

	/* DP for score */
	for(row = 1; row < m; row++){
		for(column = 1; column < n; column++){
			val1 = ARRAY(a, row-1, column-1) + sm[b2num(sequence[row-1])][b2num(reference[column-1])];
			val2 = ARRAY(a, row-1, column) + gop;
			val3 = ARRAY(a, row, column-1) + gop;

			ARRAY(a, row, column) = MAX(MAX(val1, val2), val3);
		}
	}
	//PRINT_ARRAY(a);

	/* back tracking */
	row = m-1, column = n-1;
	while((row > 0) && (column > 0)){
		s = ARRAY(a, row, column);
		s_diag = ARRAY(a, row-1, column-1);
		s_left = ARRAY(a, row, column-1);
		s_up = ARRAY(a, row-1, column);

		if(s == (s_diag + sm[b2num(sequence[row-1])][b2num(reference[column-1])])){
			ref[ri--] = reference[column-1];
			ref_len++;
			seq[si--] = sequence[row-1];
			row--;
			column--;
		}
		else if(s == (s_left+gop)){
			ref[ri--] = reference[column-1];
			ref_len++;
			seq[si--] = '-';
			column--;
		}
		else if(s == (s_up+gop)){
			ref[ri--] = '-';
			ref_len++;
			seq[si--] = sequence[row-1];
			row--;
		}
		else
			printf("Unexpected Error\n");

	}

	while(column > 0){
		ref[ri--] = reference[column-1];
		ref_len++;
		seq[si--] = '-';
		column--;
	}

	while(row > 0){
		ref[ri--] = '-';
		ref_len++;
		seq[si--] = sequence[row-1];
		row--;
	}

	//printInfo("\nref=%s\nseq=%s\n%f\n", &(ref[ri+1]), &(seq[si+1]), s);
	if(ri != si)
		printf("Error2\n");

	/* calculate the lengths */
	seq_len = ref_len;
	for(i = ri+1; i < n+m; i++){ /* eat first '-' */
		if((IS_VALID_BASE(ref[i]) == 1)||(IS_VALID_BASE(seq[i]) == 1)) break;
		ref_len--;
		seq_len--;
	}
	for(i = n+m-1; i > ri; i--){ /* eat last '-'  */
		if((IS_VALID_BASE(ref[i]) == 1)||(IS_VALID_BASE(seq[i]) == 1)) break;
		ref_len--;
		seq_len--;
	}
	if(seq_len!=ref_len)
	    printf("Unexpected Error 2.\n");
	//align.s1_len = seq_len;
	//align.s2_len = ref_len;
	//align.long_len = (ref_len > seq_len)? ref_len : seq_len;
	//align.short_len = (ref_len > seq_len)? seq_len : ref_len;
    
	bool tempflagref=false,tempflagseq=false;
	for(i = ri+1; i < n+m; i++){
		if((IS_VALID_BASE(ref[i]) == NO) && (IS_VALID_BASE(seq[i]) == NO)) continue;
		if(IS_VALID_BASE(ref[i]) == NO){
            tempflagseq=false;
            if(tempflagref==false){
                 tempflagref=true;
            }
            else{ 
                 ref_len--;
            }
        }
        if(IS_VALID_BASE(seq[i]) == NO){
            tempflagref=false;
            if(tempflagseq==false){
                 tempflagseq=true;
            }
            else{
                 seq_len--;
            }
        }      
		if(ref[i] == seq[i]){
            align.matches++;
            tempflagref=false;
            tempflagseq=false;
        }
	}
	
	//align.long_len = (ref_len > seq_len)? ref_len : seq_len;
	//align.short_len = (ref_len > seq_len)? seq_len : ref_len;
	align.s1_len = seq_len;
	align.s2_len = ref_len;
	//printInfo("Similarity(%d/%d): %f%%\n", align.matches, n+m-2-ri+1, (double)align.matches/(double)(n+m-2-ri+1)*100.0);

	FREE(a);
	FREE(ref);
	FREE(seq);
	//return (double)align.matches/(double)(n+m-ri);
	return align;
}

/*
Reference: Setubal, Meidanis, "Introduction to Computational Molecular Biology,
		   Chap. 3.3.3 AFFINE GAP PENALTY FUNCTIONS (pp. 64 - 66)
This function may get exactly same performance with program needle of EMBOSS package. */
alignment NW_alignmentEx(const char *sequence, const char *reference, alignment align, int matrix, int print)
{
	double *a[3] = {NULL, NULL, NULL}, val1 = 0, val2 = 0, val3 = 0, sm[4][4];
	int n = strlen(reference) + 1, m = strlen(sequence) + 1, row_len = strlen(reference) + 1;
	int colmax=n, colmin=1;
	bool colmaxflag=false,colminflag=false;
	int threshold = (int)((double)MAX(m,n)*0.05);
	int row = 0, column = 0, ri = n+m-1, si = n+m-1, i = 0, j = 0, type = -1;
	int ref_len = 0, seq_len = 0;
	int tempri=0, tempsi=0, tempriend=0, tempsiend=0;
	char *ref = NULL, *seq = NULL;

	memset(&align, 0x00, sizeof(alignment));
	memcpy(sm, gsm[matrix], sizeof(double)*4*4);

	MALLOC(a[0], double, n*m);
	ARRAY(a[0], 0, 0) = 0.0;
	for(i = 1; i < m; i++) ARRAY(a[0], i, 0) = -INF;
	for(j = 1; j < n; j++) ARRAY(a[0], 0, j) = -INF;
	MALLOC(a[1], double, n*m);
	for(i = 0; i < m; i++) ARRAY(a[1], i, 0) = -INF;
	for(j = 1; j < n; j++) ARRAY(a[1], 0, j) = gop + j*gec;
	MALLOC(a[2], double, n*m);
	for(i = 1; i < m; i++) ARRAY(a[2], i, 0) = gop + i*gec;;
	for(j = 0; j < n; j++) ARRAY(a[2], 0, j) = -INF;
	MALLOC(ref, char, n+m);
	MALLOC(seq, char, n+m);

	/* DP for score */
	for(row = 1; row < m; row++){
        if(m>=n)
        {
            if((row-threshold+n-m)>1)
                colminflag=true;
            if(1+row+threshold<=n)
                colmaxflag=true;    
            colmin=MAX(1,row-threshold+n-m);
            colmax=MIN(n,1+row+threshold);
        }
        if(m<n)
        {
            if((row-threshold)>1)
                colminflag=true;
            if((n-m+1+row+threshold)<=n)
                colmaxflag=true; 
            colmin=MAX(1,row-threshold);
            colmax=MIN(n,n-m+1+row+threshold);
        }    
		for(column = colmin; column < colmax; column++){
			/* update a */
			val1 = ARRAY(a[0], row-1, column-1) + sm[b2num(sequence[row-1])][b2num(reference[column-1])];
			val2 = ARRAY(a[1], row-1, column-1) + sm[b2num(sequence[row-1])][b2num(reference[column-1])];
			val3 = ARRAY(a[2], row-1, column-1) + sm[b2num(sequence[row-1])][b2num(reference[column-1])];

			ARRAY(a[0], row, column) = MAX(MAX(val1, val2), val3);

			/* update b */
			if((colminflag==true)&&(column==colmin))
			     ARRAY(a[1], row, column) = -INF;
			     
            else
            {
			val1 = ARRAY(a[0], row, column-1) + gop;
			val2 = ARRAY(a[1], row, column-1) + gec;
			val3 = ARRAY(a[2], row, column-1) + gop;

			ARRAY(a[1], row, column) = MAX(MAX(val1, val2), val3);
            }
                

			/* update c */
			if((colmaxflag==true)&&(column==colmax))
                 ARRAY(a[2], row, column) = -INF;
            
            else
            {
			val1 = ARRAY(a[0], row-1, column) + gop;
			val2 = ARRAY(a[1], row-1, column) + gop;
			val3 = ARRAY(a[2], row-1, column) + gec;

			ARRAY(a[2], row, column) = MAX(MAX(val1, val2), val3);
            }
               
		}
	}
	//PRINT_ARRAY(a[0]);
	//PRINT_ARRAY(a[1]);
	//PRINT_ARRAY(a[2]);

#define TYPE(r, c, p0, p1, p2)	(((ARRAY(a[0], (r), (c)) + p0) > (ARRAY(a[1], (r), (c)) + p1))?\
		(((ARRAY(a[0], (r), (c)) + p0) > (ARRAY(a[2], (r), (c)) + p2))? 0:2):\
		(((ARRAY(a[1], (r), (c)) + p1) > (ARRAY(a[2], (r), (c)) + p2))? 1:2))

	/* back tracking */
	row = m-1, column = n-1;
	type = TYPE(row, column, 0, 0, 0);
	while((row > 0) && (column > 0)){

		switch(type){
			case 0: /* From diagonal */
				ref[ri--] = reference[column-1];
				ref_len++;
				seq[si--] = sequence[row-1];

				row--;
				column--;

				type = TYPE(row, column, 0, 0, 0);
				continue;
			case 1: /* From left */
				ref[ri--] = reference[column-1];
				ref_len++;
				seq[si--] = '-';
			
				column--;

				type = TYPE(row, column, gop, gec, gop);
				continue;
			case 2: /* From upward */
				ref[ri--] = '-';
				ref_len++;
				seq[si--] = sequence[row-1];

				row--;

				type = TYPE(row, column, gop, gop, gec);
				continue;
			default:
				printf("Type was errorneous.\n");
		}
	}

	while(column > 0){
		ref[ri--] = reference[column-1];
		ref_len++;
		seq[si--] = '-';
		column--;
	}

	while(row > 0){
		ref[ri--] = '-';
		ref_len++;
		seq[si--] = sequence[row-1];
		row--;
	}

	//printInfo("\nref=%s\nseq=%s\n", &(ref[ri+1]), &(seq[si+1]));
	if(ri != si)
		printf("Unexpected Error\n");

	/* calculate the lengths */
	seq_len = ref_len;
	//align.s1_len = seq_len;
	//align.s2_len = ref_len;
	//align.long_len = (ref_len > seq_len)? ref_len : seq_len;
	//align.short_len = (ref_len > seq_len)? seq_len : ref_len;
	for(i = ri+1; i < n+m; i++){ /* eat first '-' */
		if(IS_VALID_BASE(ref[i]) == YES) 
        {
            tempri=i;
            break;
        }
		ref_len=MAX(0,ref_len-1);
	}
	for(i = n+m-1; i > ri; i--){ /* eat last '-' */
		if(IS_VALID_BASE(ref[i]) == YES)
        {
            tempriend=i;
            break;
        }
		ref_len=MAX(0,ref_len-1);	
	}
	for(i = si+1; i < n+m; i++){ /* eat first '-' */
		if(IS_VALID_BASE(seq[i]) == YES)
        {
            tempsi=i;
            break;
        }
		seq_len=MAX(0,seq_len-1);
	}
	for(i = n+m-1; i > si; i--){ /* eat last '-' */
		if(IS_VALID_BASE(seq[i]) == YES)
        {
            tempsiend=i;
            break;
        }
		seq_len=MAX(0,seq_len-1);	
	}    
	bool tempflagref=false,tempflagseq=false;
    for(i = MAX(tempsi, tempri); i <= MIN(tempsiend,tempriend); i++){
		if((IS_VALID_BASE(ref[i]) == NO) && (IS_VALID_BASE(seq[i]) == NO)) continue;
		if(IS_VALID_BASE(ref[i]) == NO){
            tempflagseq=false;
            if(tempflagref==false){
                 tempflagref=true;
            }
            else{ 
                 ref_len=MAX(0,ref_len-1);
                 seq_len=MAX(0,seq_len-1);           
		}
        }
        if(IS_VALID_BASE(seq[i]) == NO){
            tempflagref=false;
            if(tempflagseq==false){
                 tempflagseq=true;
            }
            else{
                 seq_len=MAX(0,seq_len-1);                 
		   ref_len=MAX(0,ref_len-1);
            }
        }      
		if(ref[i] == seq[i]){
            align.matches++;
            tempflagref=false;
            tempflagseq=false;
        }
	}
	align.s1_len = MAX(1,seq_len);
	align.s2_len = MAX(1,ref_len);
	//align.long_len = (ref_len > seq_len)? ref_len : seq_len;
	//align.short_len = (ref_len > seq_len)? seq_len : ref_len;

	//printInfo("Similarity(%d/%d): %f%%\n", align.matches, n+m-2-ri+1, (double)align.matches/(double)(n+m-2-ri+1)*100.0);

	FREE(a[0]);
	FREE(a[1]);
	FREE(a[2]);
	FREE(ref);
	FREE(seq);
	//return (double)same/(double)(n+m-ri);
	return align;
}

/*
Reference: Michael S. Waterman, "Introduction to Computational Biology - Maps, sequences and genomes,
		   Chap. 9.6 Local Alignment and Clumps (pp. 202 - 206)
This is the local alignment function. */

/* NOTE: the indexing of sequences and matrices are different. The first row and column of matrices
   are not related to any base. That is, i-th row does not mean reference[i] base, it means reference[i-1]. */
alignment SW_alignmentEx(const char *sequence /* seq 1 */, const char *reference /* seq 2 */, alignment align, int matrix, int print)
{
	double *a[3] = {NULL, NULL, NULL}, val1 = 0, val2 = 0, sm[4][4];
	int n = strlen(reference) + 1, m = strlen(sequence) + 1, row_len = strlen(reference) + 1;
	int row = 0, column = 0, ri = n+m-1, si = n+m-1, i = 0, type = -1;
	int rowMax = 0, columnMax = 0;
	unsigned int ref_len = 0, seq_len = 0;
	char *ref = NULL, *seq = NULL;

	memset(&align, 0x00, sizeof(alignment));
	memcpy(sm, gsm[matrix], sizeof(double)*4*4);
	/*
	if(matrix == 1)
		memcpy(sm, sm1, sizeof(double)*4*4);
	else
		memcpy(sm, sm3, sizeof(double)*4*4);
	*/

	MALLOC(a[0], double, n*m); // E
	MALLOC(a[1], double, n*m); // F
	MALLOC(a[2], double, n*m); // H
	MALLOC(ref, char, n+m+1);
	MALLOC(seq, char, n+m+1);

	/* DP for score */
	for(row = 1; row < m; row++){
		for(column = 1; column < n; column++){
			/* update E */
			val1 = ARRAY(a[0], row, column-1) + gec;
			val2 = ARRAY(a[2], row, column-1) + gop;

			ARRAY(a[0], row, column) = MAX(val1, val2);

			/* update F */
			val1 = ARRAY(a[1], row-1, column) + gec;
			val2 = ARRAY(a[2], row-1, column) + gop;

			ARRAY(a[1], row, column) = MAX(val1, val2);

			/* update H */
			val1 = ARRAY(a[2], row-1, column-1) + sm[b2num(sequence[row-1])][b2num(reference[column-1])];

			if((ARRAY(a[2], row, column) = MAX(MAX(val1, 0), MAX(ARRAY(a[1], row, column), ARRAY(a[0], row, column)))) >= ARRAY(a[2], rowMax, columnMax)){
				rowMax = row;
				columnMax = column;
			}
		}
	}
	//PRINT_ARRAY(a[0]);
	//PRINT_ARRAY(a[1]);
	//PRINT_ARRAY(a[2]);

#define TYPESW(r, c, p0, p1, p2)	(((ARRAY(a[1], (r), (c)) + p1) >= (ARRAY(a[2], (r-1), (c-1)) + p2))?\
		(((ARRAY(a[1], (r), (c)) + p1) >= (ARRAY(a[0], (r), (c)) + p0))? 1:0):\
		(((ARRAY(a[2], (r-1), (c-1)) + p2) >= (ARRAY(a[0], (r), (c)) + p0))? 2:0))

	/* back tracking */
	row = rowMax, column = columnMax;
	align.s1_e = rowMax-1; /* This must be the offset from the first base, so you should subtract 1. */
	align.s2_e = columnMax-1;
	type = TYPESW(row, column, 0, 0, sm[b2num(sequence[row-1])][b2num(reference[column-1])]);
	while((row > 0) && (column > 0)){
		if(ARRAY(a[2], row, column) <= 0)
			break;

		/*
		if(ARRAY(a[0], row, column) == ARRAY(a[1], row, column) ||
				ARRAY(a[0], row, column) == ARRAY(a[2], row-1, column-1)+sm[b2num(sequence[row-1])][b2num(reference[column-1])] ||
				ARRAY(a[1], row, column) == ARRAY(a[2], row-1, column-1)+sm[b2num(sequence[row-1])][b2num(reference[column-1])])
			printInfo("There is tie!\n");
		*/

		if((IS_VALID_BASE(reference[column-1]) == 0) || (IS_VALID_BASE(sequence[row-1]) == 0))
			printf("Not a valid base was found (%X:%X)\n", reference[column-1], sequence[row-1]);

		switch(type){
			case 0: /* From left */
				ref[ri--] = reference[column-1];
				ref_len++;
				seq[si--] = '-';
			
				column--;

				break;
			case 1: /* From upward */
				ref[ri--] = '-';
				ref_len++;
				seq[si--] = sequence[row-1];

				row--;

				break;
			case 2: /* From diagonal */
				ref[ri--] = reference[column-1];
				ref_len++;
				seq[si--] = sequence[row-1];

				row--;
				column--;

				break;
			default:
				printf("Type was errorneous.\n");
				break;
		}

		if((row > 0) && (column > 0))
			type = TYPESW(row, column, 0, 0, sm[b2num(sequence[row-1])][b2num(reference[column-1])]);
		if(ARRAY(a[type], row, column) != ARRAY(a[2], row, column) && ARRAY(a[2], row, column) != 0)
			printf("Weird\n");
	}
	if((row > 0) && (column > 0)) align.s1_s = row, align.s2_s = column;
	else align.s1_s = row, align.s2_s = column;
	//align.s1_s = (row > 0)? row-1:0;
	//align.s2_s = (column > 0)? column-1:0;

	//printInfo("\nref=%s\nseq=%s\nScore = %f\n", &(ref[ri+1]), &(seq[si+1]), ARRAY(a[2], rowMax, columnMax));
	if(ri != si)
		printf("Unexpected Error\n");

	/* calculate the lengths */
	seq_len = ref_len;
	for(i = ri+1; i < n+m; i++){ /* eat first '-' */
		if(IS_VALID_BASE(ref[i]) == YES) break;
		ref_len--;
	}
	for(i = n+m-1; i > ri; i--){ /* eat last '-' */
		if(IS_VALID_BASE(ref[i]) == YES) break;
		ref_len--;
	}
	for(i = si+1; i < n+m; i++){ /* eat first '-' */
		if(IS_VALID_BASE(seq[i]) == YES) break;
		seq_len--;
	}
	for(i = n+m-1; i > si; i--){ /* eat last '-' */
		if(IS_VALID_BASE(seq[i]) == YES) break;
		seq_len--;
	}
	align.s1_len = seq_len;
	align.s2_len = ref_len;
	//align.long_len = (ref_len > seq_len)? ref_len : seq_len;
	//align.short_len = (ref_len > seq_len)? seq_len : ref_len;

	for(i = ri+1; i < n+m; i++){
		if((IS_VALID_BASE(ref[i]) == NO) && (IS_VALID_BASE(seq[i]) == NO)) continue;
		if(ref[i] == seq[i]) align.matches++;
	}
	//printInfo("Similarity(%d/%d): %f%%\n", align.matches, n+m-2-ri+1, (double)align.matches/(double)(n+m-2-ri+1)*100.0);

	FREE(a[0]);
	FREE(a[1]);
	FREE(a[2]);
	FREE(ref);
	FREE(seq);
	//return (double)same/(double)(n+m-ri);
	return align;
}

int align_type(alignment align, char *_seq1, char *_seq2)
{
	unsigned int seq1_len = 0U, seq2_len = 0U;
	int left_offset = 0, right_offset = 0;
	char *seq1 = NULL, *seq2 = NULL;

	if(strlen(_seq1) > strlen(_seq2)) seq1 = _seq1, seq2 = _seq2;
	else seq1 = _seq2, seq2 = _seq1;
	seq1_len = strlen(seq1), seq2_len = strlen(seq2);

	left_offset = (int)align.s1_s - (int)align.s2_s;
	right_offset = (int)align.s1_e + ((int)seq2_len - (int)align.s2_e - 1);

	if(left_offset < 0) return LEFT_PROTRUDING;
	else if(right_offset >= seq1_len) return RIGHT_PROTRUDING;
	else return MIDDLE;
}

char* merge_align(alignment align, char *_seq1, char *_seq2, unsigned int min_len /* This might be Klen. */)
{
	unsigned int seq1_len = 0U, seq2_len = 0U;
	int left_offset = 0, right_offset = 0;
	char *merge = NULL, *seq1 = NULL, *seq2 = NULL;

	if(strlen(_seq1) < strlen(_seq2)){
		swapAlign(&align);
		seq1 = _seq2, seq2 = _seq1;
	}
	else seq1 = _seq1, seq2 = _seq2;
	seq1_len = strlen(seq1), seq2_len = strlen(seq2);

	if((seq1_len < min_len) || (seq2_len < min_len)) return strdup(seq1);
	
	if(MIN(seq1_len, seq2_len) == 831 || MIN(seq1_len, seq2_len) == 833) printAlign(align);
	//printAlign(align);

	left_offset = (int)align.s1_s - (int)align.s2_s;
	right_offset = (int)align.s1_e + ((int)seq2_len - (int)align.s2_e - 1);

	if(left_offset < 0){
		/* When longer one is better and they are not so similar to each other, avoid the edge part of shorter one. */
		if(align.s1_supp_read_den >= align.s2_supp_read_den){ /* Prefer longer one */
			/*
		                   	----================================
					============--------------
				              	^        ^
			*/
			if(align.s1_s+1 < min_len){
				align.s2_s += (min_len-align.s1_s-1);
				align.s1_s = min_len-1;
			}
			MALLOC(merge, char, align.s2_s+(seq1_len-align.s1_s)+1);
			memcpy(merge, seq2, sizeof(char) * align.s2_s);
			memcpy(&(merge[align.s2_s]), &(seq1[align.s1_s]), sizeof(char) * (seq1_len-align.s1_s));

			return merge;
		}

		/* If shorter one has higher supporting read density, choose shorter one.
		 * In this case, do not worry about min_len edge. */
		/*
		                   -------------=======================
				========================--
				              ^        ^
		*/
		MALLOC(merge, char, align.s2_e+1+(seq1_len-align.s1_e-1)+1);
		memcpy(merge, seq2, sizeof(char) * (align.s2_e+1));
		memcpy(&(merge[align.s2_e+1]), &(seq1[align.s1_e+1]), sizeof(char) * (seq1_len-align.s1_e-1));

		return merge;
	}
	else if(right_offset >= seq1_len){ /* You should adjust between `offset' and `length' */
		/* When longer one is better and they are not so similar to each other, avoid the edge part of shorter one. */
		if(align.s1_supp_read_den >= align.s2_supp_read_den){ /* Prefer longer one */
			/*
		     	=============================================--
		                             	---------------------========
		                                 	^               ^
			*/
			if(seq1_len-align.s1_e < min_len){
				align.s2_e -= (align.s1_e - (seq1_len-min_len));
				align.s1_e = seq1_len - min_len;
			}
			MALLOC(merge, char, align.s1_e+1+(seq2_len-align.s2_e-1)+1);
			memcpy(merge, seq1, sizeof(char) * (align.s1_e+1));
			memcpy(&(merge[align.s1_e+1]), &(seq2[align.s2_e+1]), sizeof(char) * (seq2_len-align.s2_e-1));

			return merge;
		}

		/* If shorter one has higher supporting read density, choose shorter one.
		 * In this case, do not worry about min_len edge. */
		/*
		     ============================-------------------
		                             ----=========================
		                                 ^               ^
		*/
		MALLOC(merge, char, align.s1_s+(seq2_len-align.s2_s)+1);
		memcpy(merge, seq1, sizeof(char) * (align.s1_s));
		memcpy(&(merge[align.s1_s]), &(seq2[align.s2_s]), sizeof(char) * (seq2_len-align.s2_s));

		return merge;
	}
	else{
		if(MIN(seq1_len, seq2_len) == 831 || MIN(seq1_len, seq2_len) == 833) printf("Error in Merge Align");
		/* Just choose longer one, if longer is better or shorter is exactly same with longer. */
		if((align.s1_supp_read_den > align.s2_supp_read_den) ||
					((double)align.matches/(double)MIN(seq1_len, seq2_len) == 1.0)) return strdup(seq1);
		if(MIN(seq1_len, seq2_len) == 831 || MIN(seq1_len, seq2_len) == 833) printf("Error in Merge Align 2");
		/*
		    ========------------------------------=========
		         ---==============================---
		            ^                            ^
		*/
		MALLOC(merge, char, align.s1_s+(align.s2_e-align.s2_s+1)+seq1_len-align.s1_e-1+1);
		memcpy(merge, seq1, sizeof(char) * (align.s1_s));
		memcpy(&(merge[align.s1_s]), &(seq2[align.s2_s]), sizeof(char) * (align.s2_e-align.s2_s+1));
		memcpy(&(merge[align.s1_s+(align.s2_e-align.s2_s+1)]), &(seq1[align.s1_e+1]), sizeof(char) * (seq1_len-align.s1_e-1));

		return merge;
	}
}

void printAlign(alignment align)
{
	printf("\n"
			"%u\t- sequence 1 length in the alignment\n"
			"%u\t- sequence 2 length in the alignment\n"
			"%u\t- sequence 1 align start offset\n"
			"%u\t- sequence 1 align end offset\n"
			"%u\t- sequence 2 align start offset\n"
			"%u\t- sequence 2 align end offset\n"
			"%u\t- number of matched bases\n"
			"%f\t- sequence 1 supporting read density\n"
			"%f\t- sequence 2 supporting read density\n"
			, align.s1_len, align.s2_len, align.s1_s, align.s1_e, align.s2_s, align.s2_e, align.matches,
			align.s1_supp_read_den, align.s2_supp_read_den);
}

void swapAlign(alignment *align)
{
	unsigned int temp;
	double dtemp;

	temp = align->s1_len;
	align->s1_len = align->s2_len;
	align->s2_len = temp;

	temp = align->s1_s;
	align->s1_s = align->s2_s;
	align->s2_s = temp;

	temp = align->s1_e;
	align->s1_e = align->s2_e;
	align->s2_e = temp;

	dtemp = align->s1_supp_read_den;
	align->s1_supp_read_den = align->s2_supp_read_den;
	align->s2_supp_read_den = dtemp;
}

char* itoa(int value, char* result, int base) 
{
	// check that the base if valid
	if (base < 2 || base > 36) { *result = '\0'; return result; }
	
	char* ptr = result, *ptr1 = result, tmp_char;
	int tmp_value;
	
	do {
		tmp_value = value;
		value /= base;
		*ptr++ = "zyxwvutsrqponmlkjihgfedcba9876543210123456789abcdefghijklmnopqrstuvwxyz" [35 + (tmp_value - value * base)];
	} while ( value );
	
	// Apply negative sign
	if (tmp_value < 0) *ptr++ = '-';
	    *ptr-- = '\0';
	while(ptr1 < ptr) {
		tmp_char = *ptr;
		*ptr--= *ptr1;
		*ptr1++ = tmp_char;
	}
	return result;
}
