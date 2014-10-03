#ifndef __ALIGNMENT_H__
#define __ALIGNMENT_H__

#include <time.h>
#include <errno.h>
#include <unistd.h> /* for getpagesize() */
#include "common.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct _alignment{
	unsigned int s1_len;
	unsigned int s2_len;
	unsigned int s1_s; /* sequence 1, might be longer */
	unsigned int s1_e;
	unsigned int s2_s; /* sequence 2 */
	unsigned int s2_e;
	unsigned int matches;
	double s1_supp_read_den; /* supporting read density */
	double s2_supp_read_den;
}alignment;
#define INIT_ALIGNMENT	(alignment){1, LARGE_INT, 0, 1, LARGE_INT, 0, 0, 0.0, 0.0}
#define LONG_ALIGN(a) ( ((a).s1_len > (a).s2_len)? (a).s1_len:(a).s2_len )
#define SHORT_ALIGN(a) ( ((a).s1_len > (a).s2_len)? (a).s2_len:(a).s1_len )
enum {LEFT_PROTRUDING, MIDDLE, RIGHT_PROTRUDING};

alignment NW_alignment(const char *sequence, const char *reference, alignment align, int matrix, int print);
alignment NW_alignmentEx(const char *sequence, const char *reference, alignment align, int matrix, int print);
alignment SW_alignmentEx(const char *sequence, const char *reference, alignment align, int matrix, int print);
int align_type(alignment align, char *seq1, char *seq2);
char* merge_align(alignment align, char *seq1, char *seq2, unsigned int min_len);
void printAlign(alignment align);
void swapAlign(alignment *align);
char* itoa(int, char*, int);

/* allocate memory */
#define MALLOC(a, b, c)	\
do{\
	(a) = (b*)malloc(sizeof(b) * (c));\
	memset((a), 0x00, sizeof(b) * (c));\
}while(0)
/* free memory */
#define FREE(a) \
{\
	if((a) != NULL){\
		free((a));\
		(a) = NULL;\
	}\
}
#define A	0x0
#define G	0x1
#define C	0x2
#define T	0x3

enum {DISABLE, ENABLE};
enum {FAIL, SUCCESS};
enum {NO, YES};
enum {FALSE, TRUE};
enum{EDNAFULL, SJ_MATRIX, NW_MATRIX};

typedef struct _errors{
	unsigned int TP; /* True Positive */
	unsigned int TN0; /* True Negative because of no-candidate (no possible assignment) */
	unsigned int TNx; /* True Negative becasue of many-candidates (ambiguous, give-up) */
	unsigned int FP; /* False Positive (False Alarm, Type I error) */
	unsigned int FN0; /* False Negative (Type II error)because of no-candidate (no possible assignment) */
	unsigned int FNx; /* False Negative (Type II error)because of many-candidates (ambiguous, give-up) */
}errors;

extern char *program;
extern int debug_level;
extern struct _node *N_table;
extern unsigned int t_size;
extern FILE *debugfp;
extern time_t gtstart, gtend; /* for TIME_INFORMATION */
extern struct _read_set set;
extern const double gidentity;

/* In this project, I used the upper case. */
#define b2num(b)	(((b) == 'A')? 0: (((b) == 'G')? 1: (((b) == 'C')? 2: 3)))
#define num2b(n)	(((n) == 0)? 'A': (((n) == 1)? 'G': (((n) == 2)? 'C': 'T')))
/* converting base-pair to their pair. */
#define bcompb(a)	((a)=='A')? 'T':((a)=='G')? 'C':((a)=='C')? 'G': 'A'
#define IS_VALID_BASE(cchh) ((cchh == 'A') || (cchh == 'G') || (cchh == 'C') || (cchh == 'T') || (cchh == 'a') || (cchh == 'g') || (cchh == 'c') || (cchh == 't') )
#define IS_EXTENDED_VALID_BASE(cchh) ((cchh == 'A') || (cchh == 'G') || (cchh == 'C') || (cchh == 'T') || \
        (cchh == 'U') || (cchh == 'R') || (cchh == 'Y') || (cchh == 'K') || \
        (cchh == 'M') || (cchh == 'S') || (cchh == 'W') || (cchh == 'B') || \
        (cchh == 'D') || (cchh == 'H') || (cchh == 'V') || (cchh == 'N'))

#ifdef __cplusplus
}
#endif
#endif//__ALIGNMENT_H__
