#ifndef MINPDPLUGIN_H
#define MINPDPLUGIN_H

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <math.h>


#include <string>
using std::string;

class MinPDPlugin : public Plugin {
   public:
	   void input(std::string);;
	   void run();;
	   void output(std::string);;

   private:
	           FILE    *fp1;
        char    *cc, *filename, *outputfile1, *outputfile2;
        char    temp[LINELIMIT], temp2[MAXLENGHT], timechar[3];
        int             i, k, n,j,Fr_Count, align, res_no =0;
    Fasta       *seqs[NUMSEQS]; /* Declare an array of Fasta structures */
        double  **dist, ***dist_frags;
        double  GOP, GEP, match, mismatch, alpha, mutrate;
        unsigned int    maxdim;
        Mintaxa *minresults[NUMSEQS];
int BitCriteria(int Num, int Bit);;/* Is bit set in Num */
int TwoCrossovers(int Num, int MaxBit, int *Begin, int *End);;/* Is bit set in Num */
char** MatrixInit(int s, int dim1, int dim2);
double TN93distance(int pT, int pC, int pG, int pA, int p1, int p2, int gaps,int matches, int al_len, double alpha, int s, int t);;
int DistOnly(Fasta **seqs, double **dist, double ***dist_frags, int maxdim, int n, double alpha, int Fr_Count);;
int Align_Dist(Fasta **seqs, double **dist, double ***dist_frags, int maxdim, int n,  double GOP, double  GEP, double match, double mismatch,double alpha, int Fr_Count, char* file3);;
int SaveMinResults(Mintaxa **minresults, int *resno, int a, int c, double dist);
void OptionalRecPrintout(FILE *fp4,  RecRes *recRes,Fasta **seqs,int *Frag_Seq,FDisNode (*recSolutions)[3][2],int i, int l, int r, int p, int q, int f, int index, int last, int BKP, int Fr_Count,int align,int fr_size);
void FillArrayInOrder(int *arrayInOrder, double ***dist_frags, int  *Frag_Seq, int s, int k, int rec_idx);
int PickSeqofLargerAvgDis(double ***dist_frags, int k_idx, int f_idx, int min_idx,  int s, int Fr_Count, int k, int f);
void PrintFragments(FILE *fp3, int  *Frag_Seq, Fasta **seqs, double **dist, double ***dist_frags, int Fr_Count, int i, int s, int align );
BKPNode* Check4Recombination(FILE *fp3, FILE *fp4, int  *Frag_Seq, int  *MinCount, Fasta **seqs, double **dist, double ***dist_frags, int StartAncestorSearch, int s, int Fr_Count, double min_d, int min_idx, int fr_size, int align);
BKPNode* SaveResultsinBKPLinkList(int  *Frag_Seq, Fasta **seqs, double ***dist_frags, int Fr_Count,int i, int align, int s, int maxdim );
int GetMinDist(Mintaxa **minresults, Fasta **seqs, double **dist, double ***dist_frags, int maxdim, int n, char *File1, char *File2, double mutrate, int Fr_Count, int *resno, int align);
BKPNode * GetRecSolutionsII(FILE *fp4,  double ***dist_frags, Fasta **seqs, int  *Frag_Seq, int rec_idx, int s, int Fr_Count, int min_seq, int fr_size, double min_seq_dist, int align);
void BuildPartialNJTrees(Fasta **seqs, Mintaxa **minresults, int resno, double **dist, char *File1, int n);

void heapify_mintaxa ( Mintaxa **list , int newnode );
void heapsort_mintaxa ( Mintaxa **list, int last );

void heapify_Fasta ( Fasta **list , int newnode );
void heapsort_Fasta ( Fasta **list, int last );

void PrintTree(FILE *fv, TNode *node, Fasta **seqs);
void AddTipOrNode(int *w,int *p,int *max_w, int chidx, double brlen,TNode **NodeStorage, int Join);
void NJTree(TNode **NodeStorage, double **DistMatrix, int UBound, int OutgroupIn0, int *seqIds, int *max_w, int *w);

};

#endif
