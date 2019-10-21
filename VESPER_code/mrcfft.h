

typedef struct{
 int t[3];
 double r[3];
 int code;
 double sco;
} TBL;

bool SearchMAPfft(MRC *,MRC *,double);
bool SearchMAPfftMT(MRC *,MRC *,double);
bool SearchMAPfftMT_OVCC(MRC *,MRC *,double,int); //Overlap or CCC
double GetScore(MRC *, MRC *, int [3]);

void PrintTbl(TBL *,int,MRC *,MRC *,double,double);
