

typedef struct{
 int t[3];
 double r[3];
 int code;
 double sco;
} TBL;

bool SearchMAPfft(MRC *,MRC *,double);
bool SearchMAPfftMT(MRC *,MRC *,double,bool);
bool SearchMAPfftMT_OVCC(MRC *,MRC *,double,int,bool); //Overlap or CCC
double GetScore(MRC *, MRC *, int [3]);
double GetScore2(MRC *, MRC *, int [3], double [5]);

void PrintTbl(TBL *,int,MRC *,MRC *,double,double);
