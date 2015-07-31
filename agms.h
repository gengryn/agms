
typedef struct AGMS_type{
  int depth;
  int buckets;
  int count;
  int * counts;
  int *test[6];
} AGMS_type;

extern AGMS_type * AGMS_Init(int, int);
extern void AGMS_Update(AGMS_type *, unsigned long, int);
extern long long AGMS_F2Est(AGMS_type *);
extern long long AGMS_InnerProd(AGMS_type *, AGMS_type *);
extern int AGMS_Subtract(AGMS_type *, AGMS_type *);
extern int AGMS_AddOn(AGMS_type *, AGMS_type *);
extern void AGMS_Destroy(AGMS_type *);
extern int AGMS_Size(AGMS_type *);

