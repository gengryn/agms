
#include <stdlib.h>
#include <stdio.h>
#include "agms.h"
#include "prng.h"
#include "massdal.h"

AGMS_type * AGMS_Init(int buckets, int depth)
{
    int i,j;
    AGMS_type * result;
    prng_type * prng;
    
    prng=prng_Init(-6321371,2);
    
    result=calloc(1,sizeof(AGMS_type));
    if (result==NULL) exit(1);
    result->depth=depth;
    result->buckets=buckets;
    result->count=0;
    for (i=0;i<6;i++)
        result->test[i]=calloc(depth,sizeof(long long));
    // create space for the hash functions
    
    result->counts=(int *) calloc(buckets*depth, sizeof(int));
    if (result->counts==NULL) exit(1);
    
    for (i=0;i<depth;i++)
    {
        for (j=0;j<6;j++)
        {
            result->test[j][i]=(long long) prng_int(prng);
            if (result->test[j][i]<0) result->test[j][i]= -result->test[j][i];
            // initialise the hash functions
            // prng_int() should return a random integer
            // uniformly distributed in the range 0..2^31
        }
    }
    prng_Destroy(prng);
    return (result);
}

void AGMS_Update(AGMS_type * agms, unsigned long item, int diff)
{
    // update the sketch
    // hash to one bucket in each row
    // then multiply by {+1, -1} chosen at random
    
    int j;
    unsigned int hash;
    int mult, offset;
    
    agms->count+=diff;
    offset=0;
    for (j=0;j<agms->depth;j++)
    {
        hash=hash31(agms->test[0][j],agms->test[1][j],item);
        hash=hash % (agms->buckets);
        mult=fourwise(agms->test[2][j],agms->test[3][j],
                      agms->test[4][j],agms->test[5][j],item);
        if ((mult&1)==1)
            agms->counts[offset+hash]+=diff;
        else
            agms->counts[offset+hash]-=diff;
        offset+=agms->buckets;
    }
}


int AGMS_Compatible(AGMS_type * a, AGMS_type * b){
    int i,j;
    // test whether two sketches have the same parameters
    // if so, then they can be added, subtracted, etc.
    
    
    if (!a || !b) return 0;
    if (a->buckets!=b->buckets) return 0;
    if (a->depth!=b->depth) return 0;
    for (i=0;i<a->depth;i++)
        for (j=0;j<6;j++)
            if (a->test[j][i]!=b->test[j][i])
                return 0;
    return 1;
}

int AGMS_Count(AGMS_type * agms, int item)
{
    // compute the estimated count of item
    
    int i;
    int offset;
    int * estimates;
    unsigned int hash;
    int mult;
    
    estimates=(int *) calloc(1+agms->depth, sizeof(int));
    offset=0;
    for (i=1;i<=agms->depth;i++)
    {
        hash=hash31(agms->test[0][i-1],agms->test[1][i-1],item);
        hash=hash % (agms->buckets);
        mult=fourwise(agms->test[2][i-1],agms->test[3][i-1],
                      agms->test[4][i-1],agms->test[5][i-1],item);
        if ((mult&1)==1)
            estimates[i]=agms->counts[offset+hash];
        else
            estimates[i]=-agms->counts[offset+hash];
        offset+=agms->buckets;
    }
    if (agms->depth==1) i=estimates[1];
    else if (agms->depth==2) i=(estimates[1]+estimates[2])/2;
    else
        i=MedSelect(1+agms->depth/2,agms->depth,estimates);
    free(estimates);
    return(i);
}

long long AGMS_F2Est(AGMS_type * agms)
{
    // estimate the F2 moment of the vector (sum of squares)
    
    int i,j, r;
    long long * estimates;
    long long result, z;
    
    estimates=(long long *) calloc(1+agms->depth, sizeof(long long));
    r=0;
    for (i=1;i<=agms->depth;i++)
    {
        z=0;
        for (j=0;j<agms->buckets;j++)
        {
            z+=((long long) agms->counts[r]* (long long) agms->counts[r]);
            r++;
        }
        estimates[i]=z;
    }
    if (agms->depth==1) result=estimates[1];
    else if (agms->depth==2) result=(estimates[1]+estimates[2])/2;
    else
        result=LLMedSelect(1+agms->depth/2,agms->depth,estimates);
    free(estimates);
    return(result);
}

long long AGMS_InnerProd(AGMS_type * a, AGMS_type * b){
    int i,j, r;
    long long * estimates;
    long long result, z;
    // estimate the innerproduct of two vectors using their sketches.
    
    if (AGMS_Compatible(a,b)==0) return 0;
    estimates=(long long *) calloc(1+a->depth, sizeof(long long));
    r=0;
    for (i=1;i<=a->depth;i++)
    {
        z=0;
        for (j=0;j<a->buckets;j++)
        {
            z+=((long long) a->counts[r]* (long long) b->counts[r]);
            r++;
        }
        estimates[i]=z;
    }
    if (a->depth==1) result=estimates[1];
    else if (a->depth==2) result=(estimates[1]+estimates[2])/2;
    else
        result=LLMedSelect(1+a->depth/2,a->depth,estimates);
    free(estimates);
    return(result);
    
}

int AGMS_AddOn(AGMS_type * dest, AGMS_type * source){
    int i,j,r;
    
    // add one sketch to another
    
    r=0;
    if (AGMS_Compatible(dest,source)==0) return 0;
    for (i=0;i<source->depth;i++)
        for (j=0;j<source->buckets;j++)
        {
            dest->counts[r]+=source->counts[r];
            r++;
        }
    return 1;
}

int AGMS_Subtract(AGMS_type * dest, AGMS_type * source){
    int i,j,r;
    
    // subtract one sketch from another
    
    r=0;
    if (AGMS_Compatible(dest,source)==0) return 0;
    for (i=0;i<source->depth;i++)
        for (j=0;j<source->buckets;j++)
        {
            dest->counts[r]-=source->counts[r];
            r++;
        }
    return 1;
}

int AGMS_Size(AGMS_type * agms){
    int size;
    
    // return the space used in bytes of the sketch
    
    size=(sizeof(int *))+(agms->buckets*agms->depth)*sizeof(int)+
    agms->depth*6*sizeof(long long)+sizeof(AGMS_type);
    return size;
}

void AGMS_Destroy(AGMS_type * agms)
{
    // destroy the data structure
    
    int i;
    
    if (agms)
    {
        for (i=0;i<6;i++)
            free(agms->test[i]);
        free(agms->counts);
        free(agms);
    }
}
