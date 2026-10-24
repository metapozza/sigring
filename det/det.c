// Independent C port of the fused push/pop integer-ring signature update, for cross-arch/compiler determinism checks.
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#ifndef WORD
#define WORD 64
#endif
#if WORD==64
typedef uint64_t R;
#define NB 8
#else
typedef unsigned __int128 R;
#define NB 16
#endif
static int D,N; static size_t off[16]; static R *t,*pw,*pl; static R binom[16][16];
static void powers(const R*v,R*buf){
  for(int i=0;i<D;i++) buf[off[0]+i]=v[i];
  for(int k=2;k<=N;k++){ size_t pn=off[k-1]-off[k-2];
    for(size_t a=0;a<pn;a++) for(int b=0;b<D;b++) buf[off[k-1]+a*D+b]=buf[off[k-2]+a]*v[b]; }
}
static void ingest(const int64_t*inc,const int64_t*out){
  R r[16],l[16];
  for(int i=0;i<D;i++){ r[i]=(R)(int64_t)inc[i]; l[i]=(R)0-(R)(int64_t)out[i]; }
  powers(r,pw);
  for(int m=N;m>=1;m--){
    size_t len=off[m]-off[m-1];
    for(size_t i=0;i<len;i++) t[off[m-1]+i]+=pw[off[m-1]+i];
    for(int k=1;k<m;k++){ size_t la=off[k]-off[k-1], lb=off[m-k]-off[m-k-1]; R c=binom[m][k];
      for(size_t a=0;a<la;a++){ R cx=c*t[off[k-1]+a]; for(size_t b=0;b<lb;b++) t[off[m-1]+a*lb+b]+=cx*pw[off[m-k-1]+b]; } }
  }
  powers(l,pl);
  for(int m=N;m>=1;m--){
    size_t len=off[m]-off[m-1];
    for(size_t i=0;i<len;i++) t[off[m-1]+i]+=pl[off[m-1]+i];
    for(int j=1;j<m;j++){ size_t la=off[j]-off[j-1], lb=off[m-j]-off[m-j-1]; R c=binom[m][j];
      for(size_t a=0;a<la;a++){ R cx=c*pl[off[j-1]+a]; for(size_t b=0;b<lb;b++) t[off[m-1]+a*lb+b]+=cx*t[off[m-j-1]+b]; } }
  }
}
static uint64_t fnv(const R*x,size_t n){ uint64_t h=1469598103934665603ULL;
  for(size_t i=0;i<n;i++){ R v=x[i]; for(int b=0;b<NB;b++){ h^=(uint64_t)(uint8_t)(v&0xff); h*=1099511628211ULL; v>>=8; } } return h; }
int main(int argc,char**argv){
  D=atoi(argv[1]);N=atoi(argv[2]);int W=atoi(argv[3]);long T=atol(argv[4]);long every=atol(argv[5]);
  size_t sz=1,acc=0; off[0]=0; for(int k=1;k<=N;k++){ sz*=D; acc+=sz; off[k]=acc; }
  t=calloc(acc,sizeof(R));pw=calloc(acc,sizeof(R));pl=calloc(acc,sizeof(R));
  for(int i=0;i<=N;i++){ binom[i][0]=1; for(int j=1;j<=i;j++) binom[i][j]=binom[i-1][j-1]+(j<=i-1?binom[i-1][j]:0); }
  int64_t*ring=calloc((size_t)W*D,8); int h=0; uint64_t s=12345; int64_t inc[16],old[16]; const int64_t amp=1000;
  for(long tk=1;tk<=T;tk++){
    for(int c=0;c<D;c++){ s=s*6364136223846793005ULL+1442695040888963407ULL; inc[c]=(int64_t)((s>>33)%(uint64_t)(2*amp+1))-amp; }
    memcpy(old,&ring[(size_t)h*D],8*D); memcpy(&ring[(size_t)h*D],inc,8*D); h=(h+1)%W;
    ingest(inc,old);
    if(tk%every==0) printf("%d %d %d %d %ld %016llx\n",D,N,W,WORD,tk,(unsigned long long)fnv(t,acc));
  }
  return 0;
}
