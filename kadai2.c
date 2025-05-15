#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BUFSIZE 1024 //ファイルから読み込む一行の最大文字数
#define MAX_SEQ_NUM 30 //一つの転写因子に対して与えられる結合部位配列の最大数
#define MAX_GENE_NUM 8 /*与えられるプロモータ領域の最大遺伝子数*/

char g_motif[MAX_SEQ_NUM][BUFSIZE]; //転写因子の結合部位配列を保存する配列

struct promoter{
  char name[BUFSIZE];
  char seq[BUFSIZE];
}g_pro[MAX_GENE_NUM]; //遺伝子のプロモータ領域を保存する構造体

//グローバル変数はローカル変数と区別するため、全部大文字にするかg_を先頭につけるのが一般的

int read_multi_seq(char* filename){
  int seq_num = 0;
  char buffer[BUFSIZE];
  FILE *fp = fopen("MATa1","r");

  if(fp == NULL){
    printf("motif_region_file open error.\n");
    exit(1); //ファイルが開けなかった場合プログラムを終了
  }

  while(fscanf(fp, "%s", buffer) != EOF){ //プログラムから一行ずつ読み込む
    if(buffer[strlen(buffer)-1]=='\n'){
      buffer[strlen(buffer)-1]='\0'; //改行を切り落とす
    }
    strcpy(g_motif[seq_num],buffer); //結合部位配列を保存
    seq_num++;
  }
  fclose(fp);
  return seq_num;
}

int read_promoter(char *filename){
  int gene_num = 0;  
  char buffer[BUFSIZE];
  FILE *fp = fopen("promoters","r");

  if(fp == NULL){
    printf("scorefile open error.\n");
    exit(1);
  }

  while(fscanf(fp, "%s", buffer) != EOF){
    if(buffer[strlen(buffer)-1]=='\n'){
      buffer[strlen(buffer)-1]='\0';
    }
    
    if(buffer[0]=='>'){
      strcpy(g_pro[gene_num].name,buffer+1); 
    }else{
      strcpy(g_pro[gene_num].seq,buffer);
      gene_num++;
    }    
  }
  fclose(fp);
  return gene_num;
}

int main(int argc, char* argv[]){
  int seq_num = read_multi_seq(argv[1]); //１番目の引数で指定した転写因子の複数の結合部位配列を読み込む

  printf("motif region:\n");
  for(int i = 0; i < seq_num; i++){
    printf("%s\n",g_motif[i]); //読み込んだ転写因子の結合部位配列を表示
  }
  printf("\n");

  int gene_num = read_promoter(argv[2]);  //２番目の引数で指定した遺伝子のプロモータ領域を読み込む
  
  printf("promoter_sequence:\n");
  for(int i = 0; i < gene_num; i++){
    printf(">%s\n", g_pro[i].name); //読み込んだプロモータ領域を表示
    printf("%s\n", g_pro[i].seq);
  }

 //頻度表の作成 
int k, l;
//配列の長さを取得
int num=strlen(g_motif[0]);
int freq[4][num];
//表の初期化
for(k=0; k<4; k++)
{
  for(l=0; l<num; l++)
  {
    freq[k][l]=0;
  }
}
//頻度を数える
for(k=0; k<seq_num; k++)
{
  for(l=0; l<num; l++)
  {
    if(g_motif[k][l]=='A')
    {
      freq[0][l]++;
    }
    if(g_motif[k][l]=='C')
    {
      freq[1][l]++;
    }
    if(g_motif[k][l]=='G')
    {
      freq[2][l]++;
    }
    if(g_motif[k][l]=='T')
    {
      freq[3][l]++;
    }
  }
}
//試しに出力
printf("%d\n",num);
for(k=0; k<4; k++)
{
  for(l=0; l<num; l++)
  {
    printf("%3d ",freq[k][l]);
  }
  printf("\n");
}


//対数オッズスコア行列の作成
float s_i[4][num];
for(k=0; k<4; k++)
{
  for(l=0; l<num; l++)
  {
    s_i[k][l]=log()
  }
}
  return 0;
}