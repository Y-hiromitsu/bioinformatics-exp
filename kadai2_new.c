#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define RANDMAX 24314410
#define BUFSIZE 1024 //ファイルから読み込む一行の最大文字数
#define MAX_SEQ_NUM 30 //一つの転写因子に対して与えられる結合部位配列の最大数
#define MAX_GENE_NUM 8 /*与えられるプロモータ領域の最大遺伝子数*/
#define N 4


char g_motif[MAX_SEQ_NUM][BUFSIZE]; //転写因子の結合部位配列を保存する配列
char motif_name[BUFSIZE];

struct promoter{
  char name[BUFSIZE];
  char seq[BUFSIZE];
}g_pro[MAX_GENE_NUM]; //遺伝子のプロモータ領域を保存する構造体

//グローバル変数はローカル変数と区別するため、全部大文字にするかg_を先頭につけるのが一般的

int read_multi_seq(char* filename){
  int seq_num = 0;
  char buffer[BUFSIZE];
  FILE *fp = fopen(filename,"r");

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
  FILE *fp = fopen(filename,"r");

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

void freqtable(int num, int seq_num, int freq[][BUFSIZE]);//頻度表の作成
void oddsscorematrix(int num, int freq[][BUFSIZE], double s_i[][BUFSIZE], int seq_num);//対数オッズスコア行列を作る
void PrintOddsscorematrix(int num, double s_i[][BUFSIZE]);//対数オッズスコア行列の出力
void searchbindingsites(int num, int num_pro, int gene_num, double s_i[][BUFSIZE], double threshold);//結合部位の探索

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

//配列の長さを取得
int num=strlen(g_motif[0]);

 //頻度表の作成 
int k, l, m, p, s, t;
int freq[N][BUFSIZE];
printf("frequency table\n");
freqtable(num, seq_num, freq);



//対数オッズスコア行列の作成
double s_i[N][BUFSIZE];
oddsscorematrix(num, freq, s_i, seq_num);
PrintOddsscorematrix(num, s_i);

//ゲノム配列上の結合部位の探索
int num_pro=strlen(g_pro[0].seq); //プロモーター配列の長さを取得
double threshold=5.0;
searchbindingsites(num, num_pro, gene_num, s_i, threshold);

  return 0;
}

//頻度表を作成する関数
void freqtable(int num, int seq_num, int freq[][BUFSIZE])
{
    int k, l;
    for(k=0; k<seq_num; k++)
    {
        for(l=0; l<num; l++)
        {
            freq[k][l]=0;
        }
    }
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
    //出力
    for(k=0; k<N; k++)
    {
        for(l=0; l<num; l++)
        {
            printf("%5d ",freq[k][l]);
        }
        printf("\n");
    }
}

//対数オッズスコア行列を作成する関数
void oddsscorematrix(int num, int freq[][BUFSIZE], double s_i[][BUFSIZE], int seq_num)
{
    double p_i[N][BUFSIZE];
    int k, l;
    double bg_total=7519429*2+4637676*2;
    double bg[N]={7519429, 4637676, 4637676, 7519429};
    double q[N]={bg[0]/bg_total, bg[1]/bg_total, bg[2]/bg_total, bg[3]/bg_total};
    double seq_num_pse=seq_num+N;
    for(k=0; k<N; k++)
    {
        for(l=0; l<num; l++)
        {
            p_i[k][l]=(freq[k][l]+1)/(seq_num_pse);
            s_i[k][l]=log(p_i[k][l]/q[k]);
        }
    }
}

//対数オッズスコア行列を出力する関数
void PrintOddsscorematrix(int num, double s_i[][BUFSIZE])
{
  int k, l;
    for(k=0; k<N; k++)
    {
        for(l=0; l<num; l++)
        {
            printf("%5.2lf ",s_i[k][l]);
        }
        printf("\n");
    }
}

//結合部位を探索する関数
void searchbindingsites(int num, int num_pro, int gene_num, double s_i[][BUFSIZE], double threshold)
{
    int k, l, x;
    double hit_gene[gene_num][BUFSIZE];
    for(k=0; k<gene_num; k++)
  {
    for(l=0; l<BUFSIZE; l++)
    {
      hit_gene[k][l]=0;
    }
  }
  //プロモーター配列上のヒット
  for (k=0; k<gene_num; k++)
  {
    printf("gene:%s\n",g_pro[k].name);
    int start;
    for(start=0; start<num_pro-num; start++)
    {
      for(l=0; l<num; l++)
      {
        if(g_pro[k].seq[start+l]=='A')
        {
          hit_gene[k][start]=hit_gene[k][start]+s_i[0][l];
        }
        else if(g_pro[k].seq[start+l]=='C')
        {
          hit_gene[k][start]=hit_gene[k][start]+s_i[1][l];
        }
        else if(g_pro[k].seq[start+l]=='G')
        {
          hit_gene[k][start]=hit_gene[k][start]+s_i[2][l];
        }
        else if(g_pro[k].seq[start+l]=='T')
        {
          hit_gene[k][start]=hit_gene[k][start]+s_i[3][l];
        } 
      }
      if(hit_gene[k][start]>=threshold)
    {
      printf("position:%d\n",start+1);
      printf("hit(");
      for(x=start; x<start+num; x++)
      {
        printf("%c",g_pro[k].seq[x]);
      }
      printf(")=%.2lf\n", hit_gene[k][start]);
      printf("\n");
    }
    }
  }
}
