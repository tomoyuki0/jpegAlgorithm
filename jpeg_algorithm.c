#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "hahuman.h"

#define N 256
#define pi 3.14

//関数化
void DCT(unsigned char *f, double *D);

void quantizer(double *D, int *q);

void zigzag(int *q, int *z);

int ENCOD(int *z, int *encod);

void DECOD(int *gshift, int *decod);

void IZ(int *decod, int *Iz);

void DQ(int *Iz, double *Dq);

void IDCT(double *Dq, double *ID);

int main(void)
{	
	//入力
	unsigned char *f=(unsigned char*) malloc(N*N*sizeof(unsigned char));
	//出力
	double *g=(double*) malloc(N*N*sizeof(double));
	//DCT
	double *D=(double*) malloc(N*N*sizeof(double));
	//quantizer
	int *q=(int*) malloc(N*N*sizeof(int));
	//zigzag
	int *z=(int*) malloc(N*N*sizeof(int));
	//coding
	int *encod=(int*) malloc(N*N*sizeof(int));
	//shift
	unsigned char *shift=(unsigned char*) malloc(N*N*sizeof(unsigned char));
	//
	unsigned char *fshift=(unsigned char*) malloc(N*N*sizeof(unsigned char));
	//gyakushift
	int *gshift=(int*) malloc(N*N*sizeof(int));
	//decoding
	int *decod=(int*) malloc(N*N*sizeof(int));
	//IDCT
	double *ID=(double*) malloc(N*N*sizeof(double));
	//dequantizer
	double *Dq=(double*) malloc(N*N*sizeof(double));
	//inversezigzag
	int *Iz=(int*) malloc(N*N*sizeof(int));

	char filename1[100], filename2[100], filename3[100], filename4[100];

	snprintf(filename1, sizeof(filename1), "lenna_uchar_%d-%d.raw", N, N);
	FILE* fp1 = fopen(filename1,"rb");
	if(fp1 == NULL)
	{
		printf("faild to open file of \"%s\".\n", filename1);
		return 1;
	}
	fread(f, sizeof(double), N*N, fp1);
	fclose(fp1);

	//DCT呼び出し
	DCT(f,D);

	//量子化呼び出し
	quantizer(D,q);

	//ジグザグスキャン呼び出し
	zigzag(q,z);

	//符号化呼び出し
	ENCOD(z,encod);

	int i,j;
	//シフト演算
	int hugougo;
	hugougo=(ENCOD(z,encod)/8)+1;
	for(j=0;j<hugougo;j++){
		shift[j]=0;
		for(i=0;i<8;i++){
			if(encod[j*8+i]==1){
				shift[j]++;
				if(i<7){
					shift[j]=shift[j]<<1;
				}
			}
			if(encod[j*8+i]==0){
				if(i<7){
					shift[j]=shift[j]<<1;
				}
			}
		}
	}

	
	//ファイル書き込み
	snprintf(filename3, sizeof(filename3), "1.0_k_%d-%d.text", N, N);
	FILE* fp3 = fopen(filename3, "wb");
	if(fp3 == NULL){
		printf("faild to open file of \"%s\".\n", filename3);
		return 1;
	}
	fwrite(shift, sizeof(unsigned char), hugougo, fp3);
	fclose(fp3);

	printf("Output \"%s\"\n",filename3);

	
	//ファイル読み込み
	snprintf(filename4, sizeof(filename4), "1.0_k_%d-%d.text", N, N);
	FILE* fp4 = fopen(filename4,"rb");
	if(fp4 == NULL)
	{
		printf("faild to open file of \"%s\".\n", filename4);
		return 1;
	}
	fread(fshift, sizeof(unsigned char), hugougo, fp4);
	fclose(fp4);
	
	
	

	
	//シフト演算復号
	for(j=0;j<hugougo;j++){
		for(i=0;i<8;i++){
			gshift[j*8+(7-i)]=fshift[j]%2;
			if(i<7){
				fshift[j]=fshift[j]>>1;
			}
	    }
	}

	
	//復号化呼び出し
	DECOD(gshift,decod);
	
	//逆ジグザグスキャン呼び出し
	IZ(decod,Iz);

	//逆量子化呼び出し
	DQ(Iz,Dq);

	//IDCT呼び出し
	IDCT(Dq,ID);
	
	snprintf(filename2, sizeof(filename2), "1.0_shift_jpeg_lenna_new_uchar_%d-%d.raw", N, N);
	FILE* fp2 = fopen(filename2, "wb");
	if(fp2 == NULL){
		printf("faild to open file of \"%s\".\n", filename2);
		return 1;
	}
	fwrite(ID, sizeof(double), N*N, fp2);
	fclose(fp2);

	printf("Output \"%s\"\n",filename2);

	if(f == NULL)
	{
		printf("failed to allocate memory of f.\n");
		return 1;
	}

	free(f);
	
	if(g == NULL)
	{
		printf("failed to allocate memory of g.\n");
		return 1;
	}

	free(g);

	return 0;
}

//DCT関数
void DCT(unsigned char *f, double *D)
{
	int i,j,m,n,k,u,v;
	
	double cu,cv;
	double sum=0;
	cu=0;
	cv=0;

	for(m=0;m<N;m+=8){
		for(n=0;n<N;n+=8){
			for(u=0;u<8;u++){
				for(v=0;v<8;v++){
					for(j=0;j<8;j++){
						for(k=0;k<8;k++){
							if(u==0){
								cu=sqrt(1.0/2.0);
							}

							else{
								cu=1.0;
							}

							if(v==0){
								cv=sqrt(1.0/2.0);
							}

							else{
								cv=1.0;
							}

							sum+=f[(j+m)*N+k+n]*cos(((2.0*j+1.0)*u*pi)/16.0)*cos(((2.0*k+1.0)*v*pi)/16.0);
						}
					}
					D[(u+m)*N+v+n]=sum*cu*cv/4.0;
					sum=0;
				}
			}
		}
	}
}

//量子化関数
void quantizer(double *D, int *q)
{
	int m,n,u,v;
	for(m=0;m<N;m+=8){
		for(n=0;n<N;n+=8){
			for(u=0;u<8;u++){
				for(v=0;v<8;v++){

					q[(u+m)*N+v+n]=round(D[(u+m)*N+v+n]/(1.0*q_table[u][v]));
					//printf("%d",q[(u+m)*N+v+n]);
				}
			}
		}
	}
}



//ジグザグ関数
void zigzag(int *q, int *z)
{
	int m,n,j,k;
	int dir = 0;
	int h=0;
	int d=0;
	j=k=0;
	for(m=0;m<N;m+=8){
		for(n=0;n<N;n+=8){
			
			while(d<=64+h-1){
				z[d] = q[(j+m)*N+k+n];
				if(dir){
					++j;
					--k;
					if(k<0 || j>=8){
						++k;

						if(j>=8){
							--j;
							++k;
						}
						dir = 0;
					}
				} 
				else { //dir==0
					--j;
					++k;
					if (j<0 || k>=8) {
						++j;

						if (k>=8) {
							++j;
							--k;
						}
						dir = 1;
					}
				}
				d++;
			}
			h+=64;
			j=0;
			k=0;
		}
	}
}

//符号化関数
int ENCOD(int *z, int *encod){
	int *z2=(int*) malloc(N*N*sizeof(int));
	int *bn=(int*) malloc(N*N*sizeof(int));
	int *bn2=(int*) malloc(N*N*sizeof(int));
	int *DIFF=(int*) malloc(N*N*sizeof(int));
	int *DIFF2=(int*) malloc(N*N*sizeof(int));
	int i,u;
	int SSSS=0;
	int SSSS2=0;
	int R=0;
	int e;
	int pre=0;
	int cnt=0;
	for(i=0;i<N*N;i++){
		//DC変換
		//char dc_code_table[12][9]
		//char dc_length_table[12]
		if(i%64==0){
			R=0;///////////////////////////////////////////////////ACのカウントの初期化
			DIFF[i]=z[i]-pre;
			pre=z[i];	
			if(DIFF[i]==0){
				SSSS=0;
				for(e=0;e<dc_length_table[SSSS];e++){
					encod[cnt]=dc_code_table[SSSS][e];//符号語SSSS
					cnt++;
				}
			}

			else{
				SSSS=log2(abs(DIFF[i]))+1;//対数10進数
				for(e=0;e<dc_length_table[SSSS];e++){
					encod[cnt]=dc_code_table[SSSS][e];//符号語SSSS
					cnt++;	
				}

				if(DIFF[i]>0){
					for(u=0;0<DIFF[i];u++){
						bn[u]=DIFF[i]%2;
	        	    	DIFF[i]=DIFF[i]/2;
	       			}
	           		while(u>0){
		           		encod[cnt]=bn[--u];
		   				cnt++;
	       			}	
				}

				if(DIFF[i]<0){
					DIFF2[i]=abs(DIFF[i]);
					for(u=0;0<DIFF2[i];u++){
						bn[u]=DIFF2[i]%2;
		   			    DIFF2[i]=DIFF2[i]/2;
	           		}
	    	       	while(u>0){
		    	   		encod[cnt]=bn[--u];
		       			if(encod[cnt]==0){
                   	   		encod[cnt]=1;
               			}
               			else{
                   	    	encod[cnt]=0;
                   		}
		       			cnt++;
	       			}
				}
			}
		}

		//AC変換
		//char ac_code_table[177][16]
		//char ac_length_table[177]
		if(i%64!=0){
			if(z[i]==0){
				R++;
			}
	
			if(z[i]!=0){
				SSSS2=log2(abs(z[i]))+1;//対数10進数
				while(R>=16){
					R-=16;
					for(e=0;e<ac_length_table[176];e++){
						encod[cnt]=ac_code_table[176][e];
						cnt++;
					}
				}

				for(e=0;e<ac_length_table[R*11+SSSS2];e++){
					encod[cnt]=ac_code_table[R*11+SSSS2][e];
					cnt++;
				}

				R=0;

				if(z[i]>0){
					for(u=0;0<z[i];u++){
						bn2[u]=z[i]%2;
		       			z[i]=z[i]/2;
	           		}
           			while(u>0){
	        	    	encod[cnt]=bn2[--u];
	        			cnt++;
	        		}
				}

				if(z[i]<0){
					z2[i]=abs(z[i]);
					for(u=0;0<z2[i];u++){
						bn2[u]=z2[i]%2;
		    	   	    z2[i]=z2[i]/2;
	        	   	}
	           		while(u>0){
		        		encod[cnt]=bn2[--u];
		        		if(encod[cnt]==0){
                   	    	encod[cnt]=1;
                   		}
                   		else{
                       		encod[cnt]=0;
                   		}
		       			cnt++;
	       			}
				}
			}

			else if(i%64==63){
				//EOB
				for(e=0;e<ac_length_table[0];e++){
					encod[cnt]=ac_code_table[0][e];
					cnt++;
				}
			}
		}
	}
	return cnt;
}

//復号化関数
void DECOD(int *gshift, int *decod){
	int i,j,k,m,u,v;
	int E=0;
	int X=0;
	int cnt2=0;
	int pre2=0;
	int kari=0;
	int kari2=0;
	int SSSS3=0;
	int SSSS4=0;
	int ZR=0;
	int I=0;
	int Y=0;
	for(u=0;u<32*32;u++){
		//DC復号
		for(i=0;i<12;i++){
			for(j=0;gshift[E+j]==dc_code_table[i][j]&&(j<9);j++){
				if(j==8||dc_code_table[i][j+1]==-1){
					SSSS3=i;//iの数だけのの後のビット数を読み込んで10進数に変換する
					break;
				}
			}
		}

		E+=dc_length_table[SSSS3];
		//最初が０なら負、１なら正
		X=E;
		if(gshift[X]==1){
          	for(v=SSSS3-1;v>=0;v--){
				decod[cnt2]+=gshift[E]*pow(2,v);
				E++;
			}
		}

		if(gshift[X]==0){
          	for(v=SSSS3-1;v>=0;v--){
          		if(gshift[E]==0){
       	    		gshift[E]=1;
           		}
           		else{
          			gshift[E]=0;    
          		}
				kari+=gshift[E]*pow(2,v);
				E++;
			}
			decod[cnt2]=-kari;
			kari=0;
		}
		decod[cnt2]+=pre2;
		pre2=decod[cnt2];
		cnt2++;

		//AC復号
		k=0;
		while(k<63){
			for(i=0;i<177;i++){
				for(j=0;gshift[E+j]==ac_code_table[i][j]&&(j<16);j++){
					if(j==15||ac_code_table[i][j+1]==-1){
						I=i;
						SSSS4=I%11;//桁数...あまり
						ZR=I/11;//ゼロラン
						break;
					}
				}
			}

			E+=ac_length_table[I];///////////////////////
			//EOBじゃないとき復号する
			if(0<I&&I<176){
				Y=E;//最初の文字の配列番号の保存
				//ZRの格納
				for(i=0;i<ZR;i++){
					decod[cnt2]=0;
					cnt2++;
				}
				k+=ZR;

				//最初が０なら負、１なら正
				if(gshift[Y]==1){
          			for(v=SSSS4-1;v>=0;v--){//読み込み
						decod[cnt2]+=gshift[E]*pow(2,v);
						E++;
					}
					k++;
				}

				if(gshift[Y]==0){
          			for(v=SSSS4-1;v>=0;v--){
          				if(gshift[E]==0){
            	    		gshift[E]=1;
           				}
           				else{
          					gshift[E]=0;    
          				}
						kari2+=gshift[E]*pow(2,v);
						E++;
					}
					k++;
					decod[cnt2]=-kari2;
					kari2=0;
				}
				cnt2++;
			}

			//0が16回続いたとき
			else if(I==176){
				for(m=0;m<16;m++){
					decod[cnt2]=0;
					cnt2++;
				}
				k+=16;
			}
			//EOBのとき
			else if(I==0){
				while(k<63){//真の時回り続ける
					decod[cnt2]=0;
					cnt2++;
					k++;
				}
			}
		}
	}
}

//逆ジグザグ関数
void IZ(int *decod, int *Iz){
	int i,j,m,n,k;
	int dir = 0;
	int h=0;
	int d=0;
	j=k=0;
	for(m=0;m<N;m+=8){
		for(n=0;n<N;n+=8){
			//printf("\n");
			while(d<=64+h-1){
				Iz[(j+m)*N+k+n]=decod[d];	
				if(dir){
					++j;
					--k;
					if(k<0 || j>=8){
						++k;

						if(j>=8){
							--j;
							++k;
						}
						dir = 0;
					}
				} 
				else { //dir==0
					--j;
					++k;
					if (j<0 || k>=8) {
						++j;

						if (k>=8) {
							++j;
							--k;
						}
						dir = 1;
					}
				}
				d++;
			}
			h+=64;
			j=0;
			k=0;
		}
	}
}

//逆量子化関数
void DQ(int *Iz, double *Dq){
	int m,n,u,v;
	for(m=0;m<N;m+=8){
		for(n=0;n<N;n+=8){
			for(u=0;u<8;u++){
				for(v=0;v<8;v++){
					Dq[(u+m)*N+v+n]=round(Iz[(u+m)*N+v+n]*(1.0*q_table[u][v]));
				}
			}
		}
	}
}

//IDCT関数
void IDCT(double *Dq, double *ID){
	int j,k,m,n,u,v;
	double cu,cv;
	double sum=0;
	cu=0;
	cv=0;
	for(m=0;m<N;m+=8){
		for(n=0;n<N;n+=8){
			for(u=0;u<8;u++){
				for(v=0;v<8;v++){
					for(j=0;j<8;j++){
						for(k=0;k<8;k++){
							if(j==0){
								cu=sqrt(1.0/2.0);
							}

							else{
								cu=1.0;
							}

							if(k==0){
								cv=sqrt(1.0/2.0);
							}

							else{
								cv=1.0;
							}
					
							sum+=cu*cv*Dq[(j+m)*N+k+n]*cos(((2.0*u+1.0)*j*pi)/16.0)*cos(((2.0*v+1.0)*k*pi)/16.0);
						}
					}
					ID[(u+m)*N+v+n]=sum/4.0;
					sum=0;
				}
			}
		}
	}	
}

