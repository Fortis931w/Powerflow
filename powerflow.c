#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#define fault 0.00001


FILE *in,*out;
int node,branch;
int PQ,PU,balance,generator;
int times=0;


//动态分配
int *Z_num1,*Z_num2,*k_y,*N_num,*N_type;
float *r,*x,*y,*k;
//线路两端节点号+电阻+电抗+电纳+是否有变压器+变比
float **y_grid_a,**y_grid_b;//阻抗矩阵，a=实部，b=虚部
float *P,*Q;
float *U1,*U2;//节点电压
float *U0a,*U0b;//节点初值
float *dP,*dQ,*dU;
float *Ua,*Ub;
float *dUa,*dUb;
float *sum1,*sum2;
float **jgrid;
float **jgrid_inv;

int i,j,n,m;//临时变量

void copyfloat2Dimension(int s, int n, float **source, float **dest){
    int i,j;
    for(i=0;i<s;i++)
    {
        for(j=0;j<n;j++)
        {
            *(*(dest+i)+j)=*(*(source+i)+j);
        }
    }
}

double deter(int dimension,float **array){
    int i,j,k,l,b;
    int flag =1;
    float sum=1;
    float temp;
    for(i=0,j;i<dimension-1;i++)
    {
        j=i;
        if(*(*(array+i)+j)==0)
        {
            b=0;
            for(k=i+1;k<dimension;k++)
            {
                if(*(*(array+k)+j)!=0)//找到一行不为0的,然后换行
                {
                    for(l=j;l<dimension;l++)
                    {
                        temp=*(*(array+k)+l);
                        *(*(array+k)+l)= *(*(array+i)+l);
                        *(*(array+i)+l)=temp;
                    }
                    flag*=-1;
                    b=1;
                    break;
                }
            }
            if(!b)
            {
                return 0;
            }
            i--;
            continue;
        }
        for(;j<dimension-1;j++)
        {
            if(*(*(array+j+1)+i)==0)
                continue;
            temp = -*(*(array+j+1)+i)/ *(*(array+i)+i);
            for(k=i;k<dimension;k++)
                *(*(array+j+1)+k)+= *(*(array+i)+k)*temp;
        }
    }

    for(i=0;i<dimension;i++)
        sum*= *(*(array+i)+i);
    return sum*flag;
}

void getCompanionMatrix(int dimension, float **array, float **companionMatrix){
    int i,j,k,l,m,n,o;
    int flag;//标志代数余子式的符号
    float **companionTemp;
    //float deter(int dimension,float **array);

    companionTemp =(float**)malloc((dimension-1)*sizeof(float*));
    for(i=0;i<dimension-1;i++)
        companionTemp[i]=(float*)malloc((dimension-1)*sizeof(float));

    for(i=0;i<dimension;i++)
    {
        for(j=0;j<dimension;j++)
        {
            flag=(i+j)%2==0?1:-1;
            for(k=0,m=0;k<dimension;k++)
            {
                if(k==i)continue;
                for(l=0,n=0;l<dimension;l++)
                {
                    if(l==j)continue;
                    *(*(companionTemp+m)+n) = *(*(array+k)+l);
                    n++;
                }
                m++;
            }
            //第i行，第j列的代数余子式 赋值给第j行，第i列
            *(*(companionMatrix+j)+i) = flag * deter(dimension-1,companionTemp);
        }
    }
}

void inverse(float **array){
    float deterValue=1;
    float **deterArray, **companionMatrix;
    int dimension;
        dimension=2*(node-1);
    deterArray=(float**)malloc(dimension*sizeof(float*));
    companionMatrix =(float**)malloc(dimension*sizeof(float*));
    jgrid_inv=(float **)malloc(dimension*sizeof(float*));


        for(i=0;i<dimension;i++)
    {
        deterArray[i]=(float *)malloc(dimension*sizeof(float));
        companionMatrix[i]=(float *)malloc(dimension*sizeof(float));
        jgrid_inv[i]=(float *)malloc(dimension*sizeof(float));
    }
    copyfloat2Dimension(dimension,dimension,array,deterArray);

    deterValue = deter(dimension,deterArray);//

    getCompanionMatrix(dimension,array,companionMatrix);//

    for(i=0;i<dimension;i++)
        for(j=0;j<dimension;j++)
        jgrid_inv[i][j]=companionMatrix[i][j]/deterValue;

    printf("\nInverse Jacobian Matrix:\n");
    for(i=0;i<dimension;i++)
        {
            for(j=0;j<dimension;j++) printf("%- 6.6f ",jgrid_inv[i][j]);
            printf("\n");
        }


}

void ygrid(){
    y_grid_a=(float **)malloc(sizeof(float *)*node);
    for(i=0;i<node;i++)
        y_grid_a[i]=(float *)malloc(sizeof(float)*node);
    y_grid_b=(float **)malloc(sizeof(float *)*node);
    for(i=0;i<node;i++)
        y_grid_b[i]=(float *)malloc(sizeof(float)*node);
    //INITIATE
    for(i=0;i<node;i++)
        for(j=0;j<node;j++)
    {
        y_grid_a[i][j]=0.0;
        y_grid_b[i][j]=0.0;
    }
    //SELF
    for(i=0;i<node;i++)
    {
        for(j=0;j<branch;j++)//y*1/2&1/Z
            {
            if(Z_num1[j]==i) y_grid_b[i][i]+=y[j];
            if(Z_num2[j]==i) y_grid_b[i][i]+=y[j];
            if(Z_num1[j]==i)
                {
                if (k_y[j])
                    {
                        y_grid_a[i][i]+=r[j]/(x[j]*x[j]+r[j]*r[j])/(k[j]*k[j]);
                        y_grid_b[i][i]-=x[j]/(x[j]*x[j]+r[j]*r[j])/(k[j]*k[j]);
                    }
                else
                    {
                        y_grid_a[i][i]+=r[j]/(x[j]*x[j]+r[j]*r[j]);
                        y_grid_b[i][i]-=x[j]/(x[j]*x[j]+r[j]*r[j]);
                    }
                }
            if(Z_num2[j]==i)
                {
                    y_grid_a[i][i]+=r[j]/(x[j]*x[j]+r[j]*r[j]);
                    y_grid_b[i][i]-=x[j]/(x[j]*x[j]+r[j]*r[j]);
                }
            }
    }
    //TRANS
    for(i=0;i<node;i++)
    {
        for(j=i+1;j<node;j++)
            {
            for(n=0;n<branch;n++)
                {
                if(!k_y[n])
                    if((Z_num1[n]==i&&Z_num2[n]==j)||(Z_num2[n]==i&&Z_num1[n]==j))
                    {
                        y_grid_a[j][i]=y_grid_a[i][j]-=r[n]/(x[n]*x[n]+r[n]*r[n]);
                        y_grid_b[j][i]=y_grid_b[i][j]+=x[n]/(x[n]*x[n]+r[n]*r[n]);
                    }
                if(k_y[n]&&(Z_num1[n]==i&&Z_num2[n]==j))
                    {
                        y_grid_a[j][i]=y_grid_a[i][j]-=r[n]/(x[n]*x[n]+r[n]*r[n])/k[n];
                        y_grid_b[j][i]=y_grid_b[i][j]+=x[n]/(x[n]*x[n]+r[n]*r[n])/k[n];
                    }
                if(k_y[n]&&(Z_num2[n]==i&&Z_num1[n]==j))
                    {
                        y_grid_a[j][i]=y_grid_a[i][j]-=r[n]/(x[n]*x[n]+r[n]*r[n])*k[n];
                        y_grid_b[j][i]=y_grid_b[i][j]+=x[n]/(x[n]*x[n]+r[n]*r[n])*k[n];
                    }
                }
            }
    }

    /*For Test Only*/
    //G&&B
 /*
    fprintf(out,"\na=\n");
    for(i=0;i<node;i++)
        {for(j=0;j<node;j++) fprintf(out,"%.4f\t",y_grid_a[i][j]);
        fprintf(out,"\n");
    }
    fprintf(out,"\nb=\n");
        for(i=0;i<node;i++)
        {for(j=0;j<node;j++) fprintf(out,"%.4f\t",y_grid_b[i][j]);
        fprintf(out,"\n");
    }
    */

    fprintf(out,"\nY=\n");
    for(i=0;i<node;i++)
        {
            for(j=0;j<node;j++) fprintf(out,"%- 6.4f+j% -6.4f\t\t",y_grid_a[i][j],y_grid_b[i][j]);
            fprintf(out,"\n");
        }
}

void delta(){


    for(i=0;i<node;i++) {sum1[i]=0;sum2[i]=0;dP[i]=0;dQ[i]=0;dU[i]=0;}
    for(i=0;i<node;i++)//?node:PQ
        if(N_type[i]!=3)
            for(j=0;j<node;j++)
        {
            sum1[i]+=y_grid_a[i][j]*Ua[j]-y_grid_b[i][j]*Ub[j];
            sum2[i]+=y_grid_a[i][j]*Ub[j]+y_grid_b[i][j]*Ua[j];
        }

    for(i=0;i<node;i++)
        if(N_type[i]!=3)
        {
            dP[i]=P[i]-Ua[i]*sum1[i]-Ub[i]*sum2[i];
            if(N_type[i]==1) dQ[i]=Q[i]-Ub[i]*sum1[i]-Ua[i]*sum2[i];
            if(N_type[i]==2) dU[i]=U1[i]*U1[i]-U2[i]*U2[i]-Ua[i]*Ua[i]-Ub[i]*Ub[i];//U is squared
            printf("dP[%d]=%- 6.6f \t",i+1,dP[i]);
            printf("dQ[%d]=%- 6.6f \t",i+1,dQ[i]);
            printf("dU[%d]=%- 6.6f \n",i+1,dU[i]);
        }

    fprintf(out,"\n");
    fprintf(out,"Times=%4d\n",times);
    for(i=0;i<node;i++)
            fprintf(out,"Node%d\tP=%- 6.4f Q=% -6.4f U=% -6.4f+j% -6.4f\n",N_num[i]+1,P[i],Q[i],U1[i],U2[i]);
    fprintf(out,"\n");

    times++;
}

void jacobi(){
    jgrid=(float **)malloc(sizeof(float *)*2*(node-1));
    for(i=0;i<2*(node-1);i++)
        jgrid[i]=(float *)malloc(sizeof(float)*2*(node-1));

    sum1=(float *)malloc(sizeof(float)*node);
    sum2=(float *)malloc(sizeof(float)*node);

    for(i=0;i<2*(node-1);i++)
        {sum1[i]=0;sum2[i]=0;
        for(j=0;j<2*(node-1);j++)
            jgrid[i][j]=0;
        }

    for(i=0;i<node;i++)
        {
            if(N_type[i]==1)
            for(j=0;j<node;j++)
            {
                sum1[i]+=y_grid_a[i][j]*Ua[j]-y_grid_b[i][j]*Ub[j];
                sum2[i]+=y_grid_a[i][j]*Ub[j]+y_grid_b[i][j]*Ua[j];
            }
        //("sum1[%d]=% -6.4f\tsum2[%d]=% -6.4f\n",i,sum1[i],i,sum2[i]);
        }

    for(i=0;i<node-1;i++)
        {
            if(N_type[i]==1)
        {
            for(j=0;j<node-1;j++)
                    if(i==j)
                    {
                       jgrid[2*i][j*2]=-sum1[i]-y_grid_a[i][j]*Ua[j]-y_grid_b[i][j]*Ub[j];//pe
                       jgrid[2*i][j*2+1]=-sum2[i]-y_grid_a[i][j]*Ub[j]+y_grid_b[i][j]*Ua[j];//pf
                       jgrid[2*i+1][j*2]=sum2[i]+y_grid_b[i][j]*Ua[j]+y_grid_a[i][j]*Ub[j];//qe
                       jgrid[2*i+1][j*2+1]=-sum1[i]+y_grid_a[i][j]*Ua[j]+y_grid_b[i][j]*Ub[j];//qf
                    }

                //printf("% -6.4f\t% -6.4f\t% -6.4f\t% -6.4f\n",y_grid_a[i][j],Ua[j],y_grid_b[i][j],Ub[j]);
            else
            {
                jgrid[i*2][j*2]=-y_grid_a[i][j]*Ua[j]-y_grid_b[i][j]*Ub[j];//pe
                jgrid[2*i][j*2+1]=-y_grid_a[i][j]*Ub[j]+y_grid_b[i][j]*Ua[j];//pf
                jgrid[i*2+1][j*2+1]=-y_grid_a[i][j]*Ua[j]-y_grid_b[i][j]*Ub[j];//qf
                jgrid[2*i+1][j*2]=-y_grid_a[i][j]*Ub[j]+y_grid_b[i][j]*Ua[j];//qe
            }
        }
            if(N_type[i]==2)
                for(j=0;j<node-1;j++)
                if(i==j)
                   {
                       jgrid[2*i][j*2]=-sum1[i]-y_grid_a[i][j]*Ua[j]-y_grid_b[i][j]*Ub[j];//pe
                       jgrid[2*i][j*2+1]=-sum2[i]-y_grid_a[i][j]*Ub[j]+y_grid_b[i][j]*Ua[j];//pf
                       jgrid[2*i+1][j*2]=-2*Ua[i];//Ue
                       jgrid[2*i+1][j*2+1]=-2*Ub[i];//Uf
                    }
            else
            {
                jgrid[i*2][j*2]=-y_grid_a[i][j]*Ua[j]-y_grid_b[i][j]*Ub[j];//pe
                jgrid[2*i][j*2+1]=-y_grid_a[i][j]*Ub[j]+y_grid_b[i][j]*Ua[j];//pf
                jgrid[i*2+1][j*2+1]=0;//uf
                jgrid[2*i+1][j*2]=0;//ue
            }
        }

    /*TEST ONLY*/
    printf("Jacobi:\n");
    for(i=0;i<2*(node-1);i++)
    {
      for(j=0;j<2*(node-1);j++)
      printf("%- 6.6f ",jgrid[i][j]);
    printf("\n");
    }
    inverse(jgrid);

}

void multi(){

        /*Test only*/
    //int dimension=2*(n-1);


    dUa=(float *)malloc(node*sizeof(float));
    dUb=(float *)malloc(node*sizeof(float));
    for(i=0;i<node;i++) {dUa[i]=0;dUb[i]=0;}
    //printf("dU[%d]=%- 6.6f+j%- 6.6f \n",i+1,dUa[i],dUb[i]);

    for(i=0;i<node-1;i++)
        for(j=0;j<node-1;j++)
    {
        if(N_type[i]==1)
        {
            //printf("dP[%d]=%- 6.6f \tdQ[%d]=%- 6.6f \tdU[%d]=%- 6.6f \n",j+1,dQ[i],j+1,dP[j],j+1,dU[j]);
            //printf("%- 6.6f \t%- 6.6f \n%- 6.6f \t%- 6.6f \n",jgrid_inv[2*i][2*j],jgrid_inv[2*i][2*j+1],jgrid_inv[2*i+1][2*j],jgrid_inv[2*i+1][2*j+1]);
            dUa[i]+=-jgrid_inv[2*i][2*j]*dP[j];
            dUa[i]+=-jgrid_inv[2*i][2*j+1]*dQ[j];
            dUb[i]+=-jgrid_inv[2*i+1][2*j]*dP[j];
            dUb[i]+=-jgrid_inv[2*i+1][2*j+1]*dQ[j];
        }
        //printf("%- 6.6f \t",jgrid_inv[i][j]);

        if(N_type[i]==2)
        {
            dUa[i]+=-jgrid_inv[2*i][2*j]*dP[j];
            dUa[i]+=-jgrid_inv[2*i][2*j+1]*dU[j];
            dUb[i]+=-jgrid_inv[2*i+1][2*j]*dP[j];
            dUb[i]+=-jgrid_inv[2*i+1][2*j+1]*dU[j];
        }
    }
    for(i=0;i<node-1;i++)
    {
        U1[i]=Ua[i]-dUa[i];
        U2[i]=Ub[i]-dUb[i];
        printf("\nU[%d]=%- 6.6f+j%- 6.6f \n",i+1,U1[i],U2[i]);
    }

    free(dUa);free(dUb);
}

void read(){
    //fscanf(in,"%d %d %d %d %d",&node,&branch,&);
    fscanf(in,"%d %d %d %d %d %d",&node,&branch,&PQ,&PU,&balance,&generator);
    /*线路阻抗*/
    Z_num1=(int *)malloc(sizeof(int)*branch);
    Z_num2=(int *)malloc(sizeof(int)*branch);
    r=(float *)malloc(sizeof(float)*branch);
    x=(float *)malloc(sizeof(float)*branch);
    y=(float *)malloc(sizeof(float)*branch);
    k_y=(int *)malloc(sizeof(int)*branch);
    k=(float *)malloc(sizeof(float)*branch);

    N_num=(int *)malloc(sizeof(int)*node);
    N_type=(int *)malloc(sizeof(int)*node);
    Q=(float *)malloc(sizeof(float)*node);
    P=(float *)malloc(sizeof(float)*node);
    U1=(float *)malloc(sizeof(float)*node);
    U2=(float *)malloc(sizeof(float)*node);
    U0a=(float *)malloc(sizeof(float)*node);
    U0b=(float *)malloc(sizeof(float)*node);

    dQ=(float *)malloc(sizeof(float)*node);//
    dP=(float *)malloc(sizeof(float)*node);
    dU=(float *)malloc(sizeof(float)*node);

    for(i=0;i<branch;i++)
        {
            fscanf(in,"%d %d %f %f %f %d %f",&Z_num1[i],&Z_num2[i],&r[i],&x[i],&y[i],&k_y[i],&k[i]);
            Z_num1[i]-=1;
            Z_num2[i]-=1;
        }

    for(i=0;i<node;i++)
        {
            fscanf(in,"%d %d %f %f %f %f %f",&N_num[i],&N_type[i],&P[i],&Q[i],&U1[i],&U0a[i],&U0b[i]);//1=PQ 2=PU 3=Balance
            N_num[i]-=1;

        }


    sum1=(float *)malloc(sizeof(float)*node);
    sum2=(float *)malloc(sizeof(float)*node);
    Ua=(float *)malloc(sizeof(float)*node);//每个节点的电压计算
    Ub=(float *)malloc(sizeof(float)*node);

    for(i=0;i<node;i++)
        {
            Ub[i]=U0b[i];Ua[i]=U0a[i];U2[i]=0;
        }



}

int main(){
    if((in=fopen("Data.txt","r"))==NULL) printf("ERROR.\n");
    if((out=fopen("Output.txt","w"))==NULL) printf("ERROR.\n");

    read();
    ygrid();
    again:delta();
    jacobi();
    multi();

    float temp[3]={dP[0],dQ[0],dU[0]};
        for(i=0;i<node;i++)
    {
        if(dP[i]>temp[0]) temp[0]=dP[i];
        if(dQ[i]>temp[1]) temp[1]=dQ[i];
        if(dU[i]>temp[2]) temp[2]=dU[i];
    }
    if(temp[0]<fault&&temp[1]<fault&&sqrt(temp[2])<fault) return(0);
    else goto again;
}
