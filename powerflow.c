#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#define delta 0.00001

FILE *in,*out;
int node,branch;
int PQ,PU,balance,generator;
int times;


//动态分配
int *Z_num1,*Z_num2,*k_y,*N_num,*N_type;
float *r,*x,*y,*k;
//线路两端节点号+电阻+电抗+电纳+是否有变压器+变比
float **y_grid_a,**y_grid_b;//阻抗矩阵，a=实部，b=虚部
float *P,*Q;
float *U;//节点电压
float *U0a,*U0b;//节点初值
float *dP,*dQ,*dU;
float *Ua,*Ub;
float *change_a,*change_b;
float *sum1,*sum2;
float **jgrid;

int i,j,n,m;//临时变量

void enough(){
    float tmp[2];
    tmp[0]=dP[0];
    tmp[1]=dQ[0];
    tmp[2]=Ua[0]*Ua[0]+Ub[0]*Ub[0];
    for(i=0;i<node;i++)
    {
        if(dP[i]>tmp[0]) tmp[0]=dP[i];
        if(dQ[i]>tmp[1]) tmp[1]=dQ[i];
        if(Ua[i]*Ua[i]+Ub[i]*Ub[i]>tmp[2]) tmp[2]=Ua[i]*Ua[i]+Ub[i]*Ub[i];
    }
     if(tmp[0]<delta&&tmp[1]<delta&&tmp[2]<delta) exit(1);
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
    U=(float *)malloc(sizeof(float)*node);
    U0a=(float *)malloc(sizeof(float)*node);
    U0b=(float *)malloc(sizeof(float)*node);

    dQ=(float *)malloc(sizeof(float)*PQ);
    dP=(float *)malloc(sizeof(float)*PQ);
    dU=(float *)malloc(sizeof(float)*PQ);

    for(i=0;i<branch;i++)
        {
            fscanf(in,"%d %d %f %f %f %d %f",&Z_num1[i],&Z_num2[i],&r[i],&x[i],&y[i],&k_y[i],&k[i]);
            Z_num1[i]-=1;
            Z_num2[i]-=1;
        }

    for(i=0;i<node;i++)
        {
            fscanf(in,"%d %d %f %f %f %f %f",&N_num[i],&N_type[i],&P[i],&Q[i],&U[i],&U0a[i],&U0b[i]);//1=PQ 2=PU 3=Balance
            N_num[i]-=1;

        }

    /*for(i=0;i<branch;i++)
        printf("%d %d %f %f %f %d %f\n",Z_num1[i],Z_num2[i],r[i],x[i],y[i],k_y[i],k[i]);
    */

}

void delta1(){

    change_a=(float *)malloc(sizeof(float)*PQ);
    change_b=(float *)malloc(sizeof(float)*PQ);
    Ua=(float *)malloc(sizeof(float)*node);//每个节点的电压计算
    Ub=(float *)malloc(sizeof(float)*node);

    for(i=0;i<node;i++)
            {
                change_a[i]=0;
                change_b[i]=0;
                dP[i]=0;
                dQ[i]=0;
                dU[i]=0;
                Ua[i]=U0a[i];
                Ub[i]=U0b[i];
            }

    for(i=0;i<node;i++)
    {
        for(j=0;j<node;j++)
            {
                change_a[i]+=y_grid_a[i][j]*Ua[j]-y_grid_b[i][j]*Ub[j];
                change_b[i]+=y_grid_a[i][j]*Ub[j]+y_grid_b[i][j]*Ua[j];
            }
        if(N_type[i]==1)//PQ
        {
            dP[i]=P[i]-Ua[i]*change_a[i]-Ub[i]*change_b[i];
            dQ[i]=Q[i]-Ub[i]*change_a[i]+Ua[i]*change_b[i];
        }
        if(N_type[i]==2)//PU
        {
            dP[i]=P[i]-U0a[i]*change_a[i]-U0b[i]*change_b[i];
            dU[i]=U0a[i]*U0a[i]+U0b[i]*U0b[i]-Ua[i]*Ua[i]-Ub[i]*Ub[i];
        }
        if(N_type[i]==3)//Balance
        {}
    fprintf(out,"Node=%d\tdP=%.5f\tdQ=%.5f\tdU=%.5f\n",i+1,dP[i],dQ[i],dU[i]);
    }

    times=0;
    enough();
    free(change_a);
    free(change_b);
}

void jacobi(){
    jgrid=(float **)malloc(sizeof(float *)*2*(node-1));
    for(i=0;i<2*(node-1);i++)
        jgrid[i]=(float *)malloc(sizeof(float *)*2*(node-1));

    sum1=(float *)malloc(sizeof(float *)*node);
    sum2=(float *)malloc(sizeof(float *)*node);
    for(i=0;i<node;i++)
    {
        sum1[i]=0;sum2[i]=0;
    }

    for(i=0;i<2*(node-1);i+=2){//Pi
            for(j=0;j<2*(node-1);j+=2)
            {
                for(n=0;n<node;n++)
                    sum1[i]+=y_grid_a[i][n]*Ua[n]-y_grid_b[i][n]*Ub[n];
                if(i!=j) jgrid[i][j]=-y_grid_a[i][j]*Ua[i]-y_grid_b[i][j]*Ub[i];
                else jgrid[i][j]=-y_grid_a[i][j]*Ua[i]-y_grid_b[i][j]*Ub[i]-sum1[i];
            }
            for(j=1;j<2*(node-1);j+=2)
            {
                 for(n=0;n<node;n++)
                    sum2[i]+=y_grid_a[i][n]*Ub[n]-y_grid_b[i][n]*Ua[n];
                if(i!=j) jgrid[i][j]=-y_grid_a[i][j]*Ub[i]+y_grid_b[i][j]*Ua[i];
                else jgrid[i][j]=-y_grid_a[i][j]*Ua[i]-y_grid_b[i][j]*Ub[i]-sum2[i];
            }
        }

    for(i=1;i<2*PQ+1;i+=2){//Qi
            for(j=0;j<2*(node-1);j+=2)
            {
                /*for(n=0;n<node;n++)
                    sum2[i]+=y_grid_a[i][n]*Ub[n]-y_grid_b[i][n]*Ua[n];*/
                if(i!=j) jgrid[i][j]=-y_grid_a[i][j]*Ub[i]+y_grid_b[i][j]*Ua[i];
                else jgrid[i][j]=-y_grid_a[i][j]*Ub[i]+y_grid_b[i][j]*Ua[i]+sum2[i];
            }
            for(j=1;j<2*(node-1);j+=2)
            {
                 /*for(n=0;n<node;n++)
                    sum1[i]+=y_grid_a[i][n]*Ua[n]-y_grid_b[i][n]*Ub[n];*/
                if(i!=j) jgrid[i][j]=-y_grid_a[i][j]*Ua[i]-y_grid_b[i][j]*Ub[i];
                else jgrid[i][j]=y_grid_a[i][j]*Ua[i]+y_grid_b[i][j]*Ub[i]-sum1[i];
            }
        }

    for(i=2*PQ+2;i<2*(node-1);i++){//Ui
            for(j=0;j<2*(node-1);j+=2)
            {
                if(i!=j) jgrid[i][j]=0;
                else jgrid[i][j]=-2*Ua[i];
            }

            for(j=1;j<2*(node-1);j+=2)
            {
                if(i!=j) jgrid[i][j]=0;
                else jgrid[i][j]=-2*Ub[i];
            }
        }
    /*TEST ONLY*/
    printf("Jacobi:\n");
    for(i=0;i<2*(node-1);i++)
        for(j=0;j<2*(node-1);j++)
            printf("%.5f\t",jgrid[i][j]);
}

void freemem(){
    free(Z_num1);
    free(Z_num2);
    free(r);
    free(x);
    free(y);
    free(k_y);
    free(k);
    free(dP);
    free(dQ);
    free(dU);
    free(N_num);
    free(N_type);
    free(P);
    free(Q);
    free(U0a);
    free(U0b);
    free(Ua);
    free(Ua);
    free(jgrid);
}

int main(){
    if((in=fopen("Data.txt","r"))==NULL) printf("ERROR.\n");
    if((out=fopen("Output.txt","w"))==NULL) printf("ERROR.\n");

    read();
    ygrid();
    fprintf(out,"\nTimes=0\n");
    delta1();
    jacobi();
    recurse();

    freemem();
    return 0;
}
