#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <conio.h>


int main()
{
    int i,j,k,k1,k2,w,p,g,r,l,t,li,co,h,elem[50][50],comu,col_est,lin_est;
    int no1,no2,x1,y1,x2,y2,z1,z2,coord[50][50],L[50][50];
    int qs,qm,qn,qe,keq;
    int modulo_lu,modulo_lv,lx,ly,lz,ux,uy,uz,vx,vy,vz,wx,wy,wz,lv,uu,vv,ww,prod_uv,prod_uw,prod_vw;
    int neq,mapano[50][6],mapael[50][12],apoio[50][50];
    double ni,nf,ea,in,ar,soma1,soma2,s;
    double carga[50][50],x[50],y[50],z[50],sec[50][10],cosseno,seno;
    double mat[50][10],rot[50][50],rig[50][50],eu,au,iu,ju,gu,ang,lu,area[50],pm[50][50],meg[50][12][12],mre[200][200],vet_for[50],u[50],iuz,iuy;
    double veint[50],defor[50],v_temp[50];
    double yi[50];
    double f2v, f8v, f6w, f12w, f3w, f9w, f5v, f11v, q1v, q2v, q1w, q2w;
    int secao_tipo;

char filename[80];
FILE *filein,*fileout;

  /*   Abertura de arquivo     */
  /*   sem argumentos na main   */
  /*   arquivo de entrada      */

  printf("arquivo de entrada:  ");
  /*scanf("%s",filename); */

  strcpy(filename,"teste1.txt");

  if((filein = fopen(filename,"rt")) == NULL)
  {
    printf("erro na abertura de arquivo:  %s\n",filename);
    system("PAUSE");
    exit(1);
  }
  else
  printf("abertura de %s OK\n",filename);

  /*   arquivo de saída   */

  printf("arquivo de saida:  ");
  /*scanf("%s",filename);*/

  strcpy(filename,"teste1_saida.txt");

  if((fileout = fopen(filename,"wt")) == NULL)
  {
    printf("erro na abertura de arquivo:  %s\n",filename);
    system("PAUSE");
    exit(1);
  }
  else printf("abertura de %s OK\n",filename);

  /* leitura no arquivo */
  /* declaração e leitura dos valores */

    printf("Digite a quantidade de nos:\n");
    fscanf(filein,"%d",&qn);
    fprintf(fileout,"qn = %d\n",qn);
    printf("%d\n",qn);

    printf("Informe as coordenadas dos nos:\n");
    for(i=0;i<qn;i++)
    {
       printf("coordenadas X do no %d:",i+1);
       fscanf(filein,"%lf",&x[i]);

       fprintf(fileout,"x[%d] = %f\n",i+1,x[i]);
       printf("x[i] = %f\n",x[i]);

       printf("coordenadas Y do no %d:",i+1);
       fscanf(filein,"%lf",&y[i]);

       fprintf(fileout,"y[%d] = %f\n",i+1,y[i]);
       printf("y[i] = %f\n",y[i]);

       printf("coordenadas Z do no %d:",i+1);
       fscanf(filein,"%lf",&z[i]);

       fprintf(fileout,"z[%d] = %f\n",i+1,z[i]);
       printf("z[i] = %f\n",z[i]);
    }

       printf("Informar o tipo de apoio: \n");

       for(j=0; j<qn; j++)
       {
         for(l=0; l<6; l++)
         {
          printf("Qual tipo de apoio do no %i do grau de liberdade %i:\n",j+1,l+1);
          fscanf(filein,"%d",&apoio[j][l]);
          fprintf(fileout,"apoio[%d][%d] = %d\n",j,l,apoio[j][l]);
          printf("apoio[%d][%d] = %d\n",j,l,apoio[j][l]);

          /* Digite o algarismo '0' para grau de liberdade livre e o algarismo '1' para grau de liberdade impedido */
         }
       }

     printf("determinar a quantidade de elementos:\n");
     fscanf(filein,"%d",&qe);
     fprintf(fileout,"qe = %d\n",qe);
     printf("%d\n",qe);

     /* Indexar a posição, as propriedades geometricas e de material a cada elemento de barra */
     for(k=0;k<qe;k++)
      {
       printf("No inicial do elemento %d:\n",k+1);
       fscanf(filein,"%d",&elem[k][0]);
       fprintf(fileout,"elem[k][0] = %d\n",elem[k][0]);
       printf("%d\n",elem[k][0]);
       elem[k][0]-=1;

       printf("No final do elemento %d:\n",k+1);
       fscanf(filein,"%d",&elem[k][1]);
       fprintf(fileout,"elem[k][1] = %d\n",elem[k][1]);
       printf("%d\n",elem[k][1]);
       elem[k][1]-=1;

       printf("Numero de identificacao do material do elemento (E) %d:\n",k+1);
       fscanf(filein,"%d",&elem[k][2]);
       fprintf(fileout,"elem[k][2] = %d\n",elem[k][2]);
       printf("%d\n",elem[k][2]);
       elem[k][2]-=1;

       printf("Numero de identificacao da secao do elemento %d:\n",k+1);
       fscanf(filein,"%d",&elem[k][3]);
       fprintf(fileout,"elem[k][3] = %d\n",elem[k][3]);
       printf("%d\n",elem[k][3]);
       elem[k][3]-=1;
      }

      printf("Impressao dos elementos\n");
      for(k=0;k<qe;k++)
       {
        printf("No inicial do elemento %d: %d \n",k+1,elem[k][0]+1);
        printf("No final   do elemento %d: %d \n",k+1,elem[k][1]+1);
        printf("Numero do material do elemento %d: %d \n",k+1,elem[k][2]+1);
        printf("Numero da secao do elemento %d: %d \n",k+1,elem[k][3]+1);
       }

        printf("Ler quantidade de materiais:\n");
        fscanf(filein,"%d",&qm);
        fprintf(fileout,"qm = %d\n",qm);
        printf("qm = %d\n",qm);

       /* Leitura das propriedades dos materiais */

       for(k=0;k<qm;k++)
        {
         printf("Modulo de elasticidade (E) do material (Pa) %i:\n",k+1);   /* E  */
         fscanf(filein,"%lf",&mat[k][0]);            /* vetor double mat*/
         fprintf(fileout,"mat[k] = %lf\n",mat[k][0]);
         printf("%lf\n",mat[k][0]);

         printf("Modulo de elasticidade ao cisalhamento (G) do material (Pa) %i: \n",k+1);    /* G  */
         fscanf(filein,"%lf",&mat[k][1]);            /* vetor double mat*/
         fprintf(fileout,"mat[k] = %lf\n",mat[k][1]);
         printf("%lf\n",mat[k][1]);
        }

         printf("Ler quantidade de secoes:\n");
         fscanf(filein,"%d",&qs);
         fprintf(fileout,"qs = %d\n",qs);
         printf("qs = %d\n",qs);

        /* Leitura das propriedades geometricas das seções */
       for(k=0;k<qs;k++)
        {

          /* "0" para secao retangular "1" para secao circular */
         fprintf(fileout,"Informe o tipo de secao = %d\n",secao_tipo);
         printf("Informe o tipo de secao = %d\n",secao_tipo);
         fscanf(filein,"%d",&secao_tipo);



         if (secao_tipo==0)
            {
         fscanf(filein,"%lf",&sec[k][4]);
         fprintf(fileout,"Inserir a dimensaão B da secao retangular (Metros) = %f\n",sec[k][4]);
         printf("Inserir a dimensaão B da secaao retangular (Metros) = %f\n",sec[k][4]);

         fscanf(filein,"%lf",&sec[k][5]);
         fprintf(fileout,"Inserir a dimensaão H da secaao retangular (Metros) = %f\n",sec[k][5]);
         printf("Inserir a dimensaão H da secaao retangular (Metros) = %f\n",sec[k][5]);

         sec[k][0]=(sec[k][4])*(sec[k][5]);

         sec[k][1]=((sec[k][4])*pow(sec[k][5],3))/12;

         sec[k][2]=((sec[k][5])*pow(sec[k][4],3))/12;

         sec[k][3]=sec[k][1]+sec[k][2];

         printf("Area da secao %i: \n",k+1);     /* A  */
         /*fscanf(filein,"%lf",&sec[k][0]); /* matriz double sec*/
         fprintf(fileout,"sec[k][0] = %f\n",sec[k][0]);
         printf("sec[k][0] = %f\n",sec[k][0]);

         printf("Inercia da secao em y %i: \n",k+1); /* Iy  */
         /*fscanf(filein,"%lf",&sec[k][1]); */
         fprintf(fileout,"sec[k][1] = %lf\n",sec[k][1]);
         printf("sec[k][1] = %lf\n",sec[k][1]);

         printf("Inercia da secao em z %i: \n",k+1); /* Iz  */
         /*fscanf(filein,"%lf",&sec[k][2]); */
         fprintf(fileout,"sec[k][2] = %lf\n",sec[k][2]);
         printf("sec[k][2] = %lf\n",sec[k][2]);

         printf("Momento de inercia polar %i:\n",k+1);    /* J  */
         /*fscanf(filein,"%lf",&sec[k][3]); /* vetor double mat*/
         fprintf(fileout,"sec[k] = %lf\n",sec[k][3]);
         printf("sec[k][3] = %lf\n",sec[k][3]);
            }

           if (secao_tipo==1)
            {

         fscanf(filein,"%lf",&sec[k][4]);
         fprintf(fileout,"Inserir o diametro (D) da secao circular (Metros) = %f\n",sec[k][4]);
         printf("Inserir o diametro (D) da secao circular (Metros) = %f\n",sec[k][4]);


         sec[k][0]=(pow(sec[k][4],2)*3.14159265358979323846)/4;

         sec[k][1]=(pow((sec[k][4]/2),4)*3.14159265358979323846)/4;

         sec[k][2]=(pow((sec[k][4]/2),4)*3.14159265358979323846)/4;

         sec[k][3]=sec[k][1]+sec[k][2];


         printf("Area da secao %i: \n",k+1);     /* A  */
         /*fscanf(filein,"%lf",&sec[k][0]); /* matriz double sec*/
         fprintf(fileout,"sec[k][0] = %f\n",sec[k][0]);
         printf("sec[k][0] = %f\n",sec[k][0]);

         printf("Inercia da secao em y %i: \n",k+1); /* Iy  */
         /*fscanf(filein,"%lf",&sec[k][1]); */
         fprintf(fileout,"sec[k][1] = %lf\n",sec[k][1]);
         printf("sec[k][1] = %lf\n",sec[k][1]);

         printf("Inercia da secao em z %i: \n",k+1); /* Iz  */
         /*fscanf(filein,"%lf",&sec[k][2]);*/
         fprintf(fileout,"sec[k][2] = %lf\n",sec[k][2]);
         printf("sec[k][2] = %lf\n",sec[k][2]);

         printf("Momento de inercia polar %i:\n",k+1);    /* J  */
         /*fscanf(filein,"%lf",&sec[k][3]); /* vetor double mat*/
         fprintf(fileout,"sec[k] = %lf\n",sec[k][3]);
         printf("sec[k][3] = %lf\n",sec[k][3]);
           }
        }


              /* Leitura das cargas externas atuantes sobre a estrutura */

       printf("defina as cargas pontuais externas atuantes na estrutura:\n");
      fprintf(fileout,"carga nos nos da  estrutura \n");
      for(w=0;w<qn;w++)
       {
         for(k=0;k<6;k++)
         {

             /* GL1 -> Fx
                GL2 -> Fy
                GL3 -> Fz
                GL4 -> Mx
                GL5 -> My
                GL6 -> Mz */

          fscanf(filein,"%lf",&carga[w][k]);
          fprintf(fileout,"carga[%d][%d] = %lf\n",w+1,k+1,carga[w][k]);
          printf("carga[%d][%d] = %lf\n",w+1,k+1,carga[w][k]);
         }
       }


        /* Leitura das cargas externas atuantes sobre a estrutura */

      fprintf(fileout,"carga distribuida nas barras da estrutura \n");

      for(w=0;w<qe;w++)
       {

          fscanf(filein,"%lf",&q1v);

          fscanf(filein,"%lf",&q2v);

          fscanf(filein,"%lf",&q1w);

          fscanf(filein,"%lf",&q2w);

          fprintf(fileout,"carga distribuida[%d] = %lf,%lf,%lf,%lf\n",w+1,q1v,q2v,q1w,q2w);
          printf("carga distribuida[%d] = %lf,%lf,%lf,%lf\n",w+1,q1v,q2v,q1w,q2w);

        no1= elem[w][0];
        no2= elem[w][1];
        x1=x[no1]; y1=y[no1];  z1=z[no1];
        x2=x[no2]; y2=y[no2];  z2=z[no2];
        lu=sqrt(pow(x2-x1,2) + pow(y2-y1,2) + pow(z2-z1,2));

        f2v = -(3*q1v+2*q2v)*(lu/10);

        f8v = -(2*q1v+3*q2v)*(lu/10);

        f6w = -((q1v/20)+(q2v/30))*pow(lu,2);

        f12w =((q1v/30)+(q2v/20))*pow(lu,2);


        f3w = -(3*q1w+2*q2w)*(lu/10);

        f9w = -(2*q1w+3*q2w)*(lu/10);

        f5v = -((q1w/20)+(q2w/30))*pow(lu,2);

        f11v =((q1w/30)+(q2w/20))*pow(lu,2);


        carga[no1][5]+=f2v;
        carga[no2][5]+=f8v;

        carga[no1][5]+=f5v;
        carga[no2][5]+=f11v;

        carga[no1][1]+=f3w;
        carga[no2][1]+=f9w;

        carga[no1][1]+=f6w;
        carga[no2][1]+=f12w;

       }

       printf("calculos: \n");
       printf("\n");
       fprintf(fileout,"\n");

       /* Montagem do mapa de equacoes por no */

       printf("Mapa de equacoes por no \n ");
       fprintf(fileout,"Mapa de equacoes por no \n ");

       neq=0;
       for(i=0;i<qn;i++)
       for(j=0;j<6;j++)
         {
           mapano[i][j]=0;
           if(apoio[i][j]==0)
           {
              neq++;
              mapano[i][j]=neq;
              printf("mapano[%d][%d] %d  ",i+1,j+1,mapano[i][j]);
              fprintf(fileout,"mapano[%d][%d] %d ",i+1,j+1,mapano[i][j]);
           }
         }

        printf("\n ");
        fprintf(fileout,"\n");

        /* Montagem do mapa de equacoes por elemento */
        printf("\n Mapa de equacoes por elemento");
        fprintf(fileout,"\n Mapa de equacoes por elemento");


    for(i=0;i<qe;i++)
    {
         for(j=0;j<12;j++)
         mapael[i][j]=0;
         no1=elem[i][0];
         no2=elem[i][1];

         for(k=0;k<6;k++)
          {
              if(mapano[no1][k]!=0)
              mapael[i][k]=mapano[no1][k];

              if(mapano[no2][k]!=0)
              mapael[i][k+6]=mapano[no2][k];
          }

          printf("\n ");
          fprintf(fileout,"\n");

         for(k=0;k<12;k++)
          {
           printf("mapael[%d][%d] %d ",i+1,k+1,mapael[i][k]);
           fprintf(fileout,"mapael[%d][%d] %d ",i+1,k+1,mapael[i][k]);
          }

    }
           printf("\n ");
           fprintf(fileout,"\n");

     printf("\n\n Montagem das matrizes de rigidez;\n");
     printf("\n ");
     fprintf(fileout,"\n");

    for(r=0;r<qe;r++)
     {
      for(li=0;li<12;li++)
       {
        for(co=0;co<12;co++)
          {
           rig[li][co]=0;
           rot[li][co]=0;
          }
       }

     }

        for(li=0;li<neq;li++)
            {
             for(co=0;co<neq;co++)
              {
               mre[li][co]=0;
              }
            }

   for(r=0;r<qe;r++)
    {
      no1= elem[r][0];
      no2= elem[r][1];
      x1=x[no1]; y1=y[no1];  z1=z[no1];
      x2=x[no2]; y2=y[no2];  z2=z[no2];

      lx=x2-x1;
      ly=y2-y1;
      lz=z2-z1;

          if((x2-x1==0)&&(y2-y1==0))
       {

     {
       if(z2>z1)
        {
         ux=0;
         uy=0;
         uz=1;

         vx=1;
         vy=0;
         vz=0;

         wx=0;
         wy=1;
         wz=0;
        }

        else
       {
        ux=0;
        uy=0;
        uz=-1;

        vx=0;
        vy=1;
        vz=0;

        wx=1;
        wy=0;
        wz=0;
       }

      }
       }

          else
        {

      modulo_lu=sqrt(pow(lx,2)+pow(ly,2)+pow(lz,2));

      ux=(lx)/modulo_lu;
      uy=(ly)/modulo_lu;
      uz=(lz)/modulo_lu;

      modulo_lv=sqrt(pow(-uz*ux,2)+pow(-uz*uy,2)+pow((1-uz*uz),2));

      vx=(-uz*ux)/modulo_lv;
      vy=(-uz*uy)/modulo_lv;
      vz=(1-uz*uz)/modulo_lv;

      wx=uy*vz-uz*vy;
      wy=uz*vx-ux*vz;
      wz=ux*vy-uy*vx;

        }

     /* Verificar se sao unitarios */


      uu = sqrt(pow(ux,2)+pow(uy,2)+pow(uz,2));

      vv = sqrt(pow(vx,2)+pow(vy,2)+pow(vz,2));

      ww = sqrt(pow(wx,2)+pow(wy,2)+pow(wz,2));


      /* Verificar se sao ortogonais */


      prod_uv=ux*vx+uy*vy+uz*vz;

      prod_uw=ux*wx+uy*wy+uz*wz;

      prod_vw=vx*wx+vy*wy+vz*wz;


      printf("%verificacao se eh unitario\n ");

      printf("%i\n ",uu);
      printf("%i\n ",vv);
      printf("%i\n ",ww);

      printf("%Verificacao se eh ortogonal\n ",uu);

      printf("%i\n ",prod_uv);
      printf("%i\n ",prod_uw);
      printf("%i\n ",prod_vw);

      printf("Valor dos coeficientes da matriz de rotacao do elemento %i\n",r+1);

      printf("valor de ux: %d \n",ux);
      fprintf(fileout,"valor de ux %d \n",ux);

      printf("valor de uy: %d \n",uy);
      fprintf(fileout,"valor de uy %d \n",uy);

      printf("valor de uz: %d \n",uz);
      fprintf(fileout,"valor de uz %d \n",uz);

      printf("valor de vx: %d \n",vx);
      fprintf(fileout,"valor de vx %d \n",vx);

      printf("valor de vy: %d \n",vy);
      fprintf(fileout,"valor de vy %d \n",vy);

      printf("valor de vz: %d \n",vz);
      fprintf(fileout,"valor de vz %d \n",vz);

      printf("valor de wx: %d \n",wx);
      fprintf(fileout,"valor de wx %d \n",wx);

      printf("valor de wy: %d \n",wy);
      fprintf(fileout,"valor de wy %d \n",wy);

      printf("valor de wz: %d \n",wz);
      fprintf(fileout,"valor de wz %d \n",wz);



       printf("\n ");
       fprintf(fileout,"\n");

       /* Propriedades do material */

       eu= mat[elem[r][2]][0];
       gu = mat[elem[r][2]][1];

       /* Propriedades geometricas da secao */

       au= sec[elem[r][3]][0];
       iuy= sec[elem[r][3]][1];
       iuz= sec[elem[r][3]][2];
       ju = sec[elem[r][3]][3];


       /* Determinacao do comprimento do elemento de barra */
       lu=sqrt(pow(x2-x1,2) + pow(y2-y1,2) + pow(z2-z1,2));

            /* Calculo dos coeficiente da matriz de rigidez */
            rig[0][0]=eu*au/lu;
            rig[0][6]=-eu*au/lu;

            rig[1][1]=(12*eu*iuz)/pow(lu,3);
            rig[1][5]=(6*eu*iuz)/pow(lu,2);
            rig[1][7]=(-12*eu*iuz)/pow(lu,3);
            rig[1][11]=(6*eu*iuz)/pow(lu,2);

            rig[2][2]=(12*eu*iuy)/pow(lu,3);
            rig[2][4]=(-6*eu*iuy)/pow(lu,2);
            rig[2][8]=(-12*eu*iuy)/pow(lu,3);
            rig[2][10]=(-6*eu*iuy)/pow(lu,2);

            rig[3][3]=(gu*ju)/lu;
            rig[3][9]=(-gu*ju)/lu;

            rig[4][2]=(-6*eu*iuy)/pow(lu,2);
            rig[4][4]=(4*eu*iuy)/lu;
            rig[4][8]=(6*eu*iuy)/pow(lu,2);
            rig[4][10]=(2*eu*iuy)/lu;

            rig[5][1]=(6*eu*iuz)/pow(lu,2);
            rig[5][5]=(4*eu*iuz)/lu;
            rig[5][7]=(-6*eu*iuz)/pow(lu,2);
            rig[5][11]=(2*eu*iuz)/lu;

            rig[6][0]=(-eu*au)/lu;
            rig[6][6]=(eu*au)/lu;

            rig[7][1]=(-12*eu*iuz)/pow(lu,3);
            rig[7][5]=(-6*eu*iuz)/pow(lu,2);
            rig[7][7]=(12*eu*iuz)/pow(lu,3);
            rig[7][11]=(-6*eu*iuz)/pow(lu,2);

            rig[8][2]=(-12*eu*iuy)/pow(lu,3);
            rig[8][4]=(6*eu*iuy)/pow(lu,2);
            rig[8][8]=(12*eu*iuy)/pow(lu,3);
            rig[8][10]=(6*eu*iuy)/pow(lu,2);

            rig[9][3]=(-gu*ju)/lu;
            rig[9][9]=(gu*ju)/lu;

            rig[10][2]=(-6*eu*iuy)/pow(lu,2);
            rig[10][4]=(2*eu*iuy)/lu;
            rig[10][8]=(6*eu*iuy)/pow(lu,2);
            rig[10][10]=(4*eu*iuy)/lu;

            rig[11][1]=(6*eu*iuz)/pow(lu,2);
            rig[11][5]=(2*eu*iuz)/lu;
            rig[11][7]=(-6*eu*iuz)/pow(lu,2);
            rig[11][11]=(4*eu*iuz)/lu;

            /* Calculo dos coeficientes da matriz de rotacao */

             rot[0][0]=ux;
             rot[0][1]=vx;
             rot[0][2]=wx;
             rot[1][0]=uy;
             rot[1][1]=vy;
             rot[1][2]=wy;
             rot[2][0]=uz;
             rot[2][1]=vz;
             rot[2][2]=wz;

             rot[3][3]=ux;
             rot[3][4]=vx;
             rot[3][5]=wx;
             rot[4][3]=uy;
             rot[4][4]=vy;
             rot[4][5]=wy;
             rot[5][3]=uz;
             rot[5][4]=vz;
             rot[5][5]=wz;

             rot[6][6]=ux;
             rot[6][7]=vx;
             rot[6][8]=wx;
             rot[7][6]=uy;
             rot[7][7]=vy;
             rot[7][8]=wy;
             rot[8][6]=uz;
             rot[8][7]=vz;
             rot[8][8]=wz;

             rot[9][9]=ux;
             rot[9][10]=vx;
             rot[9][11]=wx;
             rot[10][9]=uy;
             rot[10][10]=vy;
             rot[10][11]=wy;
             rot[11][9]=uz;
             rot[11][10]=vz;
             rot[11][11]=wz;


     /* Montagem da matriz de rigidez nas coordenadas locais do elemento */

     fprintf(fileout,"Matriz de rigidez nas coordenadas locais do elemento %i\n",r+1);
     printf("Matriz de rigidez nas coordenadas locais do elemento %i\n",r+1);

     for(li=0;li<12;li++)
      {
        for(co=0;co<12;co++)
        {
         printf("%lf  ",rig[li][co]);
         fprintf(fileout,"rig[%d][%d] = %lf  ",li+1,co+1,rig[li][co]);
        }
         printf("\n");
         fprintf(fileout,"\n");
      }
         printf("\n");
         fprintf(fileout,"\n");


     /* Montagem da matriz de rotacao do elemento local para global */

       fprintf(fileout,"Matriz de rotacao do elemento: %i\n",r+1);
       printf("Matriz de rotacao do elemento: %i\n",r+1);
       for(li=0;li<12;li++)
        {
         for(co=0;co<12;co++)
          {
          printf("%lf  ",rot[li][co]);
          fprintf(fileout,"rot[%d][%d] = %lf  ",li+1,co+1,rot[li][co]);
          }

          printf("\n");
          fprintf(fileout,"\n");
        }
          printf("\n");
          fprintf(fileout,"\n");

       /* }  */

         /* Multiplicação da matriz de rotacao pela matriz de rigidez do elemento local */

           printf(" Multiplicacao da matriz de rotacao pela matriz de rigidez do elemento local: %i\n",r+1);
           fprintf(fileout,"Multiplicacao da matriz de rotacao pela matriz de rigidez do elemento local: %i\n",r+1);

         for(li=0;li<12;li++)
          {
           for(co=0;co<12;co++)
           {
            soma1=0;
            for(comu=0;comu<12;comu++)
             {
              soma1+=rot[li][comu]*rig[comu][co];
             }

            pm[li][co]=soma1;
            printf("pm[%d][%d] = %lf  ",li+1,co+1,pm[li][co]);
            fprintf(fileout,"pm[%d][%d] = %lf  ",li+1,co+1,pm[li][co]);

           }
            printf("\n");
            fprintf(fileout,"\n");
          }
            printf("\n");
            fprintf(fileout,"\n");

           /* Multiplicacao do resultado do produto anterior pela transposta da matriz de rotacao  */

           printf("Matriz de rigidez nas coordenadas globais do elemento: %i\n",r+1);
           fprintf(fileout,"Matriz de rigidez nas coordenadas globais do elemento: %i\n",r+1);

           for(li=0;li<12;li++)
            {
             for(co=0;co<12;co++)
             {
              soma2=0;
              for(comu=0;comu<12;comu++)
               {
                soma2+=pm[li][comu]*rot[co][comu];
               }
                meg[r][li][co]=soma2;

                printf("meg[%d][%d][%d] = %lf  ",r+1,li+1,co+1,soma2);
                fprintf(fileout,"meg[%d][%d][%d] = %lf  ",r+1,li+1,co+1,soma2);

             }
                printf("\n");
                fprintf(fileout,"\n");
             }
                printf("\n");
                fprintf(fileout,"\n");

            /*Montagem da matriz global da estrutura */

        for(li=0;li<12;li++)
          {
            lin_est= mapael[r][li];
            if(lin_est!=0)
             {
              for(co=0;co<12;co++)
              {
               col_est= mapael[r][co];

               if(col_est!=0)
                {
                 mre[lin_est-1][col_est-1]+=meg[r][li][co];
                }

              }
             }
          }

    }
                printf("Impressao da matriz de rigidez da estrutura nas coordenadas globais: \n");
                fprintf(fileout,"Impressao da matriz de rigidez da estrutura nas coordenadas globais: \n");

             for(li=0;li<neq;li++)
              {
               for(co=0;co<neq;co++)
                {
                printf("mre[%d][%d] = %lf  ",li+1,co+1,mre[li][co]);
                fprintf(fileout,"mre[%d][%d] = %lf  ",li+1,co+1,mre[li][co]);
                }
                printf("\n");
                fprintf(fileout,"\n");
              }
                printf("\n");
                fprintf(fileout,"\n");


                /*Mapeamento das cargas externa atuantes na estrutura*/


              for(i=0;i<qn;i++)
               {
                 for(j=0;j<6;j++)
                 {
                  lin_est=mapano[i][j];
                  if(lin_est!=0)
                   {
                    vet_for[lin_est-1]=carga[i][j];
                   }
                 }
                }

              for(i=0;i<neq;i++)
               {
                printf("vet_for[%d] = %lf  \n",i+1,vet_for[i]);
                fprintf(fileout,"vet_for[%d] = %lf  \n",i+1,vet_for[i]);
               }


/*SOLUCAO PELO METODO DE CHOLESKY*/

/*-----------------------------------------------------------------------------
* Funcao: nrerror
* Descricao: MOSTRAR UMA MENSAGEM E SAIR
*-----------------------------------------------------------------------------*/

void nrerror(char error_text[])
{
printf("%s\n", error_text);
exit(0);
};

/*-----------------------------------------------------------------------------
* Funcao: checkSymmetric
* Descricao: VERIFIQUE SE A MATRIZ QUADRADA É SIMETRICA
*-----------------------------------------------------------------------------*/


void checkSymmetric(double **A, int n) {
int i,j;
for (i = 0; i < neq; i++) {
for (j = i; j < neq; j++)
if(A[i][j] != A[j][i])
nrerror("SEM SIMETRIA!");
}
}

/*-----------------------------------------------------------------------------
* Funcaao: showSystem
* Descricao: ENVIE O SISTEMA PARA A SAÍDA PADRAO
*-----------------------------------------------------------------------------*/


void showSystem(double **A, int n, char c) {
int i,j;
printf("\n%*c[%c]\n", 7*neq,' ', c);
for (i = 0; i < neq; i++) {
for (j = 0; j < neq; j++)
printf("|%+02.5lf|\t ", A[i][j]);
printf("\n");
}
}

/*-----------------------------------------------------------------------------
* Funcao: showArray
* Descricao: ENVIE O ARRAY PARA A SAÍDA PADRAO
*-----------------------------------------------------------------------------*/


void showArray(double *X, int n, char c) {
int i;
printf("\n");
for (i = 0; i < n; i++) {
printf("%c[%d] = %2.9lf\n",c,i, X[i]);
}
}
/*-----------------------------------------------------------------------------
* Funcao: createArray
* Descricao: RETORNE UM PONTO PARA UMA MATRIZ [N]
*-----------------------------------------------------------------------------*/


double *createArray(int N)
{
double *array;

if (( array = malloc( N*sizeof( double ))) == NULL ) {
nrerror("SEM MEMORIA!");
}

return array;
}

/*-----------------------------------------------------------------------------
* Funcao: createMatrix
* Descricao: RETORNE UM PONTEIRO PARA UMA MATRIZ [N,N]
*-----------------------------------------------------------------------------*/


double **createMatrix(int N)
{
double **matrix;
int i;

/* N E O NUMERO DE LINHAS */


if (( matrix = malloc( N*sizeof( double* ))) == NULL ) {
nrerror("SEM MEMORIA!");
}

/* N E O NUMERO DE COLUNAS */


for ( i = 0; i < neq; i++ )
{
if (( matrix[i] = malloc( neq*sizeof( double ))) == NULL )
{
nrerror("SEM MEMORIA!");
}
}

/* NxN MATRIZ */


return matrix;
}

/*-----------------------------------------------------------------------------
* Funcao: Cholesky
* Descricao: DECOMPOSICAO DE CHOLESKY
*-----------------------------------------------------------------------------*/


double **cholesky(double **A, int n)
{
int i, j, k;
double sum;
double **G;

/* CRIAR UMA MATRIZ COM O MESMO TAMANHO DE A */

G = createMatrix(neq);

for (k = 0; k < neq; k++)
for (i = 0; i <= k; i++)
{
/* A SOMA E A MESMA PARA AMBAS AS SITUACOES */


sum = 0;
for (j = 0; j < i; j++)
{
sum += G[i][j] * G[k][j];
}

/* CALCULAR O Gkk E Gjk */


if (i == k)
{
if(A[i][i] - sum <= 0)
nrerror("NAO DEFINIDA POSITIVA");
G[i][i] = sqrt(A[i][i] - sum);
}
else
G[k][i] = 1.0 / G[i][i] * (A[k][i] - sum);
}

/* RETORNE A DECOMPOSICAO */

return G;
}

/*-----------------------------------------------------------------------------
* Funcao: choleskySolucion
* Descricao: RESOLVER O SISTEMA ANTERIOR CALCULADO
*-----------------------------------------------------------------------------*/

void choleskySolucion(double **a, int n, double b[], double x[])
{
int i,j;
double sum;

/* RESOLVER O SISTEMA E ARMAZENAR A RESPOSTA Y EM X */

for(i = 0; i < neq; ++i) {
sum = 0;

for(j = 0; j < i; ++j)
sum += a[i][j] * x[j];

x[i] = (b[i] - sum) / a[i][i];
}

/* showArray(x,n,'y'); */

/* RESOLVA O SISTEMA USANDO 'Y' CALCULADO ANTES, E ARMAZENE A RESPOSTA EM X */

for (i=neq-1;i>=0;i--) {
sum=0;
for(j = i+1; j < neq; j++)
sum += a[j][i] * x[j];

x[i] = (x[i] - sum) / a[i][i];
}

/* showArray(x,n,'x'); */

}

/*-----------------------------------------------------------------------------
* Funcao: main
* Descricao:
*-----------------------------------------------------------------------------*/

int resp;
double **A; /* MATRIX A */

double **G; /* CHOLESKY RESPOSTA */

double *B; /* ARRAY B */

double *X; /* RESULTADO */

int N;  /*ORDEM DA MATRIZ */

 N=neq;

/***************************************************************************

*/


A = createMatrix(neq);
B = createArray(neq);
X = createArray(neq);


for (i=0;i<neq;i++)
  {
    for(j=0;j<neq;j++)
    {
       A[i][j]=mre[i][j];
    }
  }


for (i=0;i<neq;i++)
  {
    B[i]=vet_for[i];
  }


/***************************************************************************/

/* VERIFIQUE SE O SISTEMA LIDO E REALMENTE SIMETRICO */

checkSymmetric(A, N);

/* MOSTRE O SISTEMA LIDO */

showSystem(A,N,'A');

/* CALCULAR A MATRIZ CHOLESKY E ARMAZENAR EM G */

G = cholesky(A, N);

/* MOSTRE A MATRIZ CALCULADA COM O METODO CHOLESKY */

showSystem(G,N,'G');

/* RESOLVER O SISTEMA AX = B USANDO A MATRIZ CHOLESKY CALCULADO ANTES */

choleskySolucion(G, N, B, X);

/* MOSTRE A RESPOSTA  */

showArray(X,N,'x');


            for(i=0;i<neq;i++)
            {
              u[i]=X[i];
            }

              printf("\n");
              printf("Deslocamentos nodais da estrutura [m] \n");
              fprintf(fileout,"Deslocamentos nodais da estrutura [mm]:\n");
              printf("\n");
              for(i=0;i<neq;i++)
             {
               printf("u[%d] = %e\n",i+1,u[i]);
               fprintf(fileout,"u[%d] = %e  \n",i+1,u[i]);
             }

             /* Determinação das reacoes e esforcos internos  */

      for(r=0;r<qe;r++)
        {
          printf("elemento: %d\n",r+1);
          fprintf(fileout,"elemento: %d\n",r+1);

            for(li=0;li<12;li++)
             {
               keq=mapael[r][li];
               defor[li]=u[keq-1];
             }

        /* calculo das reacoes  */


            for(li=0;li<12;li++)
              {
                soma1=0;
                for(co=0;co<12;co++)
                  {
                    soma1+=meg[r][li][co]*defor[co];
                  }
                    veint[li]=soma1;
              }

                    no1=elem[r][0];
                    no2=elem[r][1];

                    for(i=0;i<6;i++)
                      {
                        k1=apoio[no1][i];

                        if(k1==1)

                         {
                            printf("reacao[%d][%d] = %lf \n",no1+1,i+1,veint[i]);
                            fprintf(fileout,"reacao[%d][%d] = %lf  \n",no1+1,i+1,veint[i]);

                              printf(" k1 = %d \n",k1);

                             }

                             k2=apoio[no2][i];

                               if(k2==1)
                              {
                                      printf(" k2 = %d \n",k2);
                                 k=i+7;
                              printf("reacao[%d][%d] = %lf [N] \n",no2+1,k,veint[i+6]);
                              fprintf(fileout,"reacao[%d][%d] = %lf [N] \n",no2+1,k,veint[i+6]);
                              }
                          }

        /* calculo dos esforcos internos  */
        for(li=0;li<12;li++)
        {
         soma1=0;
         for(co=0;co<12;co++)
            {
             soma1+=rot[co][li]*veint[co];
            }
            defor[li]=soma1;
        }


         for(li=0;li<12;li++)
          {
            printf("esforcos internos[%d][%d] = %lf  [N] \n",r+1,li+1,defor[li]);
            fprintf(fileout,"esforcos internos[%d][%d] = %lf  [N] \n",r+1,li+1,defor[li]);
          }

        }

             fclose(fileout);
             exit(0);


}

