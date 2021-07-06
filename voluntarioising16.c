#include <stdio.h>
#include <math.h>
#include "gsl_rng.h"

#define n 16//numero de nodos (lado)
#define Temp 5
#define npasos 1000000
#define nterm npasos/100
#define nsem 175649
#define e 2.718281828

gsl_rng *tau;

int main() //Función principal
{
extern gsl_rng *tau;
int semilla=nsem; //Definiciones para los números random
tau=gsl_rng_alloc(gsl_rng_taus); //Inicializamos el puntero
gsl_rng_set(tau,semilla); //Inicializamos la semilla

int i,j,k,dE,s[n][n],t,cu,cd,modo,puntero,w,Ek;
double b,r,p,en,E,sE,see,ee,cn,arg,mn,suma,error,ini,T,sumaf[n],vecE[nterm],vecee[nterm],vecm[nterm],mediaE,desviacionE,mediaee,desviacionee,mediam,desviacionm;
double errorE,erroree,errorm,errorc;
FILE *f1,*f2,*f3,*f4,*f5;
f1=fopen("spines.txt","w");
f2=fopen("m16.txt","w");
f3=fopen("f16.txt","w");
f4=fopen("cn16.txt","w");
f5=fopen("en16.txt","w");

cu=0;
Ek=0;
cd=0;
puntero=0;
mn=0.0;
suma=0.0;
E=0.0;
en=0.0;
cn=0.0;
ee=0.0;
see=0.0;
sE=0.0;
mediaE=0.0;
desviacionE=0.0;
mediaee=0.0;
desviacionee=0.0;
mediam=0.0;
desviacionm=0.0;

    for(i=0;i<=nterm-1;i++) //Inicializamos los nodos en +1
    {
        vecE[i]=0.0;
        vecee[i]=0.0;
        vecm[i]=0.0;
    }


    for(i=0;i<=n-1;i++) //Inicializamos los nodos en +1
    {
        sumaf[i]=0.0;
        for(j=0;j<=n-1;j++)
        {
            s[i][j]=1;
        }
    }




for(T=1.5;T<=3.5;T=T+0.2)
{
semilla=gsl_rng_uniform_int(tau,100000);
tau=gsl_rng_alloc(gsl_rng_taus); //Inicializamos el puntero
gsl_rng_set(tau,semilla); //Inicializamos la semilla

for(t=1;t<=npasos;t++)
{
    for(k=1;k<=n*n;k++)
    {

            i=gsl_rng_uniform_int(tau,n); 
            j=gsl_rng_uniform_int(tau,n);  //número aleatorio entero [0,n-1]. Tomamos así un spin aleatorio [i,j]



            if(i==0)   //Condiciones de contorno mediante if
            { 
                if(j==0)        //para i=-1-->i=N y j=-1-->j=N
                dE=2*s[i][j]*(s[i+1][j]+s[n-1][j]+s[i][j+1]+s[i][n-1]);
                if(j==n-1) //para i=-1-->i=N y j=N+1-->j=0
                dE=2*s[i][j]*(s[i+1][j]+s[n-1][j]+s[i][0]+s[i][j-1]);
                else
                dE=2*s[i][j]*(s[i+1][j]+s[n-1][j]+s[i][j+1]+s[i][j-1]);
            }


            else if(i==n-1)   //Condiciones de contorno mediante if
            { 
                if(j==0)        //para i=N+1-->i=0 y j=-1-->j=N
                dE=2*s[i][j]*(s[0][j]+s[i-1][j]+s[i][j+1]+s[i][n-1]);
                if(j==n-1) //para i=-1-->i=N y j=N+1-->j=0
                dE=2*s[i][j]*(s[0][j]+s[i-1][j]+s[i][0]+s[i][j-1]);
                else
                dE=2*s[i][j]*(s[0][j]+s[i-1][j]+s[i][j+1]+s[i][j-1]);
            }

            else
            {    //Si no hay problema en i
                if(j==0)        //para i=N+1-->i=0 y j=-1-->j=N
                dE=2*s[i][j]*(s[i+1][j]+s[i-1][j]+s[i][j+1]+s[i][n-1]);
                if(j==n-1) //para i=-1-->i=N y j=N+1-->j=0
                dE=2*s[i][j]*(s[i+1][j]+s[i-1][j]+s[i][0]+s[i][j-1]);
                else
                dE=2*s[i][j]*(s[i+1][j]+s[i-1][j]+s[i][j+1]+s[i][j-1]);
            }


    //Ya hemos calculado la dE. Tomamos un Epsilon (r) random real en [0,1] y comparamos para obtener p


    arg=-1.0*dE/T;   //argumento de la exponencial 
    //printf("dE:%i\n",dE);
    //printf("arg:%lf\n",arg);
    p=pow(e,arg);           //Calculamos la exponencial
    //printf("p:%lf\n",p);

    if(p>=1.0)          //Tomamos el mínimo entre 1 y la exp
    p=1.00;

    r=gsl_rng_uniform(tau); //Generamos el númro random entre [0,1]
    //printf("P_f:%lf\n",p);

    if(r<p) //Si r es menor que p, cambiamos el spin
    s[i][j]=-1*s[i][j];

    } //Terminamos un paso montecarlo 

    puntero=puntero+1;
    if(puntero==100) //100 pasos montecarlo
    {   
       for(i=0;i<=n-1;i++) //Calculamos la sumatoria
        {
            for(j=0;j<=n-1;j++)
            {
                suma=suma+s[i][j];
            }
        }  
        suma=sqrt(suma*suma)/(n*n); 
        mn=mn+suma;
        vecm[Ek]=suma;
        suma=0.0;
        puntero=0.0;  
        

    
        for(i=0;i<=n-1;i++)
        {
        for(j=0;j<=n-1;j++)
        {
            if(i==0)   //Condiciones de contorno mediante if
            { 
                if(j==0)      //para i=-1-->i=N y j=-1-->j=N
                E=E+s[i][j]*(s[i+1][j]+s[n-1][j]+s[i][j+1]+s[i][n-1]);
                else if(j==n-1) //para i=-1-->i=N y j=N+1-->j=0
                E=E+s[i][j]*(s[i+1][j]+s[n-1][j]+s[i][0]+s[i][j-1]);
                else
                E=E+s[i][j]*(s[i+1][j]+s[n-1][j]+s[i][j+1]+s[i][j-1]);
            }


            else if(i==n-1)   //Condiciones de contorno mediante if
            { 
                if(j==0)        //para i=N+1-->i=0 y j=-1-->j=N
                E=E+s[i][j]*(s[0][j]+s[i-1][j]+s[i][j+1]+s[i][n-1]);
                else if(j==n-1) //para i=-1-->i=N y j=N+1-->j=0
                E=E+s[i][j]*(s[0][j]+s[i-1][j]+s[i][0]+s[i][j-1]);
                else
                E=E+s[i][j]*(s[0][j]+s[i-1][j]+s[i][j+1]+s[i][j-1]);
            }

            else
            {  //Si no hay problema en i
                if(j==0)        //para i=N+1-->i=0 y j=-1-->j=N
                E=E+s[i][j]*(s[i+1][j]+s[i-1][j]+s[i][j+1]+s[i][n-1]);
                else if(j==n-1) //para i=-1-->i=N y j=N+1-->j=0
                E=E+s[i][j]*(s[i+1][j]+s[i-1][j]+s[i][0]+s[i][j-1]);
                else
                E=E+s[i][j]*(s[i+1][j]+s[i-1][j]+s[i][j+1]+s[i][j-1]);
            }
    }
    } //Ya acabada las umatoria de E
    vecE[Ek]=E; //Vector de energías para hacer la desviación
    vecee[Ek]=E*E;
    sE=sE+E;
    see=see+E*E; //E²
    E=0.0;
    Ek=Ek+1; 
 

    for(w=0;w<=n-1;w++)
    {
        for(i=0;i<=n-1;i++)//Para la correlacion
        {
            for(j=0;j<=n-1;j++)
            {
                if(i+w>=n-1)   //Condiciones de contorno mediante if
                { 
                    sumaf[w]=sumaf[w]+s[i][j]*s[i+w-n+1][j];
                }
                else
                sumaf[w]=sumaf[w]+s[i][j]*s[i+w][j];
                //printf("%lf\n",sumaf[w]);
            }
        }
    //printf("%i\t %lf\n",w,sumaf[w]);
    }
}//Cierre del if


} //Terminamos el número de pasos deseados

//Calculamos las magnitudes
mn=mn*100.0/npasos;
en=sE*100.0/(-4.0*npasos*n);
ee=see*0.25*100.0/npasos; //<E²>
cn=ee-en*en*4.0*n*n;
cn=cn/(n*n*T);

for(w=0;w<=n-1;w++)
{
    sumaf[w]=sumaf[w]/(n*n);
    sumaf[w]=sumaf[w]*100.0/npasos;
    fprintf(f3,"%lf\t %i\t %lf\n",T,w,sumaf[w]);
    sumaf[w]=0.0;
}

//Calculamos errores
for(i=0;i<=nterm-1;i++)
{
    mediaE=mediaE+vecE[i]; //sumatoria de la media
}
mediaE=mediaE*100/npasos;
for(i=0;i<=nterm-1;i++)
{
    desviacionE=desviacionE+(vecE[i]-mediaE)*(vecE[i]-mediaE); //sumatoria de la media
}
desviacionE=desviacionE*100/npasos;

for(i=0;i<=nterm-1;i++)
{
    mediaee=mediaee+vecee[i]; //sumatoria de la media
}
mediaee=mediaee*100/npasos;
for(i=0;i<=nterm-1;i++)
{
    desviacionee=desviacionee+(vecee[i]-mediaee)*(vecee[i]-mediaee); //sumatoria de la media
}
desviacionee=desviacionee*100/npasos;

for(i=0;i<=nterm-1;i++)
{
    mediam=mediam+vecm[i]; //sumatoria de la media
}
mediam=mediam*100/npasos;
for(i=0;i<=nterm-1;i++)
{
    desviacionm=desviacionm+(vecm[i]-mediam)*(vecm[i]-mediam); //sumatoria de la media
}
desviacionm=desviacionm*100/npasos;

errorE=3.0*desviacionE/sqrt(nterm);
erroree=3.0*desviacionee/sqrt(nterm);
errorm=3.0*desviacionm/sqrt(nterm);
errorc=mediaE*mediaE*mediaE*mediaE*erroree*erroree+(mediaee-2.0*mediaE)*(mediaee-2.0*mediaE)*errorE*errorE;
errorc=1.0*errorc/(1.0*nterm*nterm);

fprintf(f2,"%lf\t %lf\t %lf\n",T,mn,errorm);
fprintf(f5,"%lf\t %lf\t %lf\n",T,en,errorE/(2.0*nterm));
fprintf(f4,"%lf\t %lf\t %lf\n",T,cn,sqrt(errorc)/T);


//REseteamos
mn=0.0;
en=0.0; //<E>
ee=0.0;
cn=0.0;
E=0.0;
see=0.0;
sE=0.0;
desviacionE=0.0;
mediaE=0.0;
desviacionee=0.0;
mediaee=0.0;
desviacionm=0.0;
mediam=0.0;
Ek=0;
for(i=0;i<=nterm-1;i++)
{
    vecE[i]=0.0; //Reiniciamos
    vecee[i]=0.0;
    vecm[i]=0.0;
}
}//Cierre para el ciclo de T








return 0;
fclose(f1);
fclose(f2);
fclose(f3);
fclose(f4);
fclose(f5);
}//Cerramos el main


