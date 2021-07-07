#include <stdio.h>
#include <math.h>
#include "complex.h"
#include "gsl_rng.h" //Libreria para generación de números aleatorios
#define PI 3.14159
#define N 200 //número de puntos espaciales
#define nci 20   // número de ciclos
#define npasos 150
#define nsim 1000

gsl_rng *tau;

int main()
{
extern gsl_rng *tau; //Puntero al estado del número aleatorio
int semilla=75261265; //Semilla del generador de números aleatorios
tau=gsl_rng_alloc(gsl_rng_taus); //Inicializamos el puntero
gsl_rng_set(tau,semilla); //Inicializamos la semilla

int i,j,k,x,intervalo,mt,mi,m,keep,error,nd,condicion;
double v[N],mod[N];
double ko,s,norma_o,norma_o_par,norma_o_impar,norma,norma_par,norma_impar,chi,pd,pi,cten,lamda;
fcomplex alpha[N],Ao[N],phi[N],beta[N],b[N],xi[N];

FILE *f1,*f2,*f3;
f1=fopen("pot.txt","w");
f2=fopen("phi_t.txt","w");
f3=fopen("res200_20_v4.txt","w");



//A partir de estas condiciones iniciales calculamos Ko y con ello V, alpha, phi y Ao
keep=0;
ko=2.0*PI*nci/N;
s=1/(4*ko*ko);
alpha[N-1]=Complex(0.0,0.0);
phi[0]=Complex(0.0,0.0); //Funcion de onda nula en ambos extremos
phi[N-1]=Complex(0.0,0.0);
mod[0]=0.0;
mod[N-1]=0.0;
norma_o=0.0;
intervalo=0;
pd=0.0;
pi=0.0;
m=0;
mt=0;
mi=0;
cten=0.0;
error=0;
condicion=0;


for(i=1;i<=N-2;i++) //Conocido Ko, calculamos las phi_t0
{
    phi[i].r=cos(ko*i)*exp(-8.0*(4.0*i-N)*(4.0*i-N)/(N*N));
    phi[i].i=sin(ko*i)*exp(-8.0*(4.0*i-N)*(4.0*i-N)/(N*N));
    mod[i]=phi[i].i*phi[i].i+phi[i].r*phi[i].r;
} 

//Calculamos la norma para las phi_o con el metodo del Simpson compuesto
x=0.0;
for(k=1;k<=N/2;k++) //sumatoria impar
{
    x=(2*k-1);
    norma_o_par=norma_o_par+mod[x]; 
}

x=0.0;//reiniciamos las variables acumulativas

for(k=1;k<=(N/2-1);k++) //sumatoria par
{
      x=2*k;   
      norma_o_impar=norma_o_impar+mod[x]; 
}

//Calculamos la norma total

norma_o=(1.0/3.0)*(4.0*norma_o_impar+2.0*norma_o_par);




//printf("%lf\n",norma_o);
nd=80;
lamda=1;



for(i=0;i<=N-1;i++)
{
    if(i>=2*N/5-1)
    {
        if(i<=3*N/5-1)
            v[i]=lamda*ko*ko;
    }
    else
    v[i]=0.0;   //Ahora conocido V sacamos Ao
    Ao[i]=Complex(-2.0-v[i],2.0/s); 
    //printf("Ao:%lf\t %lf\n",Ao[i].r,Ao[i].i);
} 

for(i=1;i<=N-1;i++)
{
    //Conocido Ao, calculamos alpha de arriba a abajo
    alpha[N-1-i]=Cdiv(Complex(-1.0,0.0),Cadd(alpha[N-i],Ao[N-i]));
    //printf("alpha:%lf\t %lf\n",alpha[i].r,alpha[i].i); 
} 

/*   AQUI ACABA LA PARTE DE CONDICIONES INICIALES Y ATEMPORALES. AHORA INTRODUCIMOS UN BUCLE DONDE VARÍA EL TIEMPO */
for(m=0;m<=nsim;m++)
{

for(i=1;i<=N-2;i++) //Conocido Ko, calculamos las phi_t0
{
    phi[i].r=cos(ko*i)*exp(-8.0*(4.0*i-N)*(4.0*i-N)/(N*N));
    phi[i].i=sin(ko*i)*exp(-8.0*(4.0*i-N)*(4.0*i-N)/(N*N));
    mod[i]=phi[i].i*phi[i].i+phi[i].r*phi[i].r;
} 

//Calculamos la norma para las phi_o con el metodo del Simpson compuesto
x=0.0;
for(k=1;k<=N/2;k++) //sumatoria impar
{
    x=(2*k-1);
    norma_o_par=norma_o_par+mod[x]; 
}

x=0.0;//reiniciamos las variables acumulativas

for(k=1;k<=(N/2-1);k++) //sumatoria par
{
      x=2*k;   
      norma_o_impar=norma_o_impar+mod[x]; 
}

//Calculamos la norma total

norma_o=(1.0/3.0)*(4.0*norma_o_impar+2.0*norma_o_par);

while(keep==0 && error!=20)
{
intervalo=intervalo+1;
    //Calculamos las Beta dependientes del tiempo

    //Imponemos beta[N-1]==0 para todo T
    beta[N-1]=Complex(0.0,0.0);

    for(i=0;i<=N-1;i++)//Calculamos las b primero para cada T del bucle
    {
        b[i]=Cmul(phi[i],Complex(0.0,(4.0/s)));  
        //printf("b:%lf\t %lf\n",b[i].r,b[i].i);      
    } 


    for(i=1;i<=N-1;i++)//Calculamos las Beta para cada T
    {
        beta[N-1-i]=Cdiv(Csub(b[N-i],beta[N-i]),Cadd(alpha[N-i],Ao[N-i]));
        //printf("Beta:%lf\t %lf\n",beta[i].r,beta[i].i);
    } 

    //Conocidas Beta y alpha, imponemos Xi en 0 y N nulas y calculamos e resto de Xi
    xi[0]=Complex(0.0,0.0);
    xi[N-1]=Complex(0.0,0.0);

    for(i=1;i<=N-2;i++)//Calculamos las Xi para cada T
    {
        xi[i]=Cadd(Cmul(alpha[i-1],xi[i-1]),beta[i-1]);  
        //printf("xi:%lf\t %lf\n",xi[i].r,xi[i].i);      
    } 


    //Recalculamos las funciones de onda:
    norma=0.0;
    for(i=0;i<=N-1;i++)//Calculamos las phi para el nuevo T
    {
        phi[i]=Csub(xi[i],phi[i]); 
        mod[i]=phi[i].i*phi[i].i+phi[i].r*phi[i].r;//mod al cuadrado
    }


//YA HA EVOLUCIONADO LA ONDA. CUANDO SE ALCANCE EL INTERVALO DE ND
    if(intervalo==nd)
    {
        chi=gsl_rng_uniform(tau);  //número aleatorio real [0,1]
        for(i=0;i<=N/5;i++)//Calculamos las probabilidades
        {
            pd=pd+mod[i+4*N/5-1];//mod es al cuadrado ya
        }
        if(chi<pd)
        {
            mt=mt+1;//Se ha transmitido
            pd=0.0; //reseteamos esa variable acumulativa
            keep=1; //Hacemos que pare de evolucionar            
        }
        
        else
        {
            pd=0.0; //reseteamos esa variable acumulativa
            for(i=0;i<=N/5;i++)//Proyectamos
            {
                phi[i+4*N/5-1]=Complex(0.0,0.0);
            }
            for(i=0;i<=N-1;i++)
            {
                cten=cten+phi[i].i*phi[i].i+phi[i].r*phi[i].r; //cte de normalizacion
            }
            for(i=0;i<=N-1;i++)//Recalculamos las phi
            {
                phi[i]=Cdiv(phi[i],Complex(sqrt(cten),0.0));
                mod[i]=phi[i].i*phi[i].i+phi[i].r*phi[i].r;//mod al cuadrado
            }
            //Resetamos la cten
            cten=0.0;

            //Repetir las prpbabilidad para izquierda
            chi=gsl_rng_uniform(tau);  //número aleatorio real [0,1]
            for(i=0;i<=N/5-1;i++)//Calculamos las probabilidades
            {
                pi=pi+mod[i];//mod es al cuadrado ya
            }
            if(chi<pi)
            {
                mi=mi+1;//Se ha transmitido
                pi=0.0; //reseteamos esa variable acumulativa
                keep=1; //Hacemos que pare de evolucionar                
            }
            else
            {
                pi=0.0; //reseteamos esa variable acumulativa
                for(i=0;i<=N/5-1;i++)//Proyectamos
                {
                    phi[i]=Complex(0.0,0.0);
                }
                for(i=0;i<=N-1;i++)
                {
                    cten=cten+phi[i].i*phi[i].i+phi[i].r*phi[i].r; //cte de normalizacion
                }
                for(i=0;i<=N/5-1;i++)//Recalculamos las phi
                {
                    phi[i]=Cdiv(phi[i],Complex(sqrt(cten),0.0));
                }
                //Resetamos la cten
                cten=0.0;
                error=error+1;
            }
            
        }
        intervalo=0; //reseteamos el nd
    }


}//Cierre del while    
error=0;
keep=0;
}//Cierre del numero de simulaciones
fprintf(f3,"%lf\t %lf\t %lf\t %lf\t %i\n",lamda,1.0*mt,1.0*m,1.0*mi,nd);




fclose(f1);
fclose(f2);
fclose(f3);
return 0;
}
