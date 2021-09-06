#include "projet.h"   //inclusion du fichier header
//inclusion des bibiliotheques standards
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double **matrice2n(int a){                          //creation de la matrice 2*n par allocation dynamique
    int i;

    double** p = malloc ( 2 * sizeof ( double* ) );
for(i=0;i<2;i++){
    p[i]= malloc ( a*sizeof ( double ) );}

    return p;
}

double **matrice4n(int a){                          //creation de la matrice 4*n par allocation dynamique
    int i;

    double** p = malloc ( 4 * sizeof ( double* ) );
for(i=0;i<4;i++){
    p[i]= malloc ( a*sizeof ( double ) );}

    return p;
}

//fct utilisee pour l''initialisation des donnes de la planete
planete init_planete(){
int i=0;
planete p;

m_planete terre,lune,soleil,mercure,venus,mars,jupiter,saturne,uranus,neptune;

//attribution a chaque planete sa masse correspondante
 terre.indice=3;
 terre.masse=m_terre;

 lune.indice=9;
 lune.masse=m_lune;

 soleil.indice=0;
 soleil.masse=m_soleil;

 mercure.indice=1;
 mercure.masse=m_mercure;

 venus.indice=2;
 venus.masse=m_venus;

 mars.indice=4;
 mars.masse=m_mars;

 jupiter.indice=5;
 jupiter.masse=m_jupiter;

 saturne.indice=6;
 saturne.masse=m_saturne;

 uranus.indice=7;
 uranus.masse=m_uranus;

 neptune.indice=8;
 neptune.masse=m_neptune;

m_planete liste[10]={terre,lune,soleil,mercure,venus,mars,jupiter,saturne,uranus,neptune};

//choix de la planete
printf("choisir la planete:\n0.soleil\n1.mercure\n2.venus\n3.terre\n4.mars\n5.jupiter\n6.saturne\n7.uranus\n8.neptune\n9.lune\n");
scanf(" %d",&p.indice);

//demande des condition initiales
 printf("les conditions initiales de la planete\n");
 scanf(" %lf %lf %lf %lf",&p.x0,&p.y0,&p.vx0,&p.vy0);

 //initialisation de la masse
 for(i=0;i<10;i++){
    if(liste[i].indice==p.indice){
       p.masse=liste[i].masse;}
       }

 return p;   //recuperation en sortie une structure "planete" initialisee
 }

//fct utilisee pour etape 1
double** position1(){
int i,choix;
double r;
float h;       //pas de discretisation
int n=100000;  //nombre de points de calcul
double ** resultat = matrice4n(n);   //creation de la matrice 4*n resultat pour recuperer les donnees en sortie

//matrices utilisées pour le calcul des positions et vitesses de la lune pour n point selon x et y
double ** x1 = matrice2n(n);
double ** y1 = matrice2n(n);
//matrices utilisées pour le calcul des positions et vitesses de la terre pour n point selon x et y
double ** x2 = matrice2n(n);
double ** y2 = matrice2n(n);


//initialisation des donnees
x1[0][0]=0;           //x0 de la lune
x1[1][0]=1022;        //vx0 de la lune

y1[0][0]=384400000;   //y0 de la lune
y1[1][0]=0;           //vy0 de la lune



x2[0][0]=0;           //x0 de la terre
x2[1][0]=0;           //vx0 de la terre

y2[0][0]=0;           //y0 de la terre
y2[1][0]=0;           //vy0 de la terre


r=sqrt(pow(x1[0][0]-x2[0][0],2)+pow(y1[0][0]-y2[0][0],2));    //distance entre la lune et la terre a l''instant 0

//matrices utilsee pour la methode de range kutta 4 pour le calcul des x et des vx de la lune
double**k11 = matrice2n(n);
double**k21 = matrice2n(n);
double**k31 = matrice2n(n);
double**k41 = matrice2n(n);
//matrices utilsee pour la methode de range kutta 4 pour le calcul des y et des vy de la lune
double**k12 = matrice2n(n);
double**k22 = matrice2n(n);
double**k32 = matrice2n(n);
double**k42 = matrice2n(n);

//matrices utilsee pour la methode de range kutta 4 pour le calcul des x et des vx de la terre
double**K11 = matrice2n(n);
double**K21 = matrice2n(n);
double**K31 = matrice2n(n);
double**K41 = matrice2n(n);
//matrices utilsee pour la methode de range kutta 4 pour le calcul des y et des vy de la terre
double**K12 = matrice2n(n);
double**K22 = matrice2n(n);
double**K32 = matrice2n(n);
double**K42 = matrice2n(n);


printf("\ndonnez le pas de calcul\n");   //demande a l''utilsateur de rentrer le pas de discretisation h
scanf(" %f",&h);

printf(" choisissez la methode de calcul:\n1.euleur\n2.range kutta 2\n3.range kutta 4\n4.diff finies\n");   //choix de la methode de calcul
scanf("%d",&choix);

switch(choix)
{
//euler
case 1:

    for (i=1;i<n;i++){
    x1[0][i]=x1[0][i-1]+h*x1[1][i-1];                        //calcul de x+1 de la lune
    x1[1][i]=x1[1][i-1]+h*(-G*m_terre/pow(r,3))*x1[0][i-1];  //calcul de vx+1 de la lune

    y1[0][i]=y1[0][i-1]+h*y1[1][i-1];                        //calcul de y+1 de la lune
    y1[1][i]=y1[1][i-1]+h*(-G*m_terre/pow(r,3))*y1[0][i-1];  //calcul de vy+1 de la lune

    x2[0][i]=x2[0][i-1]+h*x2[1][i-1];                        //calcul de x+1 de la terre
    x2[1][i]=x2[0][i-1]+h*(G*m_lune/pow(r,3))*x2[0][i-1];    //calcul de vx+1 de la terre

    y2[0][i]=y2[0][i-1]+h*y2[1][i-1];                        //calcul de y+1 de la terre
    y2[1][i]=y2[0][i-1]+h*(G*m_lune/pow(r,3))*y2[0][i-1];    //calcul de vy+1 de la terre

  r=sqrt(pow(x1[0][i]-x2[0][i],2)+pow(y1[0][i]-y2[0][i],2)); //calcul de la distance entre la lune et la terre a l'instant t+i

}
break;

//range kutta 2
case 2:

    for (i=1;i<n;i++){
    x1[0][i]=x1[0][i-1]+h*(2*x1[1][i-1]-h*G*m_terre*x1[0][i-1]/pow(r,3))/2;      //calcul de x+1 de la lune
    x1[1][i]=x1[1][i-1]+h*(-G*m_terre*(2*x1[0][i-1]+h*x1[1][i-1])/pow(r,3))/2;   //calcul de vx+1 de la lune

    y1[0][i]=y1[0][i-1]+h*(2*y1[1][i-1]-h*G*m_terre*y1[0][i-1]/pow(r,3))/2;      //calcul de y+1 de la lune
    y1[1][i]=y1[1][i-1]+h*(-G*m_terre*(2*y1[0][i-1]+h*y1[1][i-1])/pow(r,3))/2;   //calcul de vy+1 de la lune


    x2[0][i]=x2[0][i-1]+h*(2*x2[1][i-1]-h*G*m_lune*x2[0][i-1]/pow(r,3))/2;       //calcul de x+1 de la terre
    x2[1][i]=x2[1][i-1]+h*(G*m_lune*(2*x2[0][i-1]+h*x2[1][i-1])/pow(r,3))/2;     //calcul de vx+1 de la terre

    y2[0][i]=y2[0][i-1]+h*(2*y2[1][i-1]-h*G*m_lune*y2[0][i-1]/pow(r,3))/2;       //calcul de y+1 de la terre
    y2[1][i]=y2[1][i-1]+h*(G*m_lune*(2*y2[0][i-1]+h*y2[1][i-1])/pow(r,3))/2;     //calcul de vy+1 de la terre

    r=sqrt(pow(x1[0][i]-x2[0][i],2)+pow(y1[0][i]-y2[0][i],2));                   //calcul de la distance entre la lune et la terre a l'instant t+i
}


break;

//range kutta 4

case 3:

for (i=0;i<n;i++){
//x1 vx1
k11[0][i]=h*x1[1][i];
k11[1][i]=-h*G*m_terre*x1[0][i]/pow(r,3);

k21[0][i]=k11[0][i]+h*k11[1][i]/2;
k21[1][i]=k11[1][i]-G*m_terre*h*k11[0][i]/(2*pow(r,3));

k31[0][i]=k11[0][i]+h*k21[1][i]/2;
k31[1][i]=k11[1][i]-G*m_terre*h*k21[0][i]/(2*pow(r,3));

k41[0][i]=k11[0][i]+h*k31[1][i];
k41[1][i]=k11[1][i]-G*m_terre*h*k31[0][i]/pow(r,3);

x1[0][i+1]=x1[0][i]+(k11[0][i]+2*k21[0][i]+2*k31[0][i]+k41[0][i])/6;  //calcul de x+1 de la lune
x1[1][i+1]=x1[1][i]+(k11[1][i]+2*k21[1][i]+2*k31[1][i]+k41[1][i])/6;  //calcul de vx+1 de la lune

//y1 vy1

k12[0][i]=h*y1[1][i];
k12[1][i]=-h*G*m_terre*y1[0][i]/pow(r,3);

k22[0][i]=k12[0][i]+h*k12[1][i]/2;
k22[1][i]=k12[1][i]-G*m_terre*h*k12[0][i]/(2*pow(r,3));

k32[0][i]=k12[0][i]+h*k22[1][i]/2;
k32[1][i]=k12[1][i]-G*m_terre*h*k22[0][i]/(2*pow(r,3));

k42[0][i]=k12[0][i]+h*k32[1][i];
k42[1][i]=k12[1][i]-G*m_terre*h*k32[0][i]/pow(r,3);

y1[0][i+1]=y1[0][i]+(k12[0][i]+2*k22[0][i]+2*k32[0][i]+k42[0][i])/6;  //calcul de y+1 de la lune
y1[1][i+1]=y1[1][i]+(k12[1][i]+2*k22[1][i]+2*k32[1][i]+k42[1][i])/6;  //calcul de vy+1 de la lune

//x2 vx2
K11[0][i]=h*x2[1][i];
K11[1][i]=h*G*m_lune*x2[0][i]/pow(r,3);

K21[0][i]=K11[0][i]+h*K11[1][i]/2;
K21[1][i]=K11[1][i]+G*m_lune*h*K11[0][i]/(2*pow(r,3));

K31[0][i]=K11[0][i]+h*K21[1][i]/2;
K31[1][i]=K11[1][i]+G*m_lune*h*K21[0][i]/(2*pow(r,3));

K41[0][i]=K11[0][i]+h*K31[1][i];
K41[1][i]=K11[1][i]+G*m_lune*h*K31[0][i]/pow(r,3);

x2[0][i+1]=x2[0][i]+(K11[0][i]+2*K21[0][i]+2*K31[0][i]+K41[0][i])/6;  //calcul de x+1 de la terre
x2[1][i+1]=x2[1][i]+(K11[1][i]+2*K21[1][i]+2*K31[1][i]+K41[1][i])/6;  //calcul de vx+1 de la terre

//y2 vy2
K12[0][i]=h*y2[1][i];
K12[1][i]=h*G*m_lune*y2[0][i]/pow(r,3);

K22[0][i]=K12[0][i]+h*K12[1][i]/2;
K22[1][i]=K12[1][i]+G*m_lune*h*K12[0][i]/(2*pow(r,3));

K32[0][i]=K12[0][i]+h*K22[1][i]/2;
K32[1][i]=K12[1][i]+G*m_lune*h*K22[0][i]/(2*pow(r,3));

K42[0][i]=K12[0][i]+h*K32[1][i];
K42[1][i]=K12[1][i]+G*m_lune*h*K32[0][i]/pow(r,3);

y2[0][i+1]=y2[0][i]+(K12[0][i]+2*K22[0][i]+2*K32[0][i]+K42[0][i])/6;  //calcul de y+1 de la terre
y2[1][i+1]=y2[1][i]+(K12[1][i]+2*K22[1][i]+2*K32[1][i]+K42[1][i])/6;  //calcul de vy+1 de la terre


r=sqrt(pow(x1[0][i]-x2[0][i],2)+pow(y1[0][i]-y2[0][i],2));   //calcul de la distance entre la lune et la terre a l'instant t+i

}
break;

//differences finies
case 4:

    x1[0][-1]=x1[0][0]-h*x1[1][0];  //calcul de x a l''instant t-1 pour la lune
    y1[0][-1]=y1[0][0]-h*y1[1][0];  //calcul de y a l''instant t-1 pour la lune

    x2[0][-1]=x2[0][0]-h*x2[1][0];  //calcul de x a l''instant t-1 pour la terre
    y2[0][-1]=y2[0][0]-h*y2[1][0];  //calcul de y a l''instant t-1 pour la terre


    for (i=0;i<n;i++){


    x1[0][i+1]=2*x1[0][i]-x1[0][i-1]+h*h*(-G*m_terre/pow(r,3))*x1[0][i];  //calcul de x+1 de la lune

    y1[0][i+1]=2*y1[0][i]-y1[0][i-1]+h*h*(-G*m_terre/pow(r,3))*y1[0][i];  //calcul de y+1 de la lune

    x2[0][i+1]=2*x2[0][i]-x2[0][i-1]+h*h*(-G*m_lune/pow(r,3))*x2[0][i];   //calcul de x+1 de la terre

    y2[0][i+1]=2*y2[0][i]-y2[0][i-1]+h*h*(-G*m_lune/pow(r,3))*y2[0][i];   //calcul de y+1 de la terre


  r=sqrt(pow(x1[0][i]-x2[0][i],2)+pow(y1[0][i]-y2[0][i],2));              //calcul de la distance entre la lune et la terre a l'instant t+i
}


break;

}

//recuperation en sortie de la position de la lune et de la terre dans la matrice resultat
for (i=0;i<n;i++){
    resultat[0][i]=x1[0][i];
    resultat[1][i]=y1[0][i];
    resultat[2][i]=x2[0][i];
    resultat[3][i]=y2[0][i];
}

return resultat;
}

//fct utilisee pour etape 2
double** position2(){
double r;
float h;  //pas de discretisation h
int i;
int n=1000000;  //nombre de points de calcul
double ** resultat = matrice4n(n);   //creation de la matrice 4*n resultat pour recuperer les donnees en sortie

//matrices utilisées pour le calcul des positions et vitesses de la lune pour n point selon x et y
double ** x1 = matrice2n(n);
double ** y1 = matrice2n(n);
//matrices utilisées pour le calcul des positions et vitesses de la terre pour n point selon x et y
double ** x2 = matrice2n(n);
double ** y2 = matrice2n(n);


//initialisation des donnees
x1[0][0]=0;           //x0 de la lune
x1[1][0]=1022;        //vx0 de la lune

y1[0][0]=384400000;   //y0 de la lune
y1[1][0]=0;           //vy0 de la lune



x2[0][0]=0;           //x0 de la terre
x2[1][0]=0;           //vx0 de la terre

y2[0][0]=0;           //y0 de la terre
y2[1][0]=0;           //vy0 de la terre


printf("\ndonnez le pas de calcul\n");   //demande a l''utilsateur de rentrer le pas de discretisation h
scanf(" %f",&h);

r=sqrt(pow(x1[0][0]-x2[0][0],2)+pow(y1[0][0]-y2[0][0],2));    //calcul de la distance entre la lune et la terre a l''instant 0

    x1[0][-1]=x1[0][0]-h*x1[1][0];    //calcul de x a l''instant t-1 pour la lune
    y1[0][-1]=y1[0][0]-h*y1[1][0];    //calcul de y a l''instant t-1 pour la lune


    x2[0][-1]=x2[0][0]-h*x2[1][0];    //calcul de x a l''instant t-1 pour la terre
    y2[0][-1]=y2[0][0]-h*y2[1][0];    //calcul de y a l''instant t-1 pour la terre

    for (i=0;i<n;i++){


    x1[0][i+1]=2*x1[0][i]-x1[0][i-1]+h*h*(-G*m_terre/pow(r,3))*(x1[0][i]-x2[0][i]);  //calcul de x+1 de la lune

    y1[0][i+1]=2*y1[0][i]-y1[0][i-1]+h*h*(-G*m_terre/pow(r,3))*(y1[0][i]-y2[0][i]);  //calcul de y+1 de la lune

    x2[0][i+1]=2*x2[0][i]-x2[0][i-1]+h*h*(-G*m_lune/pow(r,3))*(x2[0][i]-x1[0][i]);   //calcul de x+1 de la terre

    y2[0][i+1]=2*y2[0][i]-y2[0][i-1]+h*h*(-G*m_lune/pow(r,3))*(y2[0][i]-y1[0][i]);   //calcul de y+1 de la terre


  r=sqrt(pow(x1[0][i]-x2[0][i],2)+pow(y1[0][i]-y2[0][i],2));   //calcul de la distance entre la lune et la terre a l''instant t+i
}

//recuperation en sortie de la position de la lune et de la terre dans la matrice resultat
for (i=0;i<n;i++){
    resultat[0][i]=x1[0][i];
    resultat[1][i]=y1[0][i];
    resultat[2][i]=x2[0][i];
    resultat[3][i]=y2[0][i];}

    return resultat;
}

//fct utilisee pour etape 3
double **position3(planete p){
int i;
double r, masse;
float h;   //pas de discretisation h
int n=100000;  //nombre de points de calcul
double ** resultat = matrice4n(n);   //creation de la matrice 4*n resultat pour recuperer les donnees en sortie

//matrices utilisées pour le calcul des positions et vitesses de la planete pour n point selon x et y
double ** x1 = matrice2n(n);
double ** y1 = matrice2n(n);
//matrices utilisées pour le calcul des positions et vitesses du soleil pour n point selon x et y
double ** x2 = matrice2n(n);
double ** y2 = matrice2n(n);


//initialisation des donnees
x1[0][0]=p.x0;   //x0 de la planete
x1[1][0]=p.vx0;  //vx0 de la planete

y1[0][0]=p.y0;   //y0 de la planete
y1[1][0]=p.vy0;  //vy0 de la planete



x2[0][0]=0;     //x0 du soleil
x2[1][0]=0;     //vx0 du soleil

y2[0][0]=0;     //y0 du soleil
y2[1][0]=0;     //vy0 du soleil


masse=p.masse;

r=sqrt(pow(x1[0][0]-x2[0][0],2)+pow(y1[0][0]-y2[0][0],2));  //calcul de la distance entre le soleil et la planete a l''insatant 0

printf("\ndonnez le pas de calcul\n");   //demande a l''utilsateur de rentrer le pas de discretisation h
scanf(" %f",&h);

//differences finies

    x1[0][-1]=x1[0][0]-h*x1[1][0];   //calcul de x a l''instant t-1 pour la planete
    y1[0][-1]=y1[0][0]-h*y1[1][0];   //calcul de y a l''instant t-1 pour la planete

    x2[0][-1]=x2[0][0]-h*x2[1][0];   //calcul de x a l''instant t-1 pour le soleil
    y2[0][-1]=y2[0][0]-h*y2[1][0];   //calcul de y a l''instant t-1 pour le soleil

    for (i=0;i<n;i++){
    x1[0][i+1]=2*x1[0][i]-x1[0][i-1]+h*h*(-G*m_soleil/pow(r,3))*x1[0][i];   //calcul de x+1 de la planete

    y1[0][i+1]=2*y1[0][i]-y1[0][i-1]+h*h*(-G*m_soleil/pow(r,3))*y1[0][i];   //calcul de y+1 de la planete

    x2[0][i+1]=2*x2[0][i]-x2[0][i-1]+h*h*(-G*masse/pow(r,3))*x2[0][i];      //calcul de x+1 du soleil

    y2[0][i+1]=2*y2[0][i]-y2[0][i-1]+h*h*(-G*masse/pow(r,3))*y2[0][i];      //calcul de y+1 du soleil


  r=sqrt(pow(x1[0][i]-x2[0][i],2)+pow(y1[0][i]-y2[0][i],2));    //calcul de la distance entre le soleil et la planete a l''insatant t+i
}

//recuperation en sortie de la position de la planete et du soleil dans la matrice resultat
for (i=0;i<n;i++){
    resultat[0][i]=x1[0][i];
    resultat[1][i]=y1[0][i];
    resultat[2][i]=x2[0][i];
    resultat[3][i]=y2[0][i];
}
return resultat;
}

//fct utilisee pour etape 4
void **interaction(void){
int i,n;
n=100000;  //nombre de points de calcul

double h1,h2,h3,h4,h5,h6,h7,h8;  //pas de discrétisation des différentes planetes

FILE* p1; p1=fopen("4.mercure.txt", "w");     //fichier texte pour recuperer les coordonnees de mercure
FILE* p2; p2=fopen("4.venus.txt", "w");       //fichier texte pour recuperer les coordonnees de venus
FILE* p3; p3=fopen("4.terre.txt", "w");       //fichier texte pour recuperer les coordonnees de la terre
FILE* p4; p4=fopen("4.mars.txt", "w");        //fichier texte pour recuperer les coordonnees de mars
FILE* p5; p5=fopen("4.jupiter.txt", "w");     //fichier texte pour recuperer les coordonnees de jupiter
FILE* p6; p6=fopen("4.saturne.txt", "w");     //fichier texte pour recuperer les coordonnees de saturne
FILE* p7; p7=fopen("4.uranus.txt", "w");      //fichier texte pour recuperer les coordonnees de uranus
FILE* p8; p8=fopen("4.neptune.txt", "w");     //fichier texte pour recuperer les coordonnees de neptune

h1=76;         //pas de discretisation de mercure
h2=194;        //pas de discretisation de venus
h3=315;        //pas de discretisation de la terre
h4=590;        //pas de discretisation de mars
h5=3785;       //pas de discretisation de jupiter
h6=9600;       //pas de discretisation de saturne
h7=27594;      //pas de discretisation de uranus
h8=52034.4;    //pas de discretisation de neptune

//matrices utilisées pour le calcul des positions  de mercure pour n point selon x et y
double**x1=matrice2n(n);
double**y1=matrice2n(n);

//matrices utilisées pour le calcul des positions  de venus pour n point selon x et y
double**x2=matrice2n(n);
double**y2=matrice2n(n);

//matrices utilisées pour le calcul des positions  de la terre pour n point selon x et y
double**x3=matrice2n(n);
double**y3=matrice2n(n);

//matrices utilisées pour le calcul des positions  de mars pour n point selon x et y
double**x4=matrice2n(n);
double**y4=matrice2n(n);

//matrices utilisées pour le calcul des positions  de jupiter pour n point selon x et y
double**x5=matrice2n(n);
double**y5=matrice2n(n);

//matrices utilisées pour le calcul des positions  de saturne pour n point selon x et y
double**x6=matrice2n(n);
double**y6=matrice2n(n);

//matrices utilisées pour le calcul des positions  de uranus pour n point selon x et y
double**x7=matrice2n(n);
double**y7=matrice2n(n);

//matrices utilisées pour le calcul des positions  de neptune pour n point selon x et y
double**x8=matrice2n(n);
double**y8=matrice2n(n);


//variables utilisee pour le calcul des distances entre les planetes
double a1,a8,a2,a3,a4,a5,a6,a7;
double b8,b1,b2,b3,b4,b5,b6,b7;
double c8,c1,c2,c3,c4,c5,c6,c7;
double d8,d1,d2,d3,d4,d5,d6,d7;
double e8,e1,e2,e3,e4,e5,e6,e7;
double f8,f1,f2,f3,f4,f5,f6,f7;
double g8,g1,g2,g3,g4,g5,g6,g7;
double j1,j2,j3,j4,j5,j6,j7,j8;

    a1=sqrt(pow(x1[0][0],2)+pow(y1[0][0],2));                    //distane entre mercure et le soleil en 0
    a2=sqrt(pow(x1[0][0]-x2[0][0],2)+pow(y1[0][0]-y2[0][0],2));  //distane entre mercure et venus en 0
    a3=sqrt(pow(x1[0][0]-x3[0][0],2)+pow(y1[0][0]-y3[0][0],2));  //distane entre mercure et la terre en 0
    a4=sqrt(pow(x1[0][0]-x4[0][0],2)+pow(y1[0][0]-y4[0][0],2));  //distane entre mercure et mars en 0
    a5=sqrt(pow(x1[0][0]-x5[0][0],2)+pow(y1[0][0]-y5[0][0],2));  //distane entre mercure et jupiter en 0
    a6=sqrt(pow(x1[0][0]-x6[0][0],2)+pow(y1[0][0]-y6[0][0],2));  //distane entre mercure et saturne en 0
    a7=sqrt(pow(x1[0][0]-x7[0][0],2)+pow(y1[0][0]-y7[0][0],2));  //distane entre mercure et uranus en 0
    a8=sqrt(pow(x1[0][0]-x8[0][0],2)+pow(y1[0][0]-y8[0][0],2));  //distane entre mercure et neptune en 0


    b1=sqrt(pow(x2[0][0]-x1[0][0],2)+pow(y2[0][0]-y1[0][0],2));  //distance entre venus et mercure en 0
    b2=sqrt(pow(x2[0][0],2)+pow(y2[0][0],2));                    //distance entre venus et le soleil en 0
    b3=sqrt(pow(x2[0][0]-x3[0][0],2)+pow(y2[0][0]-y3[0][0],2));  //distance entre venus et le terre en 0
    b4=sqrt(pow(x2[0][0]-x4[0][0],2)+pow(y2[0][0]-y4[0][0],2));  //distance entre venus et mars en 0
    b5=sqrt(pow(x2[0][0]-x5[0][0],2)+pow(y2[0][0]-y5[0][0],2));  //distance entre venus et jupiter en 0
    b6=sqrt(pow(x2[0][0]-x6[0][0],2)+pow(y2[0][0]-y6[0][0],2));  //distance entre venus et saturne en 0
    b7=sqrt(pow(x2[0][0]-x7[0][0],2)+pow(y2[0][0]-y7[0][0],2));  //distance entre venus et uranus en 0
    b8=sqrt(pow(x2[0][0]-x8[0][0],2)+pow(y2[0][0]-y8[0][0],2));  //distance entre venus et neptune en 0


    c1=sqrt(pow(x3[0][0]-x1[0][0],2)+pow(y3[0][0]-y1[0][0],2));  //distance entre la terre et mercure en 0
    c2=sqrt(pow(x3[0][0]-x2[0][0],2)+pow(y3[0][0]-y2[0][0],2));  //distance entre la terre et venus en 0
    c3=sqrt(pow(x3[0][0],2)+pow(x3[0][0],2));                    //distance entre la terre et le soleil en 0
    c4=sqrt(pow(x3[0][0]-x4[0][0],2)+pow(y3[0][0]-y4[0][0],2));  //distance entre la terre et mars en 0
    c5=sqrt(pow(x3[0][0]-x5[0][0],2)+pow(y3[0][0]-y5[0][0],2));  //distance entre la terre et jupiter en 0
    c6=sqrt(pow(x3[0][0]-x6[0][0],2)+pow(y3[0][0]-y6[0][0],2));  //distance entre la terre et saturne en 0
    c7=sqrt(pow(x3[0][0]-x7[0][0],2)+pow(y3[0][0]-y7[0][0],2));  //distance entre la terre et uranus en 0
    c8=sqrt(pow(x3[0][0]-x8[0][0],2)+pow(y3[0][0]-y8[0][0],2));  //distance entre la terre et neptune en 0


    d1=sqrt(pow(x4[0][0]-x1[0][0],2)+pow(y4[0][0]-y1[0][0],2));  //distance entre mars et mercure en 0
    d2=sqrt(pow(x4[0][0]-x2[0][0],2)+pow(y4[0][0]-y2[0][0],2));  //distance entre mars et venus en 0
    d3=sqrt(pow(x4[0][0]-x3[0][0],2)+pow(y4[0][0]-y4[0][0],2));  //distance entre mars et la terre en 0
    d4=sqrt(pow(x4[0][0],2)+pow(y4[0][0],2));                    //distance entre mars et le soleil en 0
    d5=sqrt(pow(x4[0][0]-x5[0][0],2)+pow(y4[0][0]-y5[0][0],2));  //distance entre mars et jupiter en 0
    d6=sqrt(pow(x4[0][0]-x6[0][0],2)+pow(y4[0][0]-y6[0][0],2));  //distance entre mars et saturne en 0
    d7=sqrt(pow(x4[0][0]-x7[0][0],2)+pow(y4[0][0]-y7[0][0],2));  //distance entre mars et uranus en 0
    d8=sqrt(pow(x4[0][0]-x8[0][0],2)+pow(y4[0][0]-y8[0][0],2));  //distance entre mars et neptune en 0


    e1=sqrt(pow(x5[0][0]-x1[0][0],2)+pow(y5[0][0]-y1[0][0],2));  //distance entre jupiter et mercure en 0
    e2=sqrt(pow(x5[0][0]-x2[0][0],2)+pow(y5[0][0]-y2[0][0],2));  //distance entre jupiter et venus en 0
    e3=sqrt(pow(x5[0][0]-x3[0][0],2)+pow(y5[0][0]-y3[0][0],2));  //distance entre jupiter et la terre en 0
    e4=sqrt(pow(x5[0][0]-x4[0][0],2)+pow(y5[0][0]-y4[0][0],2));  //distance entre jupiter et mars en 0
    e5=sqrt(pow(x5[0][0],2)+pow(y5[0][0],2));                    //distance entre jupiter et le soleil en 0
    e6=sqrt(pow(x5[0][0]-x6[0][0],2)+pow(y5[0][0]-y6[0][0],2));  //distance entre jupiter et saturne en 0
    e7=sqrt(pow(x5[0][0]-x7[0][0],2)+pow(y5[0][0]-y7[0][0],2));  //distance entre jupiter et uranus en 0
    e8=sqrt(pow(x5[0][0]-x8[0][0],2)+pow(y5[0][0]-y8[0][0],2));  //distance entre jupiter et neptune en 0


    f1=sqrt(pow(x6[0][0]-x1[0][0],2)+pow(y6[0][0]-y1[0][0],2));  //distance entre saturne et mercure en 0
    f2=sqrt(pow(x6[0][0]-x2[0][0],2)+pow(y6[0][0]-y2[0][0],2));  //distance entre saturne et venus en 0
    f3=sqrt(pow(x6[0][0]-x3[0][0],2)+pow(y6[0][0]-y3[0][0],2));  //distance entre saturne et la terre en 0
    f4=sqrt(pow(x6[0][0]-x4[0][0],2)+pow(y6[0][0]-y4[0][0],2));  //distance entre saturne et mars en 0
    f5=sqrt(pow(x6[0][0]-x5[0][0],2)+pow(y6[0][0]-y5[0][0],2));  //distance entre saturne et jupiter en 0
    f6=sqrt(pow(x6[0][0],2)+pow(y6[0][0],2));                    //distance entre saturne et le soleil en 0
    f7=sqrt(pow(x6[0][0]-x7[0][0],2)+pow(y6[0][0]-y7[0][0],2));  //distance entre saturne et uranus en 0
    f8=sqrt(pow(x6[0][0]-x8[0][0],2)+pow(y6[0][0]-y8[0][0],2));  //distance entre saturne et neptune en 0


    g1=sqrt(pow(x7[0][0]-x1[0][0],2)+pow(y7[0][0]-y1[0][0],2));  //distance entre uranus et mercure en 0
    g2=sqrt(pow(x7[0][0]-x2[0][0],2)+pow(y7[0][0]-y2[0][0],2));  //distance entre uranus et venus en 0
    g3=sqrt(pow(x7[0][0]-x3[0][0],2)+pow(y7[0][0]-y3[0][0],2));  //distance entre uranus et la terre en 0
    g4=sqrt(pow(x7[0][0]-x4[0][0],2)+pow(y7[0][0]-y4[0][0],2));  //distance entre uranus et mars en 0
    g5=sqrt(pow(x7[0][0]-x5[0][0],2)+pow(y7[0][0]-y5[0][0],2));  //distance entre uranus et jupiter en 0
    g6=sqrt(pow(x7[0][0]-x6[0][0],2)+pow(y7[0][0]-y6[0][0],2));  //distance entre uranus et saturne en 0
    g7=sqrt(pow(x7[0][0],2)+pow(y7[0][0],2));                    //distance entre uranus et le soleil en 0
    g8=sqrt(pow(x7[0][0]-x8[0][0],2)+pow(y7[0][0]-y8[0][0],2));  //distance entre uranus et neptune en 0


    j1=sqrt(pow(x8[0][0]-x1[0][0],2)+pow(y8[0][0]-y1[0][0],2));  //distance entre neptune et mercure en 0
    j2=sqrt(pow(x8[0][0]-x2[0][0],2)+pow(y8[0][0]-y2[0][0],2));  //distance entre neptune et venus en 0
    j3=sqrt(pow(x8[0][0]-x3[0][0],2)+pow(y8[0][0]-y3[0][0],2));  //distance entre neptune et la terre en 0
    j4=sqrt(pow(x8[0][0]-x4[0][0],2)+pow(y8[0][0]-y4[0][0],2));  //distance entre neptune et mars en 0
    j5=sqrt(pow(x8[0][0]-x5[0][0],2)+pow(y8[0][0]-y5[0][0],2));  //distance entre neptune et jupiter en 0
    j6=sqrt(pow(x8[0][0]-x6[0][0],2)+pow(y8[0][0]-y6[0][0],2));  //distance entre neptune et saturne en 0
    j7=sqrt(pow(x8[0][0]-x7[0][0],2)+pow(y8[0][0]-y7[0][0],2));  //distance entre neptune et uranus en 0
    j8=sqrt(pow(x8[0][0],2)+pow(y8[0][0],2));                    //distance entre neptune et le soleil en 0

    x1[0][0]=0;                      //x0 de mercure
    x1[1][0]=47367.93;               //vx0 de mercure
    x1[0][-1]=x1[0][0]-h1*x1[1][0];  //calcul de x a l''instant t-1 pour mercure

    y1[0][0]=57910000000;            //y0 de mercure
    y1[1][0]=0;                      //vy0 de mercure
    y1[0][-1]=y1[0][0]-h1*y1[1][0];  //calcul de y a l''instant t-1 pour mercure

    x2[0][0]=0;                      //x0 de venus
    x2[1][0]=35025.71;               //vx0 de venus
    x2[0][-1]=x2[0][0]-h2*x2[1][0];  //calcul de x a l''instant t-1 pour venus

    y2[0][0]=108200000000;           //y0 de venus
    y2[1][0]=0;                      //vy0 de venus
    y2[0][-1]=y2[0][0]-h2*y2[1][0];  //calcul de y a l''instant t-1 pour venus

    x3[0][0]=0;                      //x0 de la terre
    x3[1][0]=29780;                  //vx0 de la terre
    x3[0][-1]=x3[0][0]-h3*x3[1][0];  //calcul de x a l''instant t-1 pour le terre

    y3[0][0]=149600000000;           //y0 de la terre
    y3[1][0]=0;                      //vy0 de la terre
    y3[0][-1]=y3[0][0]-h3*y3[1][0];  //calcul de y a l''instant t-1 pour le terre

    x4[0][0]=0;                      //x0 de mars
    x4[1][0]=24080.2;                //vx0 de mars
    x4[0][-1]=x4[0][0]-h4*x4[1][0];  //calcul de x a l''instant t-1 pour mars

    y4[0][0]=227900000000;           //y0 de mars
    y4[1][0]=0;                      //vy0 de mars
    y4[0][-1]=y4[0][0]-h4*y4[1][0];  //calcul de y a l''instant t-1 pour mars

    x5[0][0]=0;                      //x0 de jupiter
    x5[1][0]=13058.5;                //vx0 de jupiter
    x5[0][-1]=x5[0][0]-h5*x5[1][0];  //calcul de x a l''instant t-1 pour jupiter

    y5[0][0]=778500000000;           //y0 de jupiter
    y5[1][0]=0;                      //vy0 de jupiter
    y5[0][-1]=y5[0][0]-h5*y5[1][0];  //calcul de y a l''instant t-1 pour jupiter


    x6[0][0]=0;                      //x0 de saturne
    x6[1][0]=9640.7;                 //vx0 de saturne
    x6[0][-1]=x6[0][0]-h6*x6[1][0];  //calcul de x a l''instant t-1 pour saturne

    y6[0][0]=1434000000000;          //y0 de saturne
    y6[1][0]=0;                      //vy0 de saturne
    y6[0][-1]=y6[0][0]-h6*y6[1][0];  //calcul de y a l''instant t-1 pour saturne

    x7[0][0]=0;                      //x0 de uranus
    x7[1][0]=6796.7;                 //vx0 de uranus
    x7[0][-1]=x7[0][0]-h7*x7[1][0];  //calcul de x a l''instant t-1 pour uranus

    y7[0][0]=2871000000000;          //y0 de uranus
    y7[1][0]=0;                      //vy0 de uranus
    y7[0][-1]=y7[0][0]-h7*y7[1][0];  //calcul de y a l''instant t-1 pour uranus

    x8[0][0]=0;                      //x0 de neptune
    x8[1][0]=5432.48;                //vx0 de neptune
    x8[0][-1]=x8[0][0]-h8*x8[1][0];  //calcul de x a l''instant t-1 pour neptune

    y8[0][0]=4495000000000;          //y0 de neptune
    y8[1][0]=0;                      //vy0 de neptune
    y8[0][-1]=y8[0][0]-h8*y8[1][0];  //calcul de y a l''instant t-1 pour neptune




for (i=0;i<n;i++){
//calcul de x+1 de la planete mercure
    x1[0][i+1]= 2*x1[0][i]-x1[0][i-1]-h1*h1*((-G*m_soleil/pow(a1,3))*x1[0][i]+(-G*m_terre/pow(a3,3))*(x1[0][i]-x3[0][i])
+(-G*m_venus/pow(a2,3))*(x1[0][i]-x2[0][i])+(-G*m_mars/pow(a4,3))*(x1[0][i]-x4[0][i])+(-G*m_jupiter/pow(a5,3))*(x1[0][i]-x5[0][i])
+(-G*m_saturne/pow(a6,3))*(x1[0][i]-x6[0][i])+(-G*m_uranus/pow(a7,3))*(x1[0][i]-x7[0][i])+(-G*m_neptune/pow(a8,3))*(x1[0][i]-x8[0][i]));

//calcul de y+1 de la planete mercure
    y1[0][i+1]=2*y1[0][i]-y1[0][i-1]-h1*h1*((-G*m_soleil/pow(a1,3))*y1[0][i]+(-G*m_terre/pow(a3,3))*(y1[0][i]-y3[0][i])+
(-G*m_venus/pow(a2,3))*(y1[0][i]-y2[0][i])+(-G*m_mars/pow(a4,3))*(y1[0][i]-y4[0][i])+(-G*m_jupiter/pow(a5,3))*(y1[0][i]-y5[0][i])+
(-G*m_saturne/pow(a6,3))*(y1[0][i]-y6[0][i])+ (-G*m_uranus/pow(a7,3))*(y1[0][i]-y7[0][i])+(-G*m_neptune/pow(a8,3))*(y1[0][i]-y8[0][i]));


//calcul de x+1 de la planete venus
  x2[0][i+1]=2*x2[0][i]-x2[0][i-1]-h2*h2*((-G*m_soleil/pow(b2,3))*x2[0][i]+(-G*m_terre/pow(b3,3))*(x2[0][i]-x3[0][i])+
(-G*m_mercure/pow(b1,3))*(x2[0][i]-x1[0][i])+ (-G*m_mars/pow(b4,3))*(x2[0][i]-x4[0][i])+(-G*m_jupiter/pow(b5,3))*(x2[0][i]-x5[0][i])+
(-G*m_saturne/pow(b6,3))*(x2[0][i]-x6[0][i])+ (-G*m_uranus/pow(b7,3))*(x2[0][i]-x7[0][i])+(-G*m_neptune/pow(b8,3))*(x2[0][i]-x8[0][i]));

//calcul de y+1 de la planete venus
   y2[0][i+1]=2*y2[0][i]-y2[0][i-1]-h2*h2*((-G*m_soleil/pow(b2,3))*y2[0][i]+(-G*m_terre/pow(b3,3))*(y2[0][i]-y3[0][i])+
(-G*m_mercure/pow(b1,3))*(y2[0][i]-y1[0][i])+ (-G*m_mars/pow(b4,3))*(y2[0][i]-y4[0][i])+(-G*m_jupiter/pow(b5,3))*(y2[0][i]-y5[0][i])+
(-G*m_saturne/pow(b6,3))*(y2[0][i]-y6[0][i])+ (-G*m_uranus/pow(b7,3))*(y2[0][i]-y7[0][i])+(-G*m_neptune/pow(b8,3))*(y2[0][i]-y8[0][i]));


//calcul de x+1 de la planete terre
  x3[0][i+1]=2*x3[0][i]-x3[0][i-1]-h3*h3*((-G*m_soleil/pow(c3,3))*x3[0][i]+(-G*m_mercure/pow(c1,3))*(x3[0][i]-x1[0][i])+
(-G*m_venus/pow(c2,3))*(x3[0][i]-x2[0][i])+ (-G*m_mars/pow(c4,3))*(x3[0][i]-x4[0][i])+(-G*m_jupiter/pow(c5,3))*(x3[0][i]-x5[0][i])+
(-G*m_saturne/pow(c6,3))*(x3[0][i]-x6[0][i])+ (-G*m_uranus/pow(c7,3))*(x3[0][i]-x7[0][i])+(-G*m_neptune/pow(c8,3))*(x3[0][i]-x8[0][i]));

//calcul de y+1 de la planete terre
    y3[0][i+1]=2*y3[0][i]-y3[0][i-1]-h3*h3*((-G*m_soleil/pow(c3,3))*y3[0][i]+(-G*m_mercure/pow(c1,3))*(y3[0][i]-y1[0][i])+
(-G*m_venus/pow(c2,3))*(y3[0][i]-y2[0][i])+ (-G*m_mars/pow(c4,3))*(y3[0][i]-y4[0][i])+(-G*m_jupiter/pow(c5,3))*(y3[0][i]-y5[0][i])+
(-G*m_saturne/pow(c6,3))*(y3[0][i]-y6[0][i])+ (-G*m_uranus/pow(c7,3))*(y3[0][i]-y7[0][i])+(-G*m_neptune/pow(c8,3))*(y3[0][i]-y8[0][i]));


//calcul de x+1 de la planete mars
  x4[0][i+1]=2*x4[0][i]-x4[0][i-1]-h4*h4*((-G*m_soleil/pow(d4,3))*x4[0][i]+(-G*m_terre/pow(d3,3))*(x4[0][i]-x3[0][i])+
(-G*m_venus/pow(d2,3))*(x4[0][i]-x2[0][i])+ (-G*m_mercure/pow(d1,3))*(x4[0][i]-x1[0][i])+(-G*m_jupiter/pow(d5,3))*(x4[0][i]-x5[0][i])+
(-G*m_saturne/pow(d6,3))*(x4[0][i]-x6[0][i])+ (-G*m_uranus/pow(d7,3))*(x4[0][i]-x7[0][i])+(-G*m_neptune/pow(d8,3))*(x4[0][i]-x8[0][i]));

//calcul de y+1 de la planete mars
   y4[0][i+1]=2*y4[0][i]-y4[0][i-1]-h4*h4*((-G*m_soleil/pow(d4,3))*y4[0][i]+(-G*m_terre/pow(d3,3))*(y4[0][i]-y3[0][i])+
(-G*m_venus/pow(d2,3))*(y4[0][i]-y2[0][i])+ (-G*m_mercure/pow(d1,3))*(y4[0][i]-y1[0][i])+(-G*m_jupiter/pow(d5,3))*(y4[0][i]-y5[0][i])+
(-G*m_saturne/pow(d6,3))*(y4[0][i]-y6[0][i])+ (-G*m_uranus/pow(d7,3))*(y4[0][i]-y7[0][i])+(-G*m_neptune/pow(d8,3))*(y4[0][i]-y8[0][i]));


//calcul de x+1 de la planete jupiter
  x5[0][i+1]=2*x5[0][i]-x5[0][i-1]-h5*h5*((-G*m_soleil/pow(e5,3))*x5[0][i]+(-G*m_terre/pow(e3,3))*(x5[0][i]-x3[0][i])+
(-G*m_venus/pow(e2,3))*(x5[0][i]-x2[0][i])+ (-G*m_mars/pow(e4,3))*(x5[0][i]-x4[0][i])+(-G*m_mercure/pow(e1,3))*(x5[0][i]-x1[0][i])+
(-G*m_saturne/pow(e6,3))*(x5[0][i]-x6[0][i])+ (-G*m_uranus/pow(e7,3))*(x5[0][i]-x7[0][i])+(-G*m_neptune/pow(e8,3))*(x5[0][i]-x8[0][i]));

//calcul de y+1 de la planete jupiter
   y5[0][i+1]=2*y5[0][i]-y5[0][i-1]-h5*h5*((-G*m_soleil/pow(e5,3))*y5[0][i]+(-G*m_terre/pow(e3,3))*(y5[0][i]-y3[0][i])+
(-G*m_venus/pow(e2,3))*(y5[0][i]-y2[0][i])+ (-G*m_mars/pow(e4,3))*(y5[0][i]-y4[0][i])+(-G*m_mercure/pow(e1,3))*(y5[0][i]-y1[0][i])+
(-G*m_saturne/pow(e6,3))*(y5[0][i]-y6[0][i])+ (-G*m_uranus/pow(e7,3))*(y5[0][i]-y7[0][i])+(-G*m_neptune/pow(e8,3))*(y5[0][i]-y8[0][i]));

//calcul de x+1 de la planete saturne
  x6[0][i+1]=2*x6[0][i]-x6[0][i-1]-h6*h6*((-G*m_soleil/pow(f6,3))*x6[0][i]+(-G*m_terre/pow(f3,3))*(x6[0][i]-x3[0][i])+
(-G*m_venus/pow(f2,3))*(x6[0][i]-x2[0][i])+ (-G*m_mars/pow(f4,3))*(x6[0][i]-x4[0][i])+(-G*m_jupiter/pow(f5,3))*(x6[0][i]-x5[0][i])+
(-G*m_mercure/pow(f1,3))*(x6[0][i]-x1[0][i])+ (-G*m_uranus/pow(f7,3))*(x6[0][i]-x7[0][i])+(-G*m_neptune/pow(f8,3))*(x6[0][i]-x8[0][i]));

//calcul de y+1 de la planete saturne
  y6[0][i+1]=2*y6[0][i]-y6[0][i-1]-h6*h6*((-G*m_soleil/pow(f6,3))*y6[0][i]+(-G*m_terre/pow(f3,3))*(y6[0][i]-y3[0][i])+
(-G*m_venus/pow(f2,3))*(y6[0][i]-y2[0][i])+ (-G*m_mars/pow(f4,3))*(y6[0][i]-y4[0][i])+(-G*m_jupiter/pow(f5,3))*(y6[0][i]-y5[0][i])+
(-G*m_mercure/pow(f1,3))*(y6[0][i]-y1[0][i])+ (-G*m_uranus/pow(f7,3))*(y6[0][i]-y7[0][i])+(-G*m_neptune/pow(f8,3))*(y6[0][i]-y8[0][i]));

//calcul de x+1 de la planete uranus
  x7[0][i+1]=2*x7[0][i]-x7[0][i-1]-h7*h7*((-G*m_soleil/pow(g7,3))*x7[0][i]+(-G*m_terre/pow(g3,3))*(x7[0][i]-x3[0][i])+
(-G*m_venus/pow(g2,3))*(x7[0][i]-x2[0][i])+ (-G*m_mars/pow(g4,3))*(x7[0][i]-x4[0][i])+(-G*m_jupiter/pow(g5,3))*(x7[0][i]-x5[0][i])+
(-G*m_saturne/pow(g6,3))*(x7[0][i]-x6[0][i])+ (-G*m_mercure/pow(g1,3))*(x7[0][i]-x1[0][i])+(-G*m_neptune/pow(g8,3))*(x7[0][i]-x8[0][i]));

//calcul de y+1 de la planete uranus
    y7[0][i+1]=2*y7[0][i]-y7[0][i-1]-h7*h7*((-G*m_soleil/pow(g7,3))*y7[0][i]+(-G*m_terre/pow(g3,3))*(y7[0][i]-y3[0][i])+
(-G*m_venus/pow(g2,3))*(y7[0][i]-y2[0][i])+ (-G*m_mars/pow(g4,3))*(y7[0][i]-y4[0][i])+(-G*m_jupiter/pow(g5,3))*(y7[0][i]-y5[0][i])+
(-G*m_saturne/pow(g6,3))*(y7[0][i]-y6[0][i])+ (-G*m_mercure/pow(g1,3))*(y7[0][i]-y1[0][i])+(-G*m_neptune/pow(g8,3))*(y7[0][i]-y8[0][i]));


//calcul de x+1 de la planete neptune
  x8[0][i+1]=2*x8[0][i]-x8[0][i-1]-h8*h8*((-G*m_soleil/pow(j8,3))*x8[0][i]+(-G*m_terre/pow(j3,3))*(x8[0][i]-x3[0][i])+
(-G*m_venus/pow(j2,3))*(x8[0][i]-x2[0][i])+ (-G*m_mars/pow(j4,3))*(x8[0][i]-x4[0][i])+(-G*m_jupiter/pow(j5,3))*(x8[0][i]-x5[0][i])+
(-G*m_saturne/pow(j6,3))*(x8[0][i]-x6[0][i])+ (-G*m_mercure/pow(j1,3))*(x8[0][i]-x1[0][i])+(-G*m_uranus/pow(j7,3))*(x8[0][i]-x7[0][i]));

//calcul de y+1 de la planete neptune
  y8[0][i+1]=2*y8[0][i]-y8[0][i-1]-h8*h8*((-G*m_soleil/pow(j8,3))*y8[0][i]+(-G*m_terre/pow(j3,3))*(y8[0][i]-y3[0][i])+
(-G*m_venus/pow(j2,3))*(y8[0][i]-y2[0][i])+ (-G*m_mars/pow(j4,3))*(y8[0][i]-y4[0][i])+(-G*m_jupiter/pow(j5,3))*(y8[0][i]-y5[0][i])+
(-G*m_saturne/pow(j6,3))*(y8[0][i]-y6[0][i])+ (-G*m_mercure/pow(j1,3))*(y8[0][i]-y1[0][i])+(-G*m_uranus/pow(j7,3))*(y8[0][i]-y7[0][i]));



    a1=sqrt(pow(x1[0][i],2)+pow(y1[0][i],2));                    //distance entre mercure et le soleilen t+i
    a2=sqrt(pow(x1[0][i]-x2[0][i],2)+pow(y1[0][i]-y2[0][i],2));  //distance entre mercure et venus en t+i
    a3=sqrt(pow(x1[0][i]-x3[0][i],2)+pow(y1[0][i]-y3[0][i],2));  //distance entre mercure et la terre en t+i
    a4=sqrt(pow(x1[0][i]-x4[0][i],2)+pow(y1[0][i]-y4[0][i],2));  //distance entre mercure et mars  en t+i
    a5=sqrt(pow(x1[0][i]-x5[0][i],2)+pow(y1[0][i]-y5[0][i],2));  //distance entre mercure et jupiter en t+i
    a6=sqrt(pow(x1[0][i]-x6[0][i],2)+pow(y1[0][i]-y6[0][i],2));  //distance entre mercure et saturne en t+i
    a7=sqrt(pow(x1[0][i]-x7[0][i],2)+pow(y1[0][i]-y7[0][i],2));  //distance entre mercure et uranus en t+i
    a8=sqrt(pow(x1[0][i]-x8[0][i],2)+pow(y1[0][i]-y8[0][i],2));  //distance entre mercure et neptune en t+i


    b1=sqrt(pow(x2[0][i]-x1[0][i],2)+pow(y2[0][i]-y1[0][i],2));  //distance entre venus et mercure en t+i
    b2=sqrt(pow(x2[0][i],2)+pow(y2[0][i],2));                    //distance entre venus et lz soleil en t+i
    b3=sqrt(pow(x2[0][i]-x3[0][i],2)+pow(y2[0][i]-y3[0][i],2));  //distance entre venus et la terre en t+i
    b4=sqrt(pow(x2[0][i]-x4[0][i],2)+pow(y2[0][i]-y4[0][i],2));  //distance entre venus et mars en t+i
    b5=sqrt(pow(x2[0][i]-x5[0][i],2)+pow(y2[0][i]-y5[0][i],2));  //distance entre venus et jupiter en t+i
    b6=sqrt(pow(x2[0][i]-x6[0][i],2)+pow(y2[0][i]-y6[0][i],2));  //distance entre venus et saturne en t+i
    b7=sqrt(pow(x2[0][i]-x7[0][i],2)+pow(y2[0][i]-y7[0][i],2));  //distance entre venus et uranus en t+i
    b8=sqrt(pow(x2[0][i]-x8[0][i],2)+pow(y2[0][i]-y8[0][i],2));  //distance entre venus et neptune en t+i


    c1=sqrt(pow(x3[0][i]-x1[0][i],2)+pow(y3[0][i]-y1[0][i],2));  //distance entre le terre et mercure en t+i
    c2=sqrt(pow(x3[0][i]-x2[0][i],2)+pow(y3[0][i]-y2[0][i],2));  //distance entre le terre et venus en t+i
    c3=sqrt(pow(x3[0][i],2)+pow(y3[0][i],2));                    //distance entre le terre et le soleil en t+i
    c4=sqrt(pow(x3[0][i]-x4[0][i],2)+pow(y3[0][i]-y4[0][i],2));  //distance entre le terre et mars en t+i
    c5=sqrt(pow(x3[0][i]-x5[0][i],2)+pow(y3[0][i]-y5[0][i],2));  //distance entre le terre et jupiter en t+i
    c6=sqrt(pow(x3[0][i]-x6[0][i],2)+pow(y3[0][i]-y6[0][i],2));  //distance entre le terre et saturne en t+i
    c7=sqrt(pow(x3[0][i]-x7[0][i],2)+pow(y3[0][i]-y7[0][i],2));  //distance entre le terre et uranus en t+i
    c8=sqrt(pow(x3[0][i]-x8[0][i],2)+pow(y3[0][i]-y8[0][i],2));  //distance entre le terre et neptune en t+i


    d1=sqrt(pow(x4[0][i]-x1[0][i],2)+pow(y4[0][i]-y1[0][i],2));  //distance entre mars et mercure en t+i
    d2=sqrt(pow(x4[0][i]-x2[0][i],2)+pow(y4[0][i]-y2[0][i],2));  //distance entre mars et venus en t+i
    d3=sqrt(pow(x4[0][i]-x3[0][i],2)+pow(y4[0][i]-y4[0][i],2));  //distance entre mars et la terre en t+i
    d4=sqrt(pow(x4[0][i],2)+pow(y4[0][i],2));                    //distance entre mars et le soleil en t+i
    d5=sqrt(pow(x4[0][i]-x5[0][i],2)+pow(y4[0][i]-y5[0][i],2));  //distance entre mars et jupiter en t+i
    d6=sqrt(pow(x4[0][i]-x6[0][i],2)+pow(y4[0][i]-y6[0][i],2));  //distance entre mars et saturne en t+i
    d7=sqrt(pow(x4[0][i]-x7[0][i],2)+pow(y4[0][i]-y7[0][i],2));  //distance entre mars et uranus en t+i
    d8=sqrt(pow(x4[0][i]-x8[0][i],2)+pow(y4[0][i]-y8[0][i],2));  //distance entre mars et neptune en t+i


    e1=sqrt(pow(x5[0][i]-x1[0][i],2)+pow(y5[0][i]-y1[0][i],2));  //distance entre jupiter et mercure en t+i
    e2=sqrt(pow(x5[0][i]-x2[0][i],2)+pow(y5[0][i]-y2[0][i],2));  //distance entre jupiter et venus en t+i
    e3=sqrt(pow(x5[0][i]-x3[0][i],2)+pow(y5[0][i]-y3[0][i],2));  //distance entre jupiter et la terre en t+i
    e4=sqrt(pow(x5[0][i]-x4[0][i],2)+pow(y5[0][i]-y4[0][i],2));  //distance entre jupiter et mars en t+i
    e5=sqrt(pow(x5[0][i],2)+pow(y5[0][i],2));                    //distance entre jupiter et le soleil en t+i
    e6=sqrt(pow(x5[0][i]-x6[0][i],2)+pow(y5[0][i]-y6[0][i],2));  //distance entre jupiter et saturne en t+i
    e7=sqrt(pow(x5[0][i]-x7[0][i],2)+pow(y5[0][i]-y7[0][i],2));  //distance entre jupiter et uranus en t+i
    e8=sqrt(pow(x5[0][i]-x8[0][i],2)+pow(y5[0][i]-y8[0][i],2));  //distance entre jupiter et neptune en t+i


    f1=sqrt(pow(x6[0][i]-x1[0][i],2)+pow(y6[0][i]-y1[0][i],2));  //distance entre saturne et mercure en t+i
    f2=sqrt(pow(x6[0][i]-x2[0][i],2)+pow(y6[0][i]-y2[0][i],2));  //distance entre saturne et venus en t+i
    f3=sqrt(pow(x6[0][i]-x3[0][i],2)+pow(y6[0][i]-y3[0][i],2));  //distance entre saturne et la terre en t+i
    f4=sqrt(pow(x6[0][i]-x4[0][i],2)+pow(y6[0][i]-y4[0][i],2));  //distance entre saturne et mars en t+i
    f5=sqrt(pow(x6[0][i]-x5[0][i],2)+pow(y6[0][i]-y5[0][i],2));  //distance entre saturne et jupiter en t+i
    f6=sqrt(pow(x6[0][i],2)+pow(y6[0][i],2));                    //distance entre saturne et le soleil en t+i
    f7=sqrt(pow(x6[0][i]-x7[0][i],2)+pow(y6[0][i]-y7[0][i],2));  //distance entre saturne et uranus en t+i
    f8=sqrt(pow(x6[0][i]-x8[0][i],2)+pow(y6[0][i]-y8[0][i],2));  //distance entre saturne et neptune en t+i


    g1=sqrt(pow(x7[0][i]-x1[0][i],2)+pow(y7[0][i]-y1[0][i],2));  //distance entre uranus et mercure en t+i
    g2=sqrt(pow(x7[0][i]-x2[0][i],2)+pow(y7[0][i]-y2[0][i],2));  //distance entre uranus et venus en t+i
    g3=sqrt(pow(x7[0][i]-x3[0][i],2)+pow(y7[0][i]-y3[0][i],2));  //distance entre uranus et la terre en t+i
    g4=sqrt(pow(x7[0][i]-x4[0][i],2)+pow(y7[0][i]-y4[0][i],2));  //distance entre uranus et mars en t+i
    g5=sqrt(pow(x7[0][i]-x5[0][i],2)+pow(y7[0][i]-y5[0][i],2));  //distance entre uranus et jupiter en t+i
    g6=sqrt(pow(x7[0][i]-x6[0][i],2)+pow(y7[0][i]-y6[0][i],2));  //distance entre uranus et saturne en t+i
    g7=sqrt(pow(x7[0][i],2)+pow(y7[0][i],2));                    //distance entre uranus et le soleil en t+i
    g8=sqrt(pow(x7[0][i]-x8[0][i],2)+pow(y7[0][i]-y8[0][i],2));  //distance entre uranus et neptune en t+i


    j1=sqrt(pow(x8[0][i]-x1[0][i],2)+pow(y8[0][i]-y1[0][i],2));  //distance entre neptune et mercure en t+i
    j2=sqrt(pow(x8[0][i]-x2[0][i],2)+pow(y8[0][i]-y2[0][i],2));  //distance entre neptune et venus en t+i
    j3=sqrt(pow(x8[0][i]-x3[0][i],2)+pow(y8[0][i]-y3[0][i],2));  //distance entre neptune et la terre en t+i
    j4=sqrt(pow(x8[0][i]-x4[0][i],2)+pow(y8[0][i]-y4[0][i],2));  //distance entre neptune et mars en t+i
    j5=sqrt(pow(x8[0][i]-x5[0][i],2)+pow(y8[0][i]-y5[0][i],2));  //distance entre neptune et jupiter en t+i
    j6=sqrt(pow(x8[0][i]-x6[0][i],2)+pow(y8[0][i]-y6[0][i],2));  //distance entre neptune et saturne en t+i
    j7=sqrt(pow(x8[0][i]-x7[0][i],2)+pow(y8[0][i]-y7[0][i],2));  //distance entre neptune et uranus en t+i
    j8=sqrt(pow(x8[0][i],2)+pow(y8[0][i],2));                    //distance entre neptune et le soleil en t+i
    }

//ecriture des coordonnees de mercure dans le fichier p1
   for(i=0;i<n;i++){
   if(p1!=NULL)
   fprintf(p1," %lf %lf\n",x1[0][i],y1[0][i]);
    else printf("erreur");
}

//ecriture des coordonnees de venus dans le fichier p2
 for(i=0;i<n;i++){
   if(p2!=NULL)
   fprintf(p2," %lf %lf\n",x2[0][i],y2[0][i]);
    else printf("erreur");
}

//ecriture des coordonnees de la terre dans le fichier p3
 for(i=0;i<n;i++){
   if(p3!=NULL)
   fprintf(p3," %lf %lf\n",x3[0][i],y3[0][i]);
    else printf("erreur");
}

//ecriture des coordonnees de mars dans le fichier p4
 for(i=0;i<n;i++){
   if(p4!=NULL)
   fprintf(p4," %lf %lf\n",x4[0][i],y4[0][i]);
    else printf("erreur");
}

//ecriture des coordonnees de jupiter dans le fichier p5
 for(i=0;i<n;i++){
   if(p5!=NULL)
   fprintf(p5," %lf %lf\n",x5[0][i],y5[0][i]);
    else printf("erreur");
}

//ecriture des coordonnees de saturne dans le fichier p6
 for(i=0;i<n;i++){
   if(p6!=NULL)
   fprintf(p6," %lf %lf\n",x6[0][i],y6[0][i]);
    else printf("erreur");
}

//ecriture des coordonnees de uranus dans le fichier p7
 for(i=0;i<n;i++){
   if(p7!=NULL)
   fprintf(p7," %lf %lf\n",x7[0][i],y7[0][i]);
    else printf("erreur");
}

//ecriture des coordonnees de neptune dans le fichier p8
 for(i=0;i<n;i++){
   if(p8!=NULL)
   fprintf(p8," %lf %lf\n",x8[0][i],y8[0][i]);
    else printf("erreur");
}
fclose(p1); //fermeture du fichier p1
fclose(p2); //fermeture du fichier p2
fclose(p3); //fermeture du fichier p3
fclose(p4); //fermeture du fichier p4
fclose(p5); //fermeture du fichier p5
fclose(p6); //fermeture du fichier p6
fclose(p7); //fermeture du fichier p7
fclose(p8); //fermeture du fichier p8
return 0; }






