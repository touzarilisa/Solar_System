#include"projet.h"  //inclusion du fichier header
//inclusion des bibliotheques standards
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

//programme principal
 int main(){
FILE* fp1;   //fichier texte pour enregistrer les coordonnees de la premiere planete
FILE* fp2;   //fichier texte pour enregistrer les coordonnees de la deuxieme planete
int i,n;
 //les deux planetes
 planete p;

int etape;  //etape du projet a realiser
printf("quelle etape voulez vous realiser:\n1.etape 1\n2.etape 2\n3.etape 3\n4.etape 4\n");  //choix de l''etape a realiser
scanf(" %d",&etape);

//etape 1
if(etape==1){
n=100000;

fp1=fopen("resultatplanete1.txt", "w");      //ouverture du fichier fp1 en mode ecriture
fp2=fopen("resultatplanete2.txt", "w");      //ouverture du fichier fp2 en mode ecriture


double** resultat =matrice4n(n);             //creation de la matrice resultat pour recuperer les donnees de la fct position 1

resultat=position1() ;                       //appel de la fct position 1

printf(" le resultat est");
printf(" les position de la lune");

for(i=0;i<n;i++){                            //ecriture des coordonnes de la lune dans le fichier fp1
printf(" %lf %lf\n",resultat[0][i],resultat[1][i]);
if(fp1!=NULL)
fprintf(fp1," %lf %lf\n",resultat[0][i],resultat[1][i]);
else printf(" erreur");
}
fclose(fp1);

printf(" les position de la terre");
printf(" \n");


for(i=0;i<n;i++){                             //ecriture des coordonnes de la terre dans le fichier fp2
   printf(" %lf %lf\n",resultat[2][i],resultat[3][i]);
   if(fp2!=NULL)
    fprintf(fp2," %lf %lf\n",resultat[2][i],resultat[3][i]);
    else printf("erreur");
}

fclose(fp2);}


//etape 2
if(etape==2){
n=1000000;

fp1=fopen("resultatplanete1.txt", "w");            //ouverture du fichier fp1 en mode ecriture
fp2=fopen("resultatplanete2.txt", "w");            //ouverture du fichier fp2 en mode ecriture


double** resultat = matrice4n(n);                  //creation de la matrice resultat pour recuperer les donnees de la fct position 2

resultat=position2() ;                             //appel de la fct position 2

printf(" le resultat est");
printf(" les position de la lune");

for(i=0;i<n;i++){                                  //ecriture des coordonnes de la lune dans le fichier fp1
printf(" %lf %lf\n",resultat[0][i],resultat[1][i]);
if(fp1!=NULL)
fprintf(fp1," %lf %lf\n",resultat[0][i],resultat[1][i]);
else printf(" erreur");
}
fclose(fp1);

printf(" les position de la terre");
printf(" \n");


for(i=0;i<n;i++){                                 //ecriture des coordonnes de la terre dans le fichier fp2
   printf(" %lf %lf\n",resultat[2][i],resultat[3][i]);
   if(fp2!=NULL)
    fprintf(fp2," %lf %lf\n",resultat[2][i],resultat[3][i]);
    else printf("erreur");
}

fclose(fp2);}


//etape 3
if(etape==3){
n=100000;
p=init_planete();         //initialisation de la masse ,position et vitesse initiale de la planete choisie
printf(" %d %lf %lf %lf %lf %lf\n",p.indice,p.x0,p.y0,p.vx0,p.vy0,p.masse);

fp1=fopen("resultatplanete1.txt", "w");            //ouverture du fichier fp1 en mode ecriture
fp2=fopen("resultatplanete2.txt", "w");            //ouverture du fichier fp2 en mode ecriture


double** resultat =matrice4n(n);                   //creation de la matrice resultat pour recuperer les donnees de la fct position 3
resultat=position3(p) ;                            //appel de la fct position 3

printf(" le resultat est");
printf(" les position de la premiere planete");

for(i=0;i<n;i++){                                  //ecriture des coordonnes de la planete dans le fichier fp1
printf(" %lf %lf\n",resultat[0][i],resultat[1][i]);
if(fp1!=NULL)
fprintf(fp1," %lf %lf\n",resultat[0][i],resultat[1][i]);
else printf(" erreur");
}
fclose(fp1);

printf(" les position du soleil");
printf(" \n");


for(i=0;i<n;i++){                                   //ecriture des coordonnes du soleil dans le fichier fp2
   printf(" %lf %lf\n",resultat[2][i],resultat[3][i]);
   if(fp2!=NULL)
    fprintf(fp2," %lf %lf\n",resultat[2][i],resultat[3][i]);
    else printf("erreur");
}

fclose(fp2);}

//etape 4
if(etape==4){
    interaction();      //appel de la fct interaction
}

return 0;}


