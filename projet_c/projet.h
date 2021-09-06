#ifndef projet_h_         //declaration de la librairie
#define projet_h_


#define G (double)6.67408*pow(10.,-11)          //constante de gravitation
#define pi 3.14                                 //le nombre pi
#define m_terre (double)6.*pow(10.,24)          //masse de la terre
#define m_lune (double)7.36*pow(10.,22)         //masse de la lune
#define m_mercure (double)3.285*pow(10.,23)     //masse de mercure
#define m_venus (double)4.867*pow(10.,24)       //masse de venus
#define m_mars (double)6.39*pow(10.,24)         //masse de mars
#define m_jupiter (double)1.898*pow(10.,27)     //masse de jupiter
#define m_saturne (double)5.683*pow(10.,26)     //masse de saturne
#define m_uranus (double)8.681*pow(10.,25)      //masse de uranus
#define m_neptune (double)1.024*pow(10.,26)     //masse de neptune
#define m_soleil (double)1.989*pow(10.,30)      //masse du soleil

typedef struct {               //structure utilisée pour representer chaque planete
int indice;                    //necessaire pour l''attribution de la masse
double masse;                  //la masse
double x0,y0,vx0,vy0;          //postion et vitesse initiale
}planete;

typedef struct{                //structure utilisée pour l''initialisation de la masse
int indice;
double masse;}m_planete;



double** matrice2n(int a);          //creation d''une matrice 2*n
double** matrice4n(int a);          //creatiion d''une matrice 4*n

double ** position1(void);         //fonction utilisée pour la 1 ere etape
double ** position2(void);         //fonction utilisée pour la 2 eme etape
double ** position3( planete p);   //fonction utilisée pour la 3 eme etape
void **interaction(void);          //fonction utilisée pour la 4 eme etape

planete init_planete();            //fonction utilisée pour l'initialisation des donnees de la planete






# endif          //declaration de la fin de librairie
