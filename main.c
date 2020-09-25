 /*
                      Solving the TOV equations
 Based on chapters II-IV "Neutron Stars for Undergraduates" by R.R. Silbar & S.Reddy

*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double dMdx(double x, double P, double ogamma, double bb )
{
    return (bb*pow(x,2.)*pow(P,ogamma));
}
double dPdx(double x, double P, double M, double ogamma, double a )
{
    return ((-a*pow(P,ogamma)*M)/pow(x,2.));
}


int main()
{
    FILE *f = fopen("NSPR-6.txt", "w");
if (f == NULL)
{
    printf("Error opening file!\n");
    exit(1);
}
FILE *g = fopen("NSMR-6.txt", "w");
if (g == NULL)
{
    printf("Error opening file!\n");
    exit(1);
}

    double ZoA, m,  M, a;
    double Kappa;
    double P;

    double h, x;

    double k1, k2, k3, k4;
    double l1, l2, l3, l4;

//Constants
    double hbar=1.05457e-27;         //erg*s
    double c=2.99792e10;             // cm/s


    double me=9.109389e-28;          // g
    double mN=1.674928e-24;          // g

    double R0=1.477 ;                // Schwartzschild half-radius in km
    double Ms=1.989e33;              // Solar mass in g

    double gamma, grel, gnrel;
    grel = 4./3.;
    gnrel = 5./3.;

// Choose either relativistic or Non-Relativistic regime
    gamma = gnrel;

// Boundary conditions at the center of the star **/

    P = 1e-6;                       // pressure at the center (cgs units)
    h = 1e-5;                       // integration step
    x = 1e-5;                       // assume as the center of the star (for x=0 the equations diverge)
    M = 0.;                         // mass at the center


// Select the type of star
    char star;
    printf(" Select 'W' for white dwarf or 'N' for neutron star. \n");
    scanf("%c", &star);

    if (star == 'W')
    {
            ZoA = 0.5;              // atomic over mass number
            m = me;                 // reduced mass that enters to the polytrope coefficient
            printf("\n You selected White Dwarf. \n");

    }
    if (star == 'N')
    {
            ZoA = 1;                // atomic over mass number
            m = mN;                 // mass that enters to the polytrope coefficient
            printf("\n You selected Neutron Star. \n");
     }


  double kr = (3.*pow(M_PI,2.)*ZoA)/(mN*pow(c,2.));
  printf("kr = %lf", kr);

  if(gamma==gnrel)
  {
        printf("\n You are in the Non-Relativistic regime.");
 double Knrel = (pow(hbar,2.)/(15.*pow(M_PI,2.)*m))*pow(kr,gamma);  // [cm^2 * erg^(-2/3)]
        Kappa = Knrel;
  }
  else
    {
        printf("\n You are in the Relativistc Regime.");
 double Krel =  ((hbar*c)/(12.*pow(M_PI,2.)))*pow(kr,gamma);           // [cm * erg^(1/3)]
        Kappa = Krel;

  }
    printf ("\n gamma = %lf \n", gamma);
    printf("\n Enter a value for a = ");
    scanf("%lf", &a);

    printf("K = %e \n", Kappa);

//Calculate å0
    double gamma1 = gamma - 1. ;

    double eeps0 = (1./Kappa)*pow((R0/a),gamma);
    double ggamma=1./(gamma1);
    double eps0 = pow(eeps0,ggamma);                                     // [erg * cm^(-3)]
    printf("eps0 = %e \n", eps0);



// KBar calculation

    double KBar = Kappa*pow(eps0, gamma1);
    printf("\n \n KBar = %lf \n", KBar);


// Calculate â

    double ogamma= 1./gamma;
    double bb = (4.*M_PI*eps0*pow(10.,15.))/(Ms*pow(c,2.)*pow(KBar,ogamma)); //10^15 is used to regain the units of â in 1/km^3
    printf("b = %lf \n", bb);



/** Integration of TOV equations **/
/*   Runge - Kutta 4th order method
*/

      while (P>= 0)
    {
        k1 = h*dMdx(x, P, ogamma, bb);
        l1 = h*dPdx(x, P, M, ogamma, a);

        k2 = h*dMdx(x + 5e-1*h, P + 5e-1*k1, ogamma, bb);
        l2 = h*dPdx(x + 5e-1*h, P + 5e-1*l1, M, ogamma, a);


        k3 = h*dMdx(x + 5e-1*h, P + 5e-1*k2, ogamma, bb);
        l3 = h*dPdx(x + 5e-1*h, P + 5e-1*l2, M, ogamma, a);


        k4 = h*dMdx(x + h, P + k3, ogamma, bb);
        l4 = h*dPdx(x + h, P + l3, M, ogamma, a);



        M = M + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4);
        P = P + (1.0/6.0)*(l1 + 2*l2 + 2*l3 + l4);

        x = x + h;

      fprintf(f, "%lf  %e   \n", x, P);
      fprintf(g, "%lf  %lf   \n", x, M);

      printf("\n P = %e, R = %lf, M = %lf ", P, x, M);

     }
//      printf("\n P = %e, R = %lf, M = %lf ", P, x, M);




fclose(f);
fclose(g);
    //  system("pause");
    return 0;
}


