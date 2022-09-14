/*******************************************************************
                    Functions for prephoton
                        (functions.c) 

 Spectra for bremsstrahlung and pair-production,
 and angular distributions.

*****************************************************************/
double sqr(double x)
{
  return x*x;
}

#define realpower(x,y) pow(x,y)

/* \int_x^\infty dt K_{5/3}(t) */
double int_from_x_to_infty_BesselK5_3(double x)
{
  int nu=2;
  double val;
  double xpower;
  if(x<2)
    {
      val = x*(-pi/sqrt(3) + 2*dbskr3_(&x,&nu) 
	       + 2.53143825*realpower(x,2.0/3.0)*(1+0.09375000*x*x
                                                   +0.00401786*x*x*x*x
                                                   +0.00008789*x*x*x*x*x*x)
               - 1.20910963*realpower(x,4.0/3.0)*(1+0.07500000*x*x
	                                           +0.00251116*x*x*x*x
	                                           +0.00004566*x*x*x*x*x*x)
              );
    }
  else if(x<5)
    {
      val = 1.15223-0.64754*x;
      xpower=x*x;
      val=val+0.12883*xpower;
      xpower=xpower*x;
      val=val-0.00892*xpower;
    }
  else
    {
      nu=1;
      val=2*dbskr3_(&x,&nu)*(1+2*sqr(1.33333333/x));
      nu=2;
      val=val+x*dbskr3_(&x,&nu)*(1-0.5*sqr(1.33333333/x));
    }
  return(val/x);
}

/* Bremmstrahlung function */
double kappa(double x)
{
  int nu=2;
  double val;
  double xpower;
  if(x<2)
    {
      val = x*(-pi/sqrt(3) + 2*dbskr3_(&x,&nu) 
	       + 2.53143825*realpower(x,2.0/3.0)*(1+0.09375000*x*x
                                                   +0.00401786*x*x*x*x
                                                   +0.00008789*x*x*x*x*x*x)
               - 1.20910963*realpower(x,4.0/3.0)*(1+0.07500000*x*x
	                                           +0.00251116*x*x*x*x
	                                           +0.00004566*x*x*x*x*x*x)
              );
    }
  else if(x<5)
    {
      val = 1.15223-0.64754*x;
      xpower=x*x;
      val=val+0.12883*xpower;
      xpower=xpower*x;
      val=val-0.00892*xpower;
    }
  else
    {
      nu=1;
      val=2*dbskr3_(&x,&nu)*(1+2*sqr(1.33333333/x));
      nu=2;
      val=val+x*dbskr3_(&x,&nu)*(1-0.5*sqr(1.33333333/x));
    }
  return(val);
}

/*dP/(dt du), u: outgoing electron's energy fraction 
  (Klepikov)*/
double pair_production_spectrum(double v, double chi)
{
  double alpha=2/(3*(chi*v*(1-v)));
  double val;
  int i;
  i=2;
  val=dbskr3_(&alpha,&i)/(v*(1-v));
  val-=kappa(alpha)/alpha;
  return(val);
}

/*Returns a random energy fraction for the electron, with distribution
given by pair_production_spectrum()*/
double random_electron_energy_fraction(double chi)
{
  double prob_array[1025];
  double sum_prob_array[1025],at_random,eta;
  int i,min,max;
  for (i=0;i<1025;i++)
    {
      eta=(i+0.5)/1025.0;
      prob_array[i]=pair_production_spectrum(eta,chi);
    }
  sum_prob_array[0]=prob_array[0];
  for (i=1;i<1025;i++)
    {
      sum_prob_array[i]=sum_prob_array[i-1]+prob_array[i];
    }
  at_random=urandom_()*sum_prob_array[1024];
  min=0;
  max=1024;
  for(;max!=min+1;)
    {
      if (at_random<sum_prob_array[min+(max-min)/2]){max=min+(max-min)/2;}
      else {min=min+(max-min)/2;}
    }
  return((min + urandom_())/1024.0);
}

/*Erber's analytical aproximation*/
double pair_production_probability(double Energy, 
                                   double chi, double deltatime)
{
  double besselarg,result;
  int mu=1;
  besselarg=2/(3*chi);
  result=1.234E18*deltatime*(electron_mass*c*c/(Energy*1E9*e))*
               sqr(dbskr3_(&besselarg,&mu));
  return(result);
}

double bremsstrahlung_landau(double eta, double gamma)
{
  double arg=2*eta/(3*gamma*(1-eta));
  int mu=2;
  return(eta*int_from_x_to_infty_BesselK5_3(arg)+
         eta*eta*eta/(1-eta)*dbskr3_(&arg,&mu));
}

/*gamma in Erber's review*/
double bremsstrahlung_erber(double eta, double gamma)
{
  double arg=2*eta/(3*gamma*(1-eta));
  return(eta*int_from_x_to_infty_BesselK5_3(arg));
}

int photon_number_erber(double energy, double eta, double chi,
                         double Dt, double Dlogenergy)
{
  int counter=0;
  double aux;
  double number_of_intervals,interval_probability;
  double mean_photon_number=
    alpha*sqr(electron_mass*c*c)/(pi*sqrt(3)*hbar*(energy*1E9*e));
  mean_photon_number*=bremsstrahlung_erber(eta,chi)*Dt*Dlogenergy;
  number_of_intervals=floor(100*mean_photon_number)+1;
  interval_probability=mean_photon_number/number_of_intervals;
  if (mean_photon_number>1)
    {
      printf(" ** adim/E=%E; energy=%E; Dt=%E\n",
                     alpha*sqr(electron_mass*c*c)/(pi*sqrt(3)*hbar*(energy*1E9*e)),
                     energy,Dt);
    }
  for (aux=1;(aux-0.1)<number_of_intervals;aux+=1)
    {
      if (urandom_()<interval_probability){counter++;}
    }
  /*  printf("Prob=%E num=%i\n",interval_probability,counter);*/
  return(counter);
}

int photon_number_landau(double energy, double eta, double chi,
                         double Dt, double Dlogenergy)
{
  int counter=0;
  double aux;
  double number_of_intervals,interval_probability;
  double mean_photon_number=
    alpha*sqr(electron_mass*c*c)/(pi*sqrt(3)*hbar*(energy*1E9*e));
  mean_photon_number*=bremsstrahlung_landau(eta,chi)*Dt*Dlogenergy;
  number_of_intervals=floor(100*mean_photon_number)+1;
  interval_probability=mean_photon_number/number_of_intervals;
  if (mean_photon_number>1)
    {
      printf(" ** adim/E=%E; energy=%E; Dt=%E\n",
                     alpha*sqr(electron_mass*c*c)/(pi*sqrt(3)*hbar*(energy*1E9*e)),
                     energy,Dt);
    }
  for (aux=1;(aux-0.1)<number_of_intervals;aux+=1)
    {
      if (urandom_()<interval_probability){counter++;}
    }
  return(counter);
}

double conversion_probability()
{
  double chi();
  struct particle *auxiliar_photon;
  double prob_non_up_to_now=1;
  double proba, chi_now, dt;
  auxiliar_photon = (struct particle *)malloc(sizeof(struct particle));
  auxiliar_photon->particleID = first_particle->particleID;
  auxiliar_photon->energy = first_particle->energy;
  int comp;
  for(comp=0;comp<3;comp++) 
    auxiliar_photon->position[comp] = first_particle->position[comp];
  auxiliar_photon->altitude = first_particle->altitude;
  for(comp=0;comp<3;comp++) 
    auxiliar_photon->motion_direction[comp] = 
      first_particle->motion_direction[comp];
  auxiliar_photon->last_time = first_particle->last_time;
  auxiliar_photon->next = NULL;
  for (;
       auxiliar_photon->altitude >= atm_injection_altitude;
       auxiliar_photon->last_time+=dt)
    {
      chi_now=chi(auxiliar_photon);
      if (auxiliar_photon->altitude<0.1*earth_radius){dt=5E-7;}
      else {dt=1E-5;}
      proba=pair_production_probability(auxiliar_photon->energy,chi_now,dt);
      prob_non_up_to_now *= (1-proba);
      for(i=0;i<3;i++)
	auxiliar_photon->position[i]+=c*dt*first_particle->motion_direction[i];
    }
  free(auxiliar_photon);
  auxiliar_photon=NULL;
  return(1-prob_non_up_to_now);
}

double french_conversion_probability()
{
  double frenchi();
  struct particle *auxiliar_photon;
  double prob_non_up_to_now=1;
  double proba, chi_now, dt;
  auxiliar_photon = (struct particle *)malloc(sizeof(struct particle));
  auxiliar_photon->particleID = first_particle->particleID;
  auxiliar_photon->energy = first_particle->energy;
  int comp;
  for(comp=0;comp<3;comp++) 
    auxiliar_photon->position[comp] = first_particle->position[comp];
  auxiliar_photon->altitude = first_particle->altitude;
  for(comp=0;comp<3;comp++) 
    auxiliar_photon->motion_direction[comp] = 
      first_particle->motion_direction[comp];
  auxiliar_photon->last_time = first_particle->last_time;
  auxiliar_photon->next = NULL;
  for (;
       auxiliar_photon->altitude >= atm_injection_altitude;
       auxiliar_photon->last_time+=dt)
    {
      chi_now=frenchi(auxiliar_photon);
      if (auxiliar_photon->altitude<0.1*earth_radius){dt=5E-7;}
      else {dt=1E-5;}
      proba=pair_production_probability(auxiliar_photon->energy,chi_now,dt);
      prob_non_up_to_now *= (1-proba);
      for(i=0;i<3;i++)
	auxiliar_photon->position[i]+=c*dt*first_particle->motion_direction[i];
    }
  free(auxiliar_photon);
  auxiliar_photon=NULL;
  return(1-prob_non_up_to_now);
}
