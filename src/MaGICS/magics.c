/*******************************************************************

          Magnetic Gamma Interactions Computer Simulator
	                     (MaGICS_TdF)
                Photon in Earth's Magnetic Field
              For use as 'special primary' at AIRES
                         March     2003
                         June      2003
                         September 2003
                         April     2004
                         --------------
                         September 2022

Beta release September 9 2022

writen by Daniel Badagnani 
(dbadagnani@untdf.edu.ar)
Universidad Nacional de Tierra del Fuego
Ushuaia - ARGENTINA

****************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <constants.h>
#include <vars.h>
#include <prototypes.h>
#include <particle_list.c>
#include <functions.c>
#include <geometry.c>
#include <1D_erber.c>
#include <1D_erber_v2.c>
#include <1D_klepikov.c>
#include <1D_klepikov_v2.c>
#include <others.c>
#include <init.c>

int main()
{

  int irc = 215;
  int jrc;
  /* MODULE VERSION */

  int mversion0[3] = {1,5,0};
 
  int mversion, oldversion;
  int inicdone=0;
  int tswitch=0;
  int coord_sys=0;
  int TrialNumber=0;
  double peso=1;
  struct particle *this_particle;

  void speistart_(), speisetfltvar_();
  int  inic();

  /*
    Invoking AIRES interface, and general initialization.
  */


  speistart_(&shower_number, &primary_energy,
             &default_injection_position[0], &injection_depth,
             &ground_altitude, &ground_depth,
             &d_ground_inj, &shower_axis[0]);

for (;ThereWasNoConversion==1;TrialNumber++){

  if (inicdone==0)
    {
      mversion = 10000 * mversion0[0] + 100 * mversion0[1] + mversion0[2];
      printf("#\n");
      printf("# Preshower %i.%i.%i\n",mversion0[0],mversion0[1],mversion0[2]);
      printf("#\n");
      printf("# Daniel Badagnani\n");
      printf("# Universidad Nacional de La Plata\n");
      printf("# ARGENTINA\n");
      printf("#\n");
      irc = inic();
      inicdone++;
    }

  if (AllConverted==0) ThereWasNoConversion=0;

  if (debuglvl > 0) {
    printf("injection: %E %E %E\n",
                   default_injection_position[0],
                   default_injection_position[1],
                   default_injection_position[2]);
    printf("shower_axis: %E %E %E\n",
                   shower_axis[0],shower_axis[1],shower_axis[2]);
  }

  /* Starting the simulations. Setting particle list as empty, and 
     injecting primary.
  */
  first_particle=NULL;
  mytime=0;

  add_particle(PrimaryParticle,primary_energy,initial_position,
               shower_axis,mytime,5*earth_radius);

  double prob;
  prob=conversion_probability();
  /* Sendin conversion probability to AIRES kernel*/
  int    varn = 1;
  double vval;
  vval = prob;
  /*speisetfltvar_(&varn, &vval, &jrc);*/
  speisetrealvar_(&varn, &vval, &jrc);

  /*Simulation begins.
    *propagate_particle is either a xxx or xxx_v2 routine
    (set at inic.c) 
    Each call to propagate_xxx is a time step;
    Each call to propagate_xxx_v2 propagates full history for one particle*/

  while(first_particle!=NULL) {

    (*propagate_particle)();

      /*Send to AIRES if got to atmosphere*/
      struct particle *next;
      double timezero=0;
      spinjpoint_(&coord_sys,&(default_injection_position[0]),
                  &(default_injection_position[1]), 
                  &(default_injection_position[2]),
                  &tswitch,&timezero,&jrc);
      for(this_particle=first_particle;
          this_particle!=NULL;
          this_particle=next)
	{
	  next=this_particle->next;
	  if((this_particle->altitude < atm_injection_altitude)&&
             (ThereWasNoConversion==0))
	    {
	      if(debuglvl>0)
	        printf("ID=%i; Energy: %E; irc(set.inj.point)=%i;",
                               this_particle->particleID,
                               this_particle->energy,jrc);
	      spaddp0_(&(this_particle->particleID),
                       &(this_particle->energy),&(coord_sys),
                       &(this_particle->motion_direction[0]),
                       &(this_particle->motion_direction[1]),
                       &(this_particle->motion_direction[2]),&peso,&jrc);
	      if(debuglvl>0){
		printf(" irc(inj)=%i\n",jrc);
	      }
	      remove_particle(this_particle);
	    }
	  else if (this_particle->altitude < atm_injection_altitude)
	    {
	      remove_particle(this_particle);
	    }
	}
  }
  if (debuglvl>0)
  printf("Tot. Phot. Num: %i  Hard: %i\n",total_photon_number,hard_photons);

/*End of simulation*/}

  if (debuglvl>0)
    printf("number of trials: %i; ThereWasNoConversion: %i\n",
           TrialNumber,ThereWasNoConversion);

  /* Sending the number of trials to AIRES kernel. 
     If conversion is forced then returns zero.*/
  int    varn = 2;
  double vval; 
  vval = TrialNumber;
  if (AllConverted==0) vval=0;
  speisetrealvar_(&varn, &vval, &jrc);
  /* Recording the variables */
  varn = 1;
  sprecordvars_(&varn);

  speimv_(&mversion, &oldversion); /* Saving module version. */

  speiend_(&irc);
  return 0;
}
