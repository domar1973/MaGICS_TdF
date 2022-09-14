void propagate_1D_klepikov()
{
  struct particle *aux,*next;
  double dt,deltat,chi_now;/*dt: time step; Dt: aux's time flight*/
  double proba,energy1,energy2;
  int i;
  if (first_particle->altitude>0.1*earth_radius){dt=1E-5;}
  else                                          {dt=5E-7;} 
  mytime+=dt;
  for (aux=first_particle;aux!=NULL;aux=next)
    {
      next=aux->next;
      chi_now=chi(aux);
      deltat=mytime-(aux->last_time);
      deltat=mytime-aux->last_time;
      if (aux->particleID==1)
	{
	  proba=pair_production_probability(aux->energy,chi_now,deltat);
	  if ((debuglvl>0)&&(proba>0.01)) 
            printf("Pair Prod. Prob=%E\n",proba);
	  if (urandom_()<proba)
	    {
	      if (ThereWasNoConversion==1) ThereWasNoConversion=0;
	      if (FirstInteractionAlready==0)
		{
		  FirstInteractionAlready = 1;
		  int    varn = 3;
		  double vval;
		  vval = aux->altitude;
		  speisetrealvar_(&varn, &vval, &irc);
		}
	      energy1=aux->energy*random_electron_energy_fraction(chi_now);
	      energy2=aux->energy-energy1;
	      add_particle( 2,energy1,
                           aux->position,
                           aux->motion_direction,
			   mytime-urandom_()*deltat,
                           aux->altitude);
	      add_particle(-2,energy2,
                           aux->position,
                           aux->motion_direction,
			   mytime-urandom_()*deltat,
                           aux->altitude);
	      if (debuglvl>0)
		printf("Pair created at h=%E m\n",aux->altitude);
	      next=aux->next;
	      remove_particle(aux);
	    }
	  /*Inject soft photon*/
	  if ((aux!=NULL)&&(aux->energy<5E8)) {aux->altitude=atmospheric_thickness-100;}
	}
      /*If the particle is an electron or positron*/
      if ((aux->particleID==2)||(aux->particleID==-2))
	{
	  /*Next 'for' reserved in case mean_photon_number is too
	    large for the hardest photons*/
	int time_divisor=1;
	int timer;
	for(timer=0;timer<time_divisor;timer++)
	{
	  double eta;
	  int photon_number,i;
	  double energy_slice=0.9;
	  /*Radiated spectrum in energy bands (eta: photon's energy fraction)*/
	  for (eta=1;eta>1E-10;eta*=energy_slice)
	    {
	      photon_number=photon_number_landau(aux->energy,eta*sqrt(energy_slice),
                                       chi_now,deltat,-log(energy_slice));
	      for(i=0;i<photon_number;i++)
		{
		  double this_energy=eta*(aux->energy)*
                                         (energy_slice + urandom_()*(1-energy_slice));
		  add_particle(1,this_energy,aux->position,
		    aux->motion_direction,mytime-urandom_()*deltat,aux->altitude);
		  next=aux->next;
		  aux->energy-=this_energy;
		  total_photon_number++;
		  if(this_energy>5E8) /*Hard photons*/
                    {
		      hard_photons++;
		      if((debuglvl>0)&&(i==photon_number-1)&&(photon_number>1))
		      {printf("step with %i hard photons at eta=%E, time=%E\n",
                                      photon_number,eta,mytime);}
		    }
		}
	    }
	}
	}
      for (i=0;i<3;i++){aux->position[i]+=dt*c*aux->motion_direction[i];}
      if (aux!=NULL){aux->last_time=mytime;}
    }
}
