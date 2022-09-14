/*Propagating one particle at a time until it disapear or enter atmosphere*/
void propagate_1D_erberaprox_v2()
{
  double dt,chi_now,timeprod;
  double proba;
  int i,k,ready=0,photons_at_turn=0;
  double probarray[1025]; /*bremsstrahlung, by energy bands...*/
  double sumprobarray[1025] /*...and its 'integral'*/;
  double eta,P,at_random,Adim;
  double emission_time,time_left,energy_slice=0.9;
  int min,max;
  if (first_particle->particleID==1)
    {
      if (first_particle->energy<5E8) 
	{
	  first_particle->altitude=atmospheric_thickness-100;
	  ready=1;
	}
      else
	{
          for (mytime=first_particle->last_time;
               ready==0;
               first_particle->last_time+=dt)
	    {
	      chi_now=chi(first_particle);
	      if (first_particle->altitude<0.1*earth_radius){dt=5E-7;}
	      else {dt=1E-5;}
              proba=pair_production_probability(first_particle->energy,chi_now,dt);
              if ((debuglvl>0)&&(proba>0.01)) 
		{
		  printf("  ****Pair Prod. Prob=%E****\n",proba);
		  
		}
              if (urandom_()<proba)
	        {
		  if (ThereWasNoConversion==1) ThereWasNoConversion=0;
		  if (FirstInteractionAlready==0)
		    {
		      FirstInteractionAlready = 1;
		      int    varn = 3;
		      double vval;
		      vval = first_particle->altitude;
		      speisetrealvar_(&varn, &vval, &irc);
		    }
                  timeprod=mytime+urandom_()*dt;
	          add_particle( 2,first_particle->energy/2,
                               first_particle->position,
                               first_particle->motion_direction,
			       timeprod,
                               first_particle->altitude);
		  add_particle(-2,first_particle->energy/2,
                               first_particle->position,
                               first_particle->motion_direction,
			       timeprod,
                               first_particle->altitude);
		  if (debuglvl>0)
		  printf("Pair created at h=%E m\n",first_particle->altitude);
		  remove_particle(first_particle);
		  ready=1;
		}
	      if (ready==0)
		{		
		  if ((first_particle->altitude)<atmospheric_thickness){ready=1;}
		  if (ready==0) {for(i=0;i<3;i++)
                         {first_particle->position[i]+=
                                      c*dt*first_particle->motion_direction[i];};}
		}
		  
	    }
	}
    }
  /*If the particle is an electron or positron*/
  else /*if ((first_particle->particleID==2)||(first_particle->particleID==-2))*/
    {
      for (;(first_particle->altitude)>=atmospheric_thickness;)
	{
	  Adim=alpha*sqr(electron_mass*c*c)/(pi*sqrt(3)*
                 hbar*(first_particle->energy*1E9*e));
	  chi_now = chi(first_particle);
	  if ((first_particle->altitude)>=(atmospheric_thickness))
	    {
	      for (i=0;i<1025;i++)
		{
		  probarray[i] = Adim*bremsstrahlung_erber(pow(energy_slice,i)*
                                     sqr(energy_slice),chi_now)*(-log(0.9));
		}
	      sumprobarray[0]=probarray[0];
	      for (i=1;i<1025;i++)
		{
		  sumprobarray[i] = sumprobarray[i-1]+probarray[i];
		}
	      dt = 1/(sumprobarray[1024]*100);
	      P=1;
	      time_left=0;
	      for (k=0;k<3;k++) 
		{
		  time_left+=sqr(first_particle->position[k]-default_injection_position[k]);
		}
	      time_left=sqrt(time_left)/c;
	      at_random=urandom_();
		for (i=1;at_random>(1-P);i++)
		  {
		    P*=(1.0 - 1.0/100.0);
		  }
	      emission_time = i*dt;
	      if (emission_time<time_left)
		{
		  at_random=urandom_()*sumprobarray[1024];
		  min=0;
		  max=1024;
		  for (;max!=min+1;)
		    {
		      if (at_random<sumprobarray[min+(max-min)/2]){max=min+(max-min)/2;}
		      else {min=min+(max-min)/2;}
		    }
		  eta = pow(energy_slice,max)*(1 - urandom_()*(1-energy_slice));
		  double this_energy=eta*(first_particle->energy);
		  first_particle->energy-=this_energy;
		  for (i=0;i<3;i++)
		    {
		      first_particle->position[i]+=c*emission_time*
                                             first_particle->motion_direction[i];
		    }
		  first_particle->last_time+=emission_time;
		  add_particle(1,this_energy,first_particle->position,
		       first_particle->motion_direction,first_particle->last_time,
		       first_particle->altitude);
		  photons_at_turn++;
		}
	      else
		{
		  for (i=0;i<3;i++)
		    {
		      first_particle->position[i]+=c*(time_left+1E-7)*
                                             first_particle->motion_direction[i];
		    }
		  first_particle->last_time+=time_left+1E-7;
		}
	    }
	}
    }
  if (debuglvl>0) {
    printf("** %i emited photons\n",photons_at_turn);
  }
}
