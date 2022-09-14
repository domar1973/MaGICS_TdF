void propagate_test()
{
  struct particle *aux;
  for(aux=first_particle;aux!=NULL;aux=aux->next)
    {
      int i;
      for(i=0;i<3;i++){aux->position[i]+=1E-6*c*aux->motion_direction[i];}
      chi(aux);
		  if (ThereWasNoConversion==1)
		    {
		      ThereWasNoConversion=0;
		      int    varn = 3;
		      double vval;
		      vval = aux->altitude;
		      speisetrealvar_(&varn, &vval, &irc);
		    }
    }
  mytime+=1E-6;
}

void propagate_3D_klepikov()
{
		  if (ThereWasNoConversion==1) ThereWasNoConversion=0;
}

void propagate_3D_klepikov_v2()
{
  double dt,timeprod,time_left,chi_now;
  double proba,energy1,energy2;
  int i,k,ready=0,photons_at_turn=0;
  double probarray[1025]; /*bremsstrahlung, by energy bands...*/
  double sumprobarray[1025] /*...and its 'integral'*/;
  double eta,P,at_random,Adim;
  double emission_time,v,energy_slice=0.9;
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
              if (proba>0.01) 
		{
		  printf("  ****Pair Prod. Prob=%E****\n",proba);
		  
		}
              if (urandom_()<proba)
	        {
		  if (ThereWasNoConversion==1)  ThereWasNoConversion=0;
		  if (FirstInteractionAlready==0)
		    {
		      FirstInteractionAlready = 1;
		      int    varn = 3;
		      double vval;
		      vval = first_particle->altitude;
		      speisetrealvar_(&varn, &vval, &irc);
		    }
		  v=random_electron_energy_fraction(chi_now);
		  if ((v<0)||(v>1))
		    {
		      printf("v=%E !!!\n",v);
		      
		    }
		  energy1=(first_particle->energy)*v;
		  energy2=(first_particle->energy)-energy1;
		  timeprod=mytime+urandom_()*dt;
	          add_particle( 2,energy1,
                               first_particle->position,
                               first_particle->motion_direction,
			       timeprod,
                               first_particle->altitude);
		  add_particle(-2,energy2,
                               first_particle->position,
                               first_particle->motion_direction,
			       timeprod,
                               first_particle->altitude);
		  printf("Pair created at h=%E m, v=%E\n",
                                 first_particle->altitude,v);
		  remove_particle(first_particle);
		  ready=1;
		}
	      if (ready==0)
		{		
		  if ((first_particle->altitude)<atmospheric_thickness){ready=1;}
		  if (ready==0) 
		    {
		      for(i=0;i<3;i++)
			{
			  first_particle->position[i]+=
			    c*dt*first_particle->motion_direction[i];
			}
		    }
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
		  probarray[i] = Adim*bremsstrahlung_landau(pow(energy_slice,i)*
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
		  double omega=sqr(electron_mass)*c*c/(first_particle->energy*hbar)*
		                lastbfield.strength;
		  double spd[3];
		  for (i=0;i<3;i++)
		    {
		      spd[i]=c*first_particle->motion_direction[i];
		    }
		  spd[0]+=omega*c*emission_time*
			( first_particle->motion_direction[1]-lastbfield.direction[2]);
		  spd[1]+=omega*c*emission_time*
			(-first_particle->motion_direction[0]+lastbfield.direction[2]);
		  spd[2]+=omega*c*emission_time*
			( first_particle->motion_direction[0]-lastbfield.direction[1]);
		  first_particle->last_time+=emission_time;
		  double v = sqrt( sqr(spd[0])+sqr(spd[1])+sqr(spd[2]) );
		  for (i=0;i<3;i++)
		    {
		      first_particle->position[i] += emission_time*spd[i];
		      first_particle->motion_direction[i] = spd[i]/v;
		    }
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
  printf("** %i emited photons\n",photons_at_turn);
}
