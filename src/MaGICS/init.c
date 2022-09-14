/*
   Initialization routines.
*/

#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


int  rc;

double getdblglobaldef(char *varname, double defval, int *rc)
{
  /* Returning a given global variable as a double number,
     or a default value if it is not set. */

  char   auxstring[200] = "UUUUU";
  char   auxstr2[200], *auxptr;
  int    auxlen=0, vtype=0;
  double tmp;

  void   getglobalc();

                                 
  /* Trying to read the variable. */

  getglobalc(varname, &vtype, auxstring, &auxlen);
  if (auxlen <= 0) {
    /* Not set. Returning the default value. */
    *rc = 4;
    return(defval);
  }
  else {
    /* The variable seems to exist. Trying to read its value as a number.*/
    auxptr = &auxstr2[0];
    tmp = strtod(auxstring, &auxptr);
    if (strlen(auxptr) <= 0) {
      *rc = 0;
      return(tmp); /* Everything went OK. */
    }
    else {
      *rc = 8;
      return(defval);  /* Bad read. Returning the default value. */
    }
  }
}

int getintglobaldef(char *varname, int defval, int *rc)
{
  /* Returning a given global variable as an integer number,
     or a default value if it is not set. */

  char   auxstring[200];
  int    auxlen, vtype;
  int    tmp;

  void   getglobalc();

  /* Trying to read the variable. */

  *rc = 4; 
  getglobalc(varname, &vtype, auxstring, &auxlen);

  if (auxlen <= 0) {
    /* Not set. Returning the default value. */
    return(defval); 
  }
  else {
    /* The variable seems to exist. Trying to read its value as a number.*/
    if (sscanf(auxstring, "%d", &tmp)) {;
      *rc = 0;
      return(tmp); /* Everything went OK. */
    }
    else {
      *rc = 8;
      return(defval);  /* Bad read. Returning the default value. */
    }
  }
}


int inic()
{
  char *args[20];
  /* General initialization. This routine must be run immediately after
    setting up the AIRES interface.
  */

  int    getintglobaldef();
  double getdblglobaldef();
  int    err=0;


  debuglvl = getintglobaldef("PEMF_Deblevel",0, &rc);
  if (rc == 8) 
    {
      printf("Error 21: cannot read 'PEMF_Deblevel' parameter\n");
      err = 21;
    }

  /* Getting module parameters and environmental quantities. */

  speigetparsc(&parstring, &parstringlen);

  /*Reading arguments and storing them in arg[i] (array of strings)*/
  int nargs = -1;
  if (parstring[0] != ' ') {
    nargs++;
    args[nargs] = &parstring[0];
  }
  for (i = 1; i < parstringlen; i++) {
    if (parstring[i] == ' ') {
      parstring[i] = '\0';
      nargs++;
      while (parstring[i+1]==' ') i++;
      args[nargs] = &parstring[i+1];
    }
  }
  nargs++; /*This si now the number of arguments*/
  if (nargs > 2)
    {
      printf("Warning 11: too much parameters\n");
      err = 11;
    }

  croinputdata0_(&intdata[1],&realdata[1],&shprimcode[1],&paquesirve[1]);


  /*Geographical data*/
  geographical_location[0]= realdata[20];
  geographical_location[1]=-realdata[21];
  /*because latitude here is definer western-positive*/
  site_declination=realdata[24];
  sindec_0=sin(pi/180.0*site_declination);
  cosdec_0=cos(pi/180.0*site_declination);
  sintheta_0=sin(pi/180.0*geographical_location[0]); 
  costheta_0=cos(pi/180.0*geographical_location[0]);
  sinphi_0  =sin(pi/180.0*geographical_location[1]); 
  cosphi_0  =cos(pi/180.0*geographical_location[1]);

  /*Sets Site coordinates (cartesian)*/
  x_0[2]= earth_radius*sintheta_0;
  x_0[1]= earth_radius*costheta_0*cosphi_0;
  x_0[0]= earth_radius*costheta_0*sinphi_0;

  atm_injection_altitude = realdata[9];

  fdate = realdata[26];

  if (debuglvl > 0) {
    printf("PARAMS: %s Length: %i\n",parstring,parstringlen);
    printf("LAT: %g LONG: %g\n",realdata[20],realdata[21]);
    printf("bfield=%E inc=%E dec=%E\n",
                    realdata[22],realdata[23],realdata[24]);
    printf("sintheta_0=%E costheta_0=%E sinphi_0=%E cosphi_0=%E\n",
                    sintheta_0,costheta_0,sinphi_0,cosphi_0);
    printf("geog. location (cartesian): %E %E %E\n",
                    x_0[0],x_0[1],x_0[2]);
  
    }


  /***************************************************************/

  /* Still getting more quantities that come from AIRES.  */

  /*Schwinger Magnetic Field*/
  bcrit=sqr(electron_mass*c)/(e*hbar);

  /* Where are we starting from? */
  /* Setting initial possition (of primary)*/
  double defval=5*earth_radius;
  altitude_0 = getdblglobaldef("PEMF_init_altit",defval,&rc);
  if (rc == 8) 
    {
    printf("Error 22: cannot read 'PEMF_init_altit' parameter\n");
    err = 22;
    }
  if (debuglvl>0) {printf("initial_altitude= %E, rc= %i \n",
			   altitude_0,rc);}

  double k_Cartesian[3]; /*Shower axis in cartesian coords*/
  double lambda,k_CartesianDOTx_0;
  /*This solves for lambda in the vector equation:
  (x_0 - lambda*k_Cartesian) = (earth_radius+altitude_0)^2
  thus -lambda*shower_axis is the initial photon position*/
  rotAIRES2Cartesian(shower_axis,k_Cartesian);
  k_CartesianDOTx_0 = k_Cartesian[0]*x_0[0]
                     +k_Cartesian[1]*x_0[1]
                     +k_Cartesian[2]*x_0[2];
  lambda = sqrt(sqr(k_CartesianDOTx_0) +
                altitude_0*(2*earth_radius+altitude_0) )
           +k_CartesianDOTx_0; 
  /*Setting initial photon position*/
  int i;
  for(i=0;i<3;i++){
    initial_position[i] = default_injection_position[i] -
                          lambda*shower_axis[i];}
  if (debuglvl > 0){
    printf("lambda=%E init pos: %E %E %E\n",
                   lambda,initial_position[0],
                   initial_position[1],
	           initial_position[2]); }

  /* Do we convert all photons? */
  AllConverted = 0;

  /* Prototypes of available algorithms. */

  void propagate_1D_erberaprox();
  void propagate_1D_erberaprox_v2();
  void propagate_1D_klepikov();
  void propagate_1D_klepikov_v2();    /*Default*/
  void propagate_3D_klepikov();
  void propagate_3D_klepikov_v2();
  void propagate_test(); /*The particle moves without interaction*/

  /* Reading algorithm in use */
  if (nargs==0)
    {
      propagate_particle = propagate_1D_klepikov_v2;
      if (debuglvl > 0) {
         printf("Using default 1D_klepikov_v2, photon\n");
      }
    }
if (nargs>0)
{
  if (strcmp(args[0],"1Derberaprox")==0)
    {
      propagate_particle = propagate_1D_erberaprox;
      if (debuglvl > 0) {
         printf("Using 1D_erberaprox, photon\n");
      }
    }
  else if (strcmp(args[0],"1Derberaprox+")==0)
    {
      propagate_particle = propagate_1D_erberaprox;
      PrimaryParticle=2;
      if (debuglvl > 0) {
         printf("Using 1D_erberaprox, positron\n");
      }
    }
  else if (strcmp(args[0],"1Derberaprox-")==0)
    {
      propagate_particle = propagate_1D_erberaprox;
      PrimaryParticle=-2;
      if (debuglvl > 0) {
         printf("Using 1D_erberaprox, electron\n");
      }
    }
  else if(strcmp(args[0],"1Derberaprox_v2")==0)  
    {          
      propagate_particle = propagate_1D_erberaprox;
      if (debuglvl > 0) {      
	printf("Using 1D_erberaprox_v2, photon\n");
      }
    }
  else if(strcmp(args[0],"1Derberaprox_v2+")==0)  
    {          
      propagate_particle = propagate_1D_erberaprox;
      PrimaryParticle=2;
      if (debuglvl > 0) {      
	printf("Using 1D_erberaprox_v2, positron\n");
      }
    }
  else if(strcmp(args[0],"1Derberaprox_v2-")==0)  
    {          
      propagate_particle = propagate_1D_erberaprox;
      PrimaryParticle=-2;
      if (debuglvl > 0) {      
	printf("Using 1D_erberaprox_v2, electron\n");
      }
    }
  else if(strcmp(args[0],"1Dklepikov")==0)  
    {          
      propagate_particle = propagate_1D_klepikov;
      if (debuglvl > 0) {      
	printf("Using 1D_klepikov, photon\n");
      }
    }
  else if(strcmp(args[0],"1Dklepikov+")==0)  
    {          
      propagate_particle = propagate_1D_klepikov;
      PrimaryParticle=2;
      if (debuglvl > 0) {      
	printf("Using 1D_klepikov, positron\n");
      }
    }
  else if(strcmp(args[0],"1Dklepikov-")==0)  
    {          
      propagate_particle = propagate_1D_klepikov;
      PrimaryParticle=-2;
      if (debuglvl > 0) {      
	printf("Using 1D_klepikov, elecron\n");
      }
    }
  else if(strcmp(args[0],"1Dklepikov_v2")==0)  
    {
      propagate_particle = propagate_1D_erberaprox;
      if (debuglvl > 0) {      
	printf("Using 1D_klepikov_v2, photon\n");
      }
    }
  else if(strcmp(args[0],"1Dklepikov_v2+")==0)  
    {
      propagate_particle = propagate_1D_erberaprox;
      PrimaryParticle=2;
     if (debuglvl > 0) {      
	printf("Using 1D_klepikov_v2, positron\n");
      }
    }
  else if(strcmp(args[0],"1Dklepikov_v2-")==0)  
    {
      propagate_particle = propagate_1D_erberaprox;
      PrimaryParticle=-2;
     if (debuglvl > 0) {      
	printf("Using 1D_klepikov_v2, electron\n");
      }
    }
  else if(strcmp(args[0],"3Dklepikov")==0)  
    {
      propagate_particle = propagate_3D_klepikov;
      if (debuglvl > 0) {      
	printf("Using 3D_klepikov, photon\n");
      }
    }
  else if(strcmp(args[0],"3Dklepikov+")==0)  
    {
      propagate_particle = propagate_3D_klepikov;
      PrimaryParticle=2;
      if (debuglvl > 0) {      
	printf("Using 3D_klepikov, positron\n");
      }
    }
  else if(strcmp(args[0],"3Dklepikov-")==0)  
    {
      propagate_particle = propagate_3D_klepikov;
      PrimaryParticle=-2;
      if (debuglvl > 0) {      
	printf("Using 3D_klepikov, electron\n");
      }
    }
  else if(strcmp(args[0],"3Dklepikov_v2")==0)  
    {
      propagate_particle = propagate_3D_klepikov_v2;
      if (debuglvl > 0) {      
	printf("Using 3D_klepikov_v2, photon\n");
      }
    }
  else if(strcmp(args[0],"3Dklepikov_v2+")==0)  
    {
      propagate_particle = propagate_3D_klepikov_v2;
      PrimaryParticle=2;
      if (debuglvl > 0) {      
	printf("Using 3D_klepikov_v2, positron\n");
      }
    }
  else if(strcmp(args[0],"3Dklepikov_v2-")==0)  
    {
      propagate_particle = propagate_3D_klepikov_v2;
      PrimaryParticle=-2;
      if (debuglvl > 0) {      
	printf("Using 3D_klepikov_v2, electron\n");
      }
    }
  else if(strcmp(args[0],"UnconvertedPhoton")==0)  
    {
      propagate_particle = propagate_test;
      if (debuglvl > 0) {      
	printf("Using UnconvertedPhoton");
      }
    }
  else
    {
      propagate_particle = propagate_1D_klepikov_v2;
      if (strcmp(args[0],"")!=0)
	{      
	  printf("Error 23: bad propagation algorythm, using default\n");
	  err = 23;
	}
      if (debuglvl > 0) {
        printf("nargs=%i, args[0]=%s\n",nargs,args[0]);
      }
    }
}
if (nargs>1) {
if (strcmp(args[0],"UnconvertedPhoton")!=0)
{
  if (strcmp(args[1],"ForceConversion")==0)
    {
      AllConverted=1;
      if (debuglvl > 0) printf("Forcing Conversion\n");
    }
  else if(strcmp(args[1],"NoForceConversion")==0)
    {
      AllConverted=0;
    }
  else  
    {
      printf("Error 24: bad second argument, using default NoForceConversion\n");
      err = 24;
    }
}
else if (strcmp(args[1],"ForceConversion")==0)
{
  printf("Error 25: cannot force conversion");
  err = 25;
}
else if (strcmp(args[1],"NoForceConversion")!=0)
{
  printf("Error 24: bad second argument, using default NoForceConversion\n");
  err = 24;
 }}


for (i=0;i<nargs;i++) printf("# argument %i: %s\n",i,args[i]);
printf("#\n");

return(err);
}

