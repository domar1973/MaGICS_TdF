/* Calculation of \sqrt(x+1)-1 */
double heighfactor(double x)
{
  double aux[6];
  int i;
  if(x<0.6)
    {
      aux[0]=x;
      for (i=1;i<5;i++){aux[i]=aux[i-1]*x;}
      aux[5]=  aux[0]/2.0
             - aux[1]/8.0
             + aux[2]/16.0
             - aux[3]*5.0/128.0
             + aux[4]*7.0/256.0;
    }
  else {aux[5]=sqrt(1+x)-1;}
  return(aux[5]);
}


   /* spheric_coords ******************************************** 
   Next function converts cartesian coordinates into spheric ones.
   WARNING: the longitude is defined opposite to AIRES convention
   Cartesian coordinates:
      Origin: center of Earth
      z axis: Earth's rotation axis (+: north)
      y axis: from origin to 0 deg lat, 0 deg long
      x axis: from origin to West
   Spheric coordinates:
      r    : the obvious one
      theta: latitude (positive on northern hemisphere)
      phi  : longitude (positive on western hemisphere) 
   The unit for x, y, z are just the same that those of r. */
void spheric_coords(double x, double y, double z,
                    double *r, double *theta, double *phi,
                    double *sintheta,double *costheta,
                    double *sinphi, double *cosphi)
{
  double xplusysquared,zsquared,rsquared,rcostheta;
  zsquared      = z*z;
  xplusysquared = x*x + y*y;
  rsquared      = xplusysquared + zsquared;
  *sintheta     = sqrt(zsquared/rsquared);
  if (z<0) {*sintheta=-*sintheta;}
  *costheta     = sqrt(xplusysquared/rsquared);
  *r            = sqrt(rsquared);
  rcostheta     = (*r)*(*costheta);
  *sinphi       = x/rcostheta;
  *cosphi       = y/rcostheta;
  *theta        = 180.0/pi*asin(*sintheta);
  if (cosphi>0) {*phi = 180.0/pi*      asin(*sinphi);}
  else          {*phi = 180.0/pi*(pi - asin(*sinphi));};
}



    /* cartesian_magnetic_field *****************************************
    inc = inclination, dec = declination (bfield angles in tangent plane) 
    This function returns B in cartesian (geocentric) coords */
void cartesian_magnetic_field(double latitude, double longitude, double altitude,
                              double sintheta, double costheta, 
                              double sinphi, double cosphi,
                              double *bx, double *by, double *bz)
{
  double b, inc, dec, sininc, cosinc, sindec, cosdec;
  /*b: field strength; inc: inclination; dec: declination*/
  void   geomagnetic_();
  double theta, phi;
  theta=latitude;
  phi=-longitude;
  geomagnetic_(&theta,&phi,&altitude,&fdate,
	       &b,&inc,&dec);
  sininc = sin(pi*inc/180);
  cosinc = cos(pi*inc/180);
  sindec = sin(pi*dec/180);
  cosdec = cos(pi*dec/180);

  /*double sintheta, costheta,sinphi,cosphi;
  sintheta = sin(pi*latitude/180);
  costheta = cos(pi*latitude/180);
  sinphi   = sin(pi*longitude/180);
  cosphi   = cos(pi*longitude/180);*/

  *bx = b*(-sintheta*sinphi*(cosinc*cosdec) + cosphi*(cosinc*sindec) - costheta*sinphi*(sininc));
  *by = b*(-sintheta*cosphi*(cosinc*cosdec) - sinphi*(cosinc*sindec) - costheta*cosphi*(sininc));
  *bz = b*( costheta*       (cosinc*cosdec)                          - sintheta*       (sininc));
}


   /* Rotation from AIRES to cartesian coordinates (SIGN BUG CORRECTED)
   vector_AIRES =     [V_theta, V_phi, V_r]
   vector_Cartesian = [C_x, C_y, C_z]
   (only ROTATION. Origin remains at site)*/
void rotAIRES2Cartesian(double *vector_AIRES, double *vector_Cartesian)
{
  double north_pointing[3];
  /*Rotating from 'local magnetic' North to geographical North*/
  north_pointing[2]=vector_AIRES[2];
  north_pointing[1]=  -vector_AIRES[0]*sindec_0 + vector_AIRES[1]*cosdec_0;
  north_pointing[0]=   vector_AIRES[0]*cosdec_0 + vector_AIRES[1]*sindec_0;
  /*Rotating to Cartesian coords (as defined at comment in spheric_coords)*/
  vector_Cartesian[0]=-north_pointing[0]*sintheta_0*sinphi_0 
                      +north_pointing[1]*cosphi_0 
                      +north_pointing[2]*sinphi_0*costheta_0;

  vector_Cartesian[1]=-north_pointing[0]*sintheta_0*cosphi_0 
                      -north_pointing[1]*sinphi_0 
                      +north_pointing[2]*cosphi_0*costheta_0;

  vector_Cartesian[2]= north_pointing[2]*sintheta_0 
                      +north_pointing[0]*costheta_0;
 
}

void rotCartesian2AIRES(double *vector_Cartesian, double *vector_AIRES)
{
  double local_north_pointing[3];
  /*Rotating from Cartesian coords (as defined at comment in spheric_coords)*/
  local_north_pointing[0]=-vector_Cartesian[0]*sintheta_0*sinphi_0 
                          -vector_Cartesian[1]*sintheta_0*cosphi_0 
                          +vector_Cartesian[2]*costheta_0;

  local_north_pointing[1]= vector_Cartesian[0]*cosphi_0 
                          -vector_Cartesian[1]*sinphi_0 ;
                 
  local_north_pointing[2]= vector_Cartesian[0]*costheta_0*sinphi_0 
                          +vector_Cartesian[1]*costheta_0*cosphi_0
                          +vector_Cartesian[2]*sintheta_0;
  /*Rotating from geographical North'local magnetic' North to 'local magnetic' North*/
  vector_AIRES[2]=   local_north_pointing[2];
  vector_AIRES[1]=   local_north_pointing[0]*sindec_0 + local_north_pointing[1]*cosdec_0;
  vector_AIRES[0]=   local_north_pointing[0]*cosdec_0 - local_north_pointing[1]*sindec_0;
}


   /* chi *********************************************************** 
   Evaluation of parameter \chi=(E/m)*(B/B_{cr}).
   As a byproduct, it records the altitude in this_particle->altitude */
double chi(struct particle *this_particle)
{
  double r[3];     /*Cartesian (geocentric) coordinates for this_particle*/
  double k[3];     /*Cartesian direction of motion of this_particle*/
  double bx,by,bz; /*Cartesian geomagnetic field at location*/
  double latitude,longitude,altitude; /*Geographic coords*/
  double sintheta,costheta,sinphi,cosphi; /*lat=theta, long=phi*/
  double bDOTk,bperp,auxchi;
  double vsquared,vDOTgeog_site,equis;
  /*Rotating possition to cartesian coords*/  
  rotAIRES2Cartesian(&(this_particle->position[0]),&r[0]);
  vsquared      = r[0]*r[0]  +  r[1]*r[1]  +  r[2]*r[2];
  vDOTgeog_site = r[0]*x_0[0] + r[1]*x_0[1] + r[2]*x_0[2];
  equis=(vsquared+2*vDOTgeog_site)/sqr(earth_radius);

  /*Shifting possition (changing origin to Earth's center))*/
  r[0]+=x_0[0];
  r[1]+=x_0[1];
  r[2]+=x_0[2];

  /*Rotating direction of motion to cartesian coords*/  
  rotAIRES2Cartesian(&(this_particle->motion_direction[0]),&k[0]);

  /*Getting geographical coords of this_particle*/
  spheric_coords(r[0],r[1],r[2],
                 &altitude, &latitude, &longitude,
                 &sintheta,&costheta,&sinphi,&cosphi);
  altitude=earth_radius*heighfactor(equis);

  cartesian_magnetic_field(latitude, longitude, altitude,
                              sintheta, costheta, 
                              sinphi, cosphi,
                              &bx, &by, &bz);

  /*Returning magnetic field in AIRES coordinates*/
  double magnetic_Cartesian[3];
  magnetic_Cartesian[0]=bx;
  magnetic_Cartesian[1]=by;
  magnetic_Cartesian[2]=bz;
  rotCartesian2AIRES(magnetic_Cartesian,lastbfield.direction);
  lastbfield.strength=sqrt( sqr(lastbfield.direction[0])
                           +sqr(lastbfield.direction[1])
                           +sqr(lastbfield.direction[2]) );
  int j;
  for(j=0;j<3;j++){lastbfield.direction[j]/=lastbfield.strength;}
  lastbfield.strength/=bcrit;

  /*Returning particle's altitude*/
  this_particle->altitude=altitude;

  /*getting transversal (to motion) magnetic field*/
  bDOTk =  k[0]*bx + k[1]*by + k[2]*bz;
  bperp = sqrt(sqr(bx-k[0]*bDOTk)
	     + sqr(by-k[1]*bDOTk)
	     + sqr(bz-k[2]*bDOTk) );

  /*Returning \chi*/
  auxchi=(this_particle->energy*e/(2*electron_mass*c*c))*(bperp/bcrit);
  return(auxchi);
}



   /* frenchi ********************************************************
   Evaluation of parameter \chi=(E/m)*(B/B_{cr})
   under the approx. of dipolar geomagnetic field.
   As a byproduct, it records the altitude in this_particle->altitude */
double frenchi(struct particle *this_particle)
{
  double r[3];     /*Cartesian (geocentric) coordinates for this_particle*/
  double k[3];     /*Cartesian direction of motion of this_particle*/
  double bx,by,bz; /*Cartesian geomagnetic field at location*/
  double latitude,longitude,altitude; /*Geographic coords*/
  double sintheta,costheta,sinphi,cosphi; /*lat=theta, long=phi*/
  double bDOTk,bperp,auxchi;
  double vsquared,vDOTgeog_site,equis;
  void   geomagneticr_();
  /*Rotating possition to cartesian coords*/  
  rotAIRES2Cartesian(&(this_particle->position[0]),&r[0]);
  vsquared      = r[0]*r[0]  +  r[1]*r[1]  +  r[2]*r[2];
  vDOTgeog_site = r[0]*x_0[0] + r[1]*x_0[1] + r[2]*x_0[2];
  equis=(vsquared+2*vDOTgeog_site)/sqr(earth_radius);

  /*Shifting possition (changing origin to Earth's center))*/
  r[0]+=x_0[0];
  r[1]+=x_0[1];
  r[2]+=x_0[2];

  /*Rotating direction of motion to cartesian coords*/  
  rotAIRES2Cartesian(&(this_particle->motion_direction[0]),&k[0]);

  /*Getting geographical coords of this_particle*/
  spheric_coords(r[0],r[1],r[2],
                 &altitude, &latitude, &longitude,
                 &sintheta,&costheta,&sinphi,&cosphi);
  altitude=earth_radius*heighfactor(equis);

  cartesian_magnetic_field(latitude, longitude, altitude,
                              sintheta, costheta, 
                              sinphi, cosphi,
                              &bx, &by, &bz);

  /*Returning magnetic field in AIRES coordinates*/
  double magnetic_Cartesian[3];
  magnetic_Cartesian[0]=bx;
  magnetic_Cartesian[1]=by;
  magnetic_Cartesian[2]=bz;
  rotCartesian2AIRES(magnetic_Cartesian,lastbfield.direction);
  lastbfield.strength=sqrt( sqr(lastbfield.direction[0])
                           +sqr(lastbfield.direction[1])
                           +sqr(lastbfield.direction[2]) );
  int j;
  for(j=0;j<3;j++){lastbfield.direction[j]/=lastbfield.strength;}
  lastbfield.strength/=bcrit;

  /*Returning particle's altitude*/
  this_particle->altitude=altitude;

  /*getting transversal (to motion) magnetic field*/
  bDOTk =  k[0]*bx + k[1]*by + k[2]*bz;
  bperp = sqrt(sqr(bx-k[0]*bDOTk)
	     + sqr(by-k[1]*bDOTk)
	     + sqr(bz-k[2]*bDOTk) );

  /*Returning \chi*/
  auxchi=(this_particle->energy*e/(2*electron_mass*c*c))*(bperp/bcrit);
  return(auxchi);
}
