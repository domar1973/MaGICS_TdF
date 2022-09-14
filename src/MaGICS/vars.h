
/* Global internal vars */

int    debuglvl;  /*if==0 => save debug data*/
int    AllConverted; /*if==1 => all photons are converted*/
int    ThereWasNoConversion=1,FirstInteractionAlready=0;
int    PrimaryParticle=1;

double bcrit;      /* Shwinger magnetic field (set at inic.c) */
double geographical_location[2];     /* [latitude theta, longitude phi] */
double sintheta_0;
double costheta_0;
double sinphi_0;
double cosphi_0;           
double site_declination,sindec_0,cosdec_0;
         /* the _0 indicate values for geographical location */
double x_0[3];
         /* Cartesian coordinates of geographical site */
double initial_position[3];
double injection_altitude;
double fdate;
double mytime = 0;
double primary_energy;                   /* Incident Photon */
int total_photon_number=0;
int hard_photons=0;
double atm_injection_altitude;  /* Altitude of atmospheric injection */
double altitude_0;              /* Altitude of Primary injection */
char *args[20];
