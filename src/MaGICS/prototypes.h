void (*propagate_particle)(); 
/*Pointer to propagator
  (1D_erberaprox, 1D_klepikov, etc */

double    urandom_(), urandomt_(), zfromdepth_();    /*__*/

void      speistart_();

void      speiend_();
void      speigetparsc();
void speitaskc();
void sprimnamec();
void getglobalc();
void speimv_();
void spaddp0_();
void spinjpoint_();
void croinputdata0_();
void speisetrealvar_();
void sprecordvars_();



void   geomagnetic_();

double dbskr3_();     /*BesselK[nu/3,x]=BSKR3_(&x,&nu), nu=-2,-1,1,2*/
double sqr();
