struct particle
{
  int    particleID;  /* The same code as AIRES */
  double energy;
  double motion_direction[3];
  double position[3];
  double altitude;
  double last_time; /*Last time of tracking*/
  struct particle *next;
};

struct bfield
{
  double strength;     /* b/B_c */
  double direction[3]; /* in AIRES coordinates */
};

struct particle *first_particle;

struct bfield lastbfield;

/* Variables taken from AIRES */
       int       shower_number;
       double    primary_energy;
       double    default_injection_position[3];
       double    injection_depth, ground_depth;
       double    ground_altitude, d_ground_inj;
       double    shower_axis[3];

       char      parstring[200], auxstring[200], task[120], spname[17];
       double realdata[99], paquesirve[99];
       int intdata[99], shprimcode[99];
       int       parstringlen, auxlen, tasklen, splen;
       int       task_version;

       double    th0 = 0.01;
       int       pipluscode = 11, piminuscode = -11;
       int       gammacode  = 1, protoncode = 31;

       int       csys, t0beta;
       double    t0, beta, wt;

       int       irc, newmv, oldmv;
       int       i, nrp;

/* Manages the linked list */
void add_particle(int ID, double energy, double position[3], 
                     double motion_direction[3], double creation_time, double altitude)
  {
  int i;
  struct particle *new_particle;
  struct particle *aux;
  new_particle=(struct particle *)malloc(sizeof(struct particle));
  new_particle->particleID=ID;
  new_particle->energy=energy;
  for (i=0;i<3;i++){new_particle->motion_direction[i]=motion_direction[i];};
  for (i=0;i<3;i++){new_particle->position[i]=position[i];}
  new_particle->last_time=creation_time;
  new_particle->altitude=altitude;
  new_particle->next=(NULL);
  if (first_particle==NULL){first_particle=new_particle;}  
  else {for(aux=first_particle;aux->next!=NULL;){aux=aux->next;};
  aux->next=new_particle;}
   }

void remove_particle(struct particle *this_particle)
{
  struct particle *aux;
  if (this_particle==first_particle)
    {
      aux=first_particle;
      first_particle=aux->next;
      free(aux);
      aux=NULL;
    }
  else
    {
    for(aux=first_particle;(aux->next!=this_particle)&&(aux!=NULL);)
            {aux=aux->next;}
    if (aux->next==this_particle)
      {
	aux->next=this_particle->next;
	free(this_particle);
      }
    }
  this_particle=NULL;
}
