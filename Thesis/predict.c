/*
  This section of code turns off the 
  kick for rigid disk particles. The kick 
  is handled by separate code.
*/


#ifndef SIDM_FREEZE
#ifdef EPSILON_MASS
if (P[i].Type != 1){
#endif // EPSILON_MASS
  #ifdef RIGID_PARTICLE_DISK
  if(P[i].Type != 2){
#endif // RIGID_PARTICLE_DISK
    for(j = 0; j < 3; j++)
      P[i].Pos[j] += P[i].Vel[j] * dt_drift;
    #ifdef RIGID_PARTICLE_DISK
  }
  #endif
  #ifdef EPSILON_MASS
 }

#endif // EPSILON_MASS
#endif
