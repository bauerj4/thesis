
/*
  Don't kick the rigid disk. 
*/

for(j = 0; j < 3; j++)	/* do the kick, only collisionless particles */
  {
    dvel[j] = P[i].GravPM[j] * dt_gravkick;
#ifndef SIDM_FREEZE
#ifdef EPSILON_MASS
            if(P[i].Type != 1){
#endif
#ifdef RIGID_PARTICLE_DISK
	      if(P[i].Type != 2){
#endif
	    P[i].Vel[j] += dvel[j];
	    P[i].dp[j] += P[i].Mass * dvel[j];
#ifdef RIGID_PARTICLE_DISK
	      }
#endif
#ifdef EPSILON_MASS
            }
#endif
#endif
  }
