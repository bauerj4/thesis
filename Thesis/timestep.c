/*
  Segment of code that handles forcing particles into same timestep bin.

  diskDVars is a differencec array that will be defined later. All other
  variables are explained in allvars.h/c.
*/


#elif defined(RIGID_PARTICLE_DISK)
double rigidDiskAcc;
      rigidDiskAcc = diskDVars[3] * diskDVars[3] + diskDVars[4] *\
	diskDVars[4] +						 \
	diskDVars[5] * diskDVars[5];
rigidDiskAcc = pow(rigidDiskAcc, 0.5);
if(P[p].Type == 2){
  dt = sqrt(2 * All.ErrTolIntAccuracy *\
	    All.cf_atime * All.SofteningTable[P[p].Type] / \
	    rigidDiskAcc)/TIMESTEP_REDUCTION;
 }

 else
   dt = sqrt(2 * All.ErrTolIntAccuracy * All.cf_atime *\
	     All.SofteningTable[P[p].Type] / ac);

#else
