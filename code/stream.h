#ifndef STREAM_H  /* header guard */
#define STREAM_H

#include "mpc.h"


inline void canonize(double Lx, double * pos);
inline double find_impact_time(double coeff[5], double dt);
inline void compute_acc(double g, double * acc);

inline int TEST_particle_in_lumen(Geometry cylinder, double * pos);
inline int TEST_all_particles_in_lumen(int n_part, Geometry cylinder, double ** pos);


void stream(double dt, const int n_part, const double g, const Geometry cylinder, double ** pos, double ** vel, double ** acc);

#endif	/* header guard */
