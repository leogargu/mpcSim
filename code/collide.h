#ifndef COLLIDE_H /* header guard */
#define COLLIDE_H
//#pragma once

#include "mpc.h"
#include <gsl/gsl_rng.h>

inline int cell_coord2cell_idx(int * cell_coord, const int * n_cells_dim);
inline int pos2cell_idx(const Geometry geometry, const double * pos, const double * shift );

void collide(const int n_part, const double T, const double m, const double m_inv, Geometry geometry, double * shift, gsl_rng * r, 
	int * c_p, double * cell_mass ,double ** cell_vel, double ** cell_rnd_vel, double ** pos, double ** vel);


#endif /* header guard */
