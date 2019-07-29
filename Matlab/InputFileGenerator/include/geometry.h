#ifndef GEOMETRY_H_INCLUDED
#define GEOMETRY_H_INCLUDED

void
GEOMETRY_initialize();

void GEOMETRY_particle_type(const double &x, const double &y, const double &z);

void GEOMETRY_particle_type_column(const double &x, const double &y, const double &z);

void GEOMETRY_particle_type_duct(const double &x, const double &y, const double &z);

void GEOMETRY_particle_type_vortex(const double &x, const double &y, const double &z);

void GEOMETRY_particle_type_3D_Bucket(const double &x, const double &y, const double &z);

void GEOMETRY_particle_type_3D_Simple(const double &x, const double &y, const double &z);

void GEOMETRY_particle_type_3D_spinningBall(const double &x, const double &y, const double &z);

void GEOMETRY_particle_type_3D_CircularPoiseuilleFlow(const double &x, const double &y, const double &z);

void GEOMETRY_particle_type_3D_StirringTank(const double &x, const double &y, const double &z);

#endif // GEOMETRY_H_INCLUDED
