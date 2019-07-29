#ifndef STRUCT_H_INCLUDED
#define STRUCT_H_INCLUDED

#include <vector>

using namespace std;

struct structParameter{
	char nameOfInputFile[100];

	char filename_data[100];	         // filename of datafile to output
	char filename_grid[100];	         // filename of gridfile to output
	char filename_concentration[100];    // filename of concentrationfile to output
	char filename_forcedMotion[100];     // filename of forcedMotion file to output
	char filename_bubble[100];           // filename of forcedMotion file to output
	// simulation parameters
	double dt;                           // length of time step in simulation [s]
	double finishTime;                   // length of whole simulation [s]
	double l_0;		                     // average distance between particles
	double intf;		                 // depth of fluid in tank
	double omega;		                 // stirring speed [deg/s]
	double c_0;                          // initial impurities [mol/m^3]
	// physical properties
	double liquid_density;               // density of the liquid [kg/m^3]
	double gas_density;                  // density of the gas [kg/m^3]
	double liquid_viscosity;             // viscosity of the liquid [Pa s]
	double temperature;                  // temperature of the tank [K]
	double pressure;                     // head pressure [Pa]
	double diffusion_coeff;              // diffusion_coeff for diffusion of impurities [m^2/s]
	double surface_tension;              // surface tension of liquid [N/m]
	double henry_coeff;                  // henry coeff of gas [Pa m^3/mol]
	double molecular_weight;             // the molecular weight of gas [g/mol]
	double theta;                        // sesshokukaku of gas bubble against wall [.]
	// fitting parameter?
	int number_of_Particle_Types;        // number of different Particle Types
	int number_of_Liquid_Types;          // number of different Liquid Types
	double Diffusion_Coefficient;        // Diffusion Coeffiecient (D)
	double max_Velocity;                 // Max Velocity of Poiseuille Flow

	int wallParticleType;                 
	int dummyWallParticleType;
	int forcedMotionParticleType;

	double eps;
	double width;
	double length;
	int max_num_particles;

	int type;
	double velocity[3];

	char flagOfDomainAutoSetting[100];
	double domainUpperLimit[3];
	double domainLowerLimit[3];

	char flagOfPeriodic[100];
	char flagOfInterfaceCorrection[100];
	char flagOfSingleVortex[100];
	char flagOfForcedMotion[100];

	double gravity;
};

struct structGeometry{

	double minCor[3];
	double maxCor[3];
	double minBucket[3];
	double maxBucket[3];
	double center[3];
	double width[3];
	double diameter;
	
	int dimensions;
	int geometryType;
};

extern struct structParameter     parameter;
extern struct structGeometry      geometry;

#endif // STRUCT_H_INCLUDED
