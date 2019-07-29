#include <iostream>
#include <fstream>
#include <math.h>
#include <ctime>
#include <string>
#include <sstream>
#include <cmath>
#include <iomanip>


#include "struct.h"
#include "file.h"
#include "geometry.h"

using namespace std;


void
FILE_readDataFile(void){
    ifstream inputfile;
    ifstream ifile(parameter.nameOfInputFile);
    if (!ifile) {
        cout << "Input File with the name " << parameter.nameOfInputFile << " does not exist." << endl;
        exit(1);
    }

    inputfile.open (parameter.nameOfInputFile);
    char tmp[100];

    if (inputfile.is_open())
  {
      inputfile >> tmp;
      inputfile >> tmp;
      inputfile >> tmp;
      inputfile >> parameter.filename_data;
      inputfile >> tmp;
      inputfile >> parameter.filename_grid;
      inputfile >> tmp;
	  inputfile >> parameter.filename_forcedMotion;
	  inputfile >> tmp;
	  inputfile >> parameter.filename_concentration;
	  inputfile >> tmp;
	  inputfile >> tmp;
	  inputfile >> parameter.dt;
	  inputfile >> tmp;
	  inputfile >> parameter.finishTime;
	  inputfile >> tmp;
	  inputfile >> parameter.l_0;
	  inputfile >> tmp;
	  inputfile >> parameter.intf;
	  inputfile >> tmp;
	  inputfile >> parameter.omega;
	  inputfile >> tmp;
	  inputfile >> parameter.c_0;
	  inputfile >> tmp;
	  inputfile >> tmp;
	  inputfile >> parameter.liquid_density;
	  inputfile >> tmp;
	  inputfile >> parameter.gas_density;
	  inputfile >> tmp;
	  inputfile >> parameter.liquid_viscosity;
	  inputfile >> tmp;
	  inputfile >> parameter.temperature;
	  inputfile >> tmp;
	  inputfile >> parameter.pressure;
	  inputfile >> tmp;
	  inputfile >> parameter.diffusion_coeff;
	  inputfile >> tmp;
	  inputfile >> parameter.surface_tension;
	  inputfile >> tmp;
	  inputfile >> parameter.henry_coeff;
	  inputfile >> tmp;
	  inputfile >> parameter.molecular_weight;
	  inputfile >> tmp;
	  inputfile >> parameter.theta;
	  inputfile >> tmp;
	  inputfile >> parameter.number_of_Particle_Types;
	  parameter.number_of_Particle_Types +=2;
	  inputfile >> tmp;
	  inputfile >> parameter.Diffusion_Coefficient;
	  inputfile >> tmp;
	  inputfile >> parameter.max_Velocity;
	  inputfile >> tmp;
	  inputfile >> tmp;
	  inputfile >> geometry.geometryType;

      inputfile.close();
  }
  else cout << "Unable to open file";
}


void
FILE_showInput(void){

cout << "creating grid file for stirred tank...\n"
			<< "--------------------------------------------------------\n"
			<< "(output filenames)\n"
			<< "    data filename              :  " << parameter.filename_data << "\n"
			<< "    grid filename              :  " << parameter.filename_grid << "\n"
			<< "    concentration filename     :  " << parameter.filename_concentration << "\n"
			<< "    forcedMotion filename      :  " << parameter.filename_forcedMotion << "\n"
			<< "(simulation parameters)\n"
			<< "    dt [s]                     :  " << parameter.dt << "\n"
			<< "    finishTime [s]             :  " << parameter.finishTime << "\n"
			<< "    l_0 [m]                    :  " << parameter.l_0 << "\n"
			<< "    intf [m]                   :  " << parameter.intf << "\n"
			<< "    omega [deg/s]              :  " << parameter.omega << "\n"
			<< "    c_0 [mol/m^3]              :  " << parameter.c_0 << "\n"
			<< "(physical properties)\n"
			<< "    liquid_density [kg/m^3]    :  " << parameter.liquid_density << "\n"
			<< "    gas_density [kg/m^3]       :  " << parameter.gas_density << "\n"
			<< "    liquid_viscosity [Pa s]    :  " << parameter.liquid_viscosity << "\n"
			<< "    temperature [K]            :  " << parameter.temperature << "\n"
			<< "    pressure [Pa]              :  " << parameter.pressure << "\n"
			<< "    diffusion_coeff [m^2/s]    :  " << parameter.diffusion_coeff << "\n"
			<< "    surface_tension [N/m]      :  " << parameter.surface_tension << "\n"
			<< "    henry_coeff [Pam^3/mol]    :  " << parameter.henry_coeff << "\n"
			<< "    molecular_weight [g/mol]   :  " << parameter.molecular_weight << "\n"
			<< "    theta [.]                  :  " << parameter.theta << "\n"
			<< "(fitting parameter?)\n"
			<< "    number_of_Particle_Types   :  " << parameter.number_of_Particle_Types << "\n"
			<< "    Diffusion_Coefficient      :  " << parameter.Diffusion_Coefficient << "\n"
			<< "(geometry)\n"
			<< "    typeOfGeometry             :  " << geometry.geometryType << "\n"
			<< "--------------------------------------------------------\n"
			<< endl;

}




void
FILE_output_data_file(){
	int iType;
	ofstream ofs_data(parameter.filename_data);

	ofs_data << "#######--FREQUENTLY-USED-DATA-#####################################" << "\n"
		<< "#--------ParticleSize----------------------------------------------" << "\n"
		<< "averageDistanceBetweenParticles(m)                 " << parameter.l_0 << "\n"
		<< "#--------Time------------------------------------------------------" << "\n"
		<< "finishTime(sec)                                    " << parameter.finishTime << "\n"
		<< "inititalDt(sec)                                    " << parameter.dt << "\n"
		<< "#--------File------------------------------------------------------" << "\n"
		<< "TypeOfOutputInterval(simulationTime/timeStep)      simulationTime" << "\n"
		<< "--simulationTime---interval(sec)                   1.0e-1" << "\n"
		<< "--timeStep---------interval(timeStep)              1" << "\n"
		<< "#--------Domain----------------------------------------------------" << "\n"
		<< "numberOfDimensions(2or3)                           " << geometry.dimensions << "\n"
		<< "autoSettingOfDomainSize(on/off)                    " << parameter.flagOfDomainAutoSetting << "\n"
		<< "--on---upperMarginRatio[XDIM](ratio)               0.0" << "\n"
		<< "--on---upperMarginRatio[YDIM](ratio)               0.0" << "\n"
		<< "--on---upperMarginRatio[ZDIM](ratio)               0.0" << "\n"
		<< "--on---lowerMarginRatio[XDIM](ratio)               0.0" << "\n"
		<< "--on---lowerMarginRatio[YDIM](ratio)               0.0" << "\n"
		<< "--on---lowerMarginRatio[ZDIM](ratio)               0.0" << "\n"
		//<< "--off--upperLimit[XDIM](m)                         2.01" << "\n"
		<< "--off--upperLimit[XDIM](m)                         " << parameter.domainUpperLimit[0] << "\n"
		<< "--off--upperLimit[YDIM](m)                         " << parameter.domainUpperLimit[1] << "\n"
		<< "--off--upperLimit[ZDIM](m)                         " << parameter.domainUpperLimit[2] << "\n"
		<< "--off--lowerLimit[XDIM](m)                         " << parameter.domainLowerLimit[0] << "\n"
		<< "--off--lowerLimit[YDIM](m)                         " << parameter.domainLowerLimit[1] << "\n"
		<< "--off--lowerLimit[ZDIM](m)                         " << parameter.domainLowerLimit[2] << "\n"
		<< "#--------RadiusOfInteractionDomain--------------------------------" << "\n"
		<< "radiusOfParticleNumberDensity(ratio)               2.1" << "\n"
		<< "radiusOfGradient(ratio)                            2.1" << "\n"
		<< "radiusOfLaplacianForViscosity(ratio)               3.1" << "\n"
		<< "radiusOfLaplacianForPressure(ratio)                3.1" << "\n"
		<< "radiusOfLaplacianForDiffusion(ratio)               3.1" << "\n"
		<< "########                                           ***" << "\n"
		<< "#######--LESS-FREQUENTLY-USED-DATA-###############################" << "\n"
		<< "#--------Bucket(subDomain)----------------------------------------" << "\n"
		<< "autoSettingOfBucketCapacity                        on" << "\n"
		<< "--on---optimizationOfMemorySize(on/off)            on" << "\n"
		<< "--on-----off--marginRatio(ratio)                   2.0" << "\n"
		<< "--off--capacityOfBucket                            40" << "\n"
		<< "#--------TableOfNeighborParticles---------------------------------" << "\n"
		<< "theNumberOfNeighborTables(1or2)                    2" << "\n"
		<< "autoSettingOfCapacityOfNeighborTable               on" << "\n"
		<< "--on---marginRatioOfLargeTable(ratio)              2.5" << "\n"
		<< "--on---marginRatioOfSmallTable(ratio)              2.5" << "\n"
		<< "--off--capacityOfLargeNeighborTable                100" << "\n"
		<< "--off--capacityOfSmallNeighborTable                100" << "\n"
		<< "#--------ProfFile-------------------------------------------------" << "\n"
		<< "divisionOfProfFile(on/off)                         on" << "\n"
		<< "--on---compressionOfProfFile                       off" << "\n"
		<< "--on---upperLimitOfNumberOfFiles                   3001" << "\n"
		<< "#--------PressureFile---------------------------------------------" << "\n"
		<< "outputPressureFile(on/off)                         off" << "\n"
		<< "--on---outputPressureOfAllWallParticle(on/off)     off" << "\n"
		<< "---------off--nameOfInputFile                      input.press" << "\n"
		<< "--on---divisionOfPressureFile(on/off)              off" << "\n"
		<< "---------on---upperLimitOfNumberOfFiles            30000" << "\n"
		<< "---------on---nameOfOutputFile                     output_" << "\n"
		<< "---------off--nameOfOutputFile                     output.press" << "\n"
		<< "#--------TorqueFile---------------------------------------------" << "\n"
		<< "outputTorqueFile(on/off)                           off" << "\n"
		<< "--on--nameOfOutputFile                             output.torque" << "\n"
		<< "numberOfRotations(/s)                              " << parameter.omega / 360.0 << "\n"
		<< "#--------TypeOfParticle-------------------------------------------" << "\n"
		<< "numberOfParticleTypes                              " << parameter.number_of_Particle_Types << "\n"
		<< "--typeNumberOfWallParticle                         " << parameter.wallParticleType << "\n"
		<< "--typeNumberOfDummyWallParticle                    " << parameter.dummyWallParticleType << "\n"
		<< "#--------MassDensity----------------------------------------------" << endl;
	for (iType = 0; iType <= parameter.number_of_Particle_Types - 1; iType++){
		ofs_data << "massDensityOfType" << iType << "(kg/m^3)                         " << parameter.liquid_density << endl;
	}
	ofs_data << "#--------Compressibility------------------------------------------" << endl;
	for (iType = 0; iType <= parameter.number_of_Particle_Types - 1; iType++){
		ofs_data << "compressibilityOfType" << iType << "                             1.0e-7" << endl;
	}
	ofs_data << "#--------Viscosity------------------------------------------------" << "\n"
		<< "flagOfViscosityCalculation(on/off)                 on" << "\n"
		<< "--on--flagOfHighViscosityCalculation(on/off)       on" << "\n"
		<< "--on--kinematicViscosity(m^2/s)                    "
		<< parameter.liquid_viscosity / parameter.liquid_density << "\n"
		<< "#--------Bubble---------------------------------------------------" << "\n"
		<< "flagOfBubbleCalculationy(on/off)                   off" << "\n"
		<< "--on--nameOfInputfile                              input.bubble" << "\n"
		<< "--on--numperOfParticleForCalculatingBeta           1344" << "\n"
		<< "--on--massDensityOfBubble(kg/m^3)                  " << parameter.gas_density << "\n"
		<< "--on--gassConstant(J/K*mol)                        8.3144621" << "\n"
		<< "--on--temperature(K)                               " << parameter.temperature << "\n"
		<< "--on--headPressure(Pa)                             " << parameter.pressure << "\n"
		<< "#-----------------------------------------------------------------" << "\n"
		<< "########                                           ***" << "\n"
		<< "########                                           ***" << "\n"
		<< "########                                           ***" << "\n"
		<< "########                                           ***" << "\n"
		<< "#--------Gravity--------------------------------------------------" << "\n"
		<< "gravity[XDIM](m/s^2)                               0.0" << "\n"
		<< "gravity[YDIM](m/s^2)                               " << parameter.gravity << "\n"
		<< "gravity[ZDIM](m/s^2)                               0.0" << "\n"
		<< "#--------TimeDifference-------------------------------------------" << "\n"
		<< "courantNumber                                      0.2" << "\n"
		<< "diffusionNumber                                    0.2" << "\n"
		<< "maxDt(ratio)                                       1.0" << "\n"
		<< "minDt(ratio)                                       0.001" << "\n"
		<< "upperLimitOfChangeRateOfDt(ratio)                  1.2" << "\n"
		<< "#--------SolvelerOfSimultaneousEquations--------------------------" << "\n"
		<< "upperLimitOfIterationNumber                        2000" << "\n"
		<< "lowerLimitOfIterationNumber                        10" << "\n"
		<< "smallNumberForCheckingConvergence                  1.0e-9" << "\n"
		<< "#--------CollisionBetweenParticles--------------------------------" << "\n"
		<< "collisionDistance(ratio)                           0.5" << "\n"
		<< "collisionCoefficient                               0.2" << "\n"
		<< "#--------DirichletBoundaryCondition-------------------------------" << "\n"
		<< "thresholdOfParticleNumberDensity(ratio)            0.97" << "\n"
		<< "#--------Other----------------------------------------------------" << "\n"
		<< "relaxationCoefficientOfLambda                      0.2" << "\n"
		<< "finishTimeStep(timeStep)                           1000000" << "\n"
		<< "######---CustomizedFunctions---###################################" << "\n"
		<< "#--------RigidBody------------------------------------------------" << "\n"
		<< "forcedMotionOfRigidBody(on/off)                    " << parameter.flagOfForcedMotion << "\n"
		<< "--on---typeNumberOfRigidParticle                   " << parameter.forcedMotionParticleType << "\n"
		<< "--on---fileNameOfSamplingDataForRigidBody          input.forcedMotion" << "\n"
		<< "#--------TANAKAandMASUNAGAModel-----------------------------------" << "\n"
		<< "flagOfTanakaAndMasunagaModel(on/off)               on" << "\n"
		<< "--on--valueOfGamma                                 0.01" << "\n"
		<< "--on--valueOfC                                     1.0" << "\n"
		<< "########                                           ***" << "\n"
		<< "########                                           ***" << "\n"
		<< "#--------ConcentrationCalculation---------------------------------" << "\n"
		<< "flagOfConcentrationCalculation(on/off)             on" << "\n"
		<< "DiffusionCoefficient                               " << parameter.Diffusion_Coefficient << "\n"
		<< "--on---nameofInputfile                             input.concentration" << "\n"
		<< "#---------InterfaceCalculation------------------------------------" << "\n"
		<< "flagOfInterfaceCalculation(on/off)                 on" << "\n"
		<< "flagOfInterfaceCorrection(on/off)                  " << parameter.flagOfInterfaceCorrection << "\n"
		<< "flagOfInterfaceCompare(on/off)                     on" << "\n"
		<< "flagOfPeriodicX(on/off)                            " << parameter.flagOfPeriodic << "\n"
		<< "#---------ParticleRedistribution----------------------------------" << "\n"
		<< "flagOfSingleVortex(on/off)                         " << parameter.flagOfSingleVortex << "\n"
		<< "alphaZero                                          0.0001" << "\n"
		<< "########                                           ***" << "\n"
		<< "########                                           ***" << "\n"
		<< "########                                           ***" << "\n"
		<< "########                                           ***" << "\n"
		<< "########                                           ***" << "\n"
		<< "########                                           ***" << "\n"
		<< "########                                           ***" << "\n"
		<< "#-----------------------------------------------------------------" << "\n"
		<< "########                                           ***" << "\n"
		<< "########                                           ***" << "\n"
		<< "########                                           ***" << "\n"
		<< "#-----------------------------------------------------------------" << "\n"
		<< "########                                           ***" << "\n"
		<< "#-----------------------------------------------------------------" << "\n"
		<< "########                                           ***" << "\n"
		<< "########                                           ***" << "\n"
		<< "########                                           ***" << "\n"
		<< "########                                           ***" << "\n"
		<< "########                                           ***" << "\n"
		<< "########                                           ***" << "\n"
		<< "########                                           ***" << endl;
}


void
FILE_output_grid_bubble_file(){
	int Nu = parameter.number_of_Particle_Types;
	int cnt = 0;               // number of particles added so far
	int cnt_fluid = 0;         // number of fluid (0) particles added so far

	//--------------
	// make grid
	//--------------
	stringstream ss;
	
	// input origin first (for numperOfParticleForCalculatingBeta in input.data)
	GEOMETRY_particle_type(0.0, 0.0, 0.0);
	if (parameter.type >= 0) {
		++cnt;  // count all particles
		if (parameter.type >= 0 && parameter.type < 2 - parameter.eps) {
			++cnt_fluid;  // count fluid particles
		}
		ss << parameter.type << " " << 0.0 << " " << 0.0 << " " << 0.0 << " "
			<< parameter.velocity[0] << " " << parameter.velocity[1] << " " << parameter.velocity[2] << " " << 0.0 << " " << 0.0 << "\n";
	}
	for (double x = geometry.minCor[0]; x <= geometry.maxCor[0] + parameter.eps; x += parameter.l_0) {
		for (double y = geometry.minCor[1]; y <= geometry.maxCor[1] + parameter.eps; y += parameter.l_0) {
			for (double z = geometry.minCor[2]; z <= geometry.maxCor[2] + parameter.eps; z += parameter.l_0) {
				if (fabs(x)<0.000001) x=0;
				if (fabs(y)<0.000001) y=0;
				if (fabs(z)<0.000001) z=0;
				if (x == 0 && y == 0 && z == 0)
					continue;
				GEOMETRY_particle_type(x, y, z);
				if (parameter.type >= 0) {
					++cnt;  // count all particles
					if (parameter.type >= 0 && parameter.type < parameter.number_of_Liquid_Types - parameter.eps) {
						++cnt_fluid;  // count fluid particles
					}
					ss << parameter.type << " " << x << " " << y << " " << z << " "
						<< parameter.velocity[0] << " " << parameter.velocity[1] << " " << parameter.velocity[2] << " " << 0.0 << " " << 0.0 << "\n";
				}
			}
		}
	}
cout << cnt << endl << cnt_fluid << endl << cnt-cnt_fluid << endl;

	//---------------------------------
	// output .grid and .bubble file
	//---------------------------------
	ofstream ofs_grid(parameter.filename_grid);
	ofstream ofs_concentration(parameter.filename_concentration);
	// TODO check
	cout << "number of particles :             " << cnt << "\n"
		<< "number of fluid particles :       " << cnt_fluid << endl;
	// ofs_grid << showpoint;
	ofs_grid << "0.000000\n" << cnt << "\n";
	for (int i = 0; i < cnt; ++i) {
		int type;
		double x, y, z, vx, vy, vz, p, n;
		ss >> type >> x >> y >> z >> vx >> vy >> vz >> p >> n;
		ofs_grid << type << " " << x << " " << y << " " << z << " "
			<< vx << " " << vy << " " << vz << " " << p << " " << n << "\n";

		if (type == 0)
			ofs_concentration << 1 << " " << 0 << " " << 0 << "\n";
		else if (type == 1 && (Nu == 4 || Nu == 5))
			ofs_concentration << 0 << " " << 1 << " " << 0 << "\n";
		else if (type == 2 && Nu == 6)
			ofs_concentration << 0 << " " << 0 << " " << 1 << "\n";
		else
			ofs_concentration << -1 << " " << -1 << " " << -1 << "\n";

		//if (type == Nu - 1)
		//	ofs_concentration << 0 << " " << 0 << " " << 0 << "\n";
		//else if (y <= 0.3 * H && Nu == 6)
		//	ofs_concentration << 0 << " " << 0 << " " << 1 << "\n";
		//else if (y <= 0.7 * H && Nu == 6)
		//	ofs_concentration << 0 << " " << 1 << " " << 0 << "\n";
		//else if (y <= 0.3 * H && Nu == 5)
		//	ofs_concentration << 0 << " " << 1 << " " << 0 << "\n";
		//else /*if (y <= H)*/
		//	ofs_concentration << 1 << " " << 0 << " " << 0 << "\n";
	}
}


void
FILE_output_forcedMotion_file(){
	ofstream ofs_forced(parameter.filename_forcedMotion);
	// the header
	ofs_forced << "#--------CenterOfRotation-----------------------------------------------------" << "\n"
		<< "position(centerOfGravity/centerOfWidth/directInput)             centerOfWidth" << "\n"
		<< "---directInput---centerOfRotation[XDIM](m)                      0.000000" << "\n"
		<< "---directInput---centerOfRotation[YDIM](m)                      0.000000" << "\n"
		<< "---directInput---centerOfRotation[ZDIM](m)                      0.000000" << "\n"
		<< "#--------SamplingData---------------------------------------------------------" << "\n"
		<< "numberOfSamplingData                                            " << (int)parameter.finishTime / parameter.dt << "\n"
		<< "startingTime(sec)                                            0.000000" << "\n"
		<< "SamplingTime(sec)-----X(m)-----Y(m)-----Z(m)-----thetaX(degree)---thetaY(degree)---thetaZ(degree)" << "\n";

	// data part
	// ofs_forced << "0 0 0 0 0 0 0\n";
	// double prev_theta_y = 0;
	// for (int k = 1; k < finishTime/dt; ++k) {
	//     const double diff_theta_y = min(k * omega * dt * dt, omega * dt);
	//     ofs_forced << dt * k << " 0 0 0 0 " << prev_theta_y + diff_theta_y << " 0\n";
	//     prev_theta_y += diff_theta_y;
	// }

	ofs_forced << fixed << setprecision(5) << "0 0 0 0 0 0 0\n";
	for (int k = 1; k < parameter.finishTime / parameter.dt; ++k) {
		const double t = parameter.dt * k;
		const double theta_y = (t < 1) ? (0.5*parameter.omega*t*t) : (parameter.omega*(t - 1) + 0.5*parameter.omega);
		ofs_forced << t << " 0 0 0 0 " << theta_y << " 0\n";
	}
}
