#include <Eigen/Core>
#include <SymEigsSolver.h>
#include <iostream>
#include "Super.h"

using namespace Spectra;

#ifndef SUPERENERGY_H
#define SUPERENERGY_H
class SuperEnergy
{
public:
	QWave wave;


	SuperEnergy(Parameter&para,Super& sup)
	{
		wave = sup.Wave;
		SymEigsSolver<double, SMALLEST_ALGE, Super> eigs(&sup, 1, 6);
		eigs.init();
		eigs.compute();
		if (eigs.info() == SUCCESSFUL)
		{
			wave.f2Wave(eigs.eigenvectors(1));
			para.Energy = eigs.eigenvalues()(0);
			//std::cout << eigs.num_iterations() << std::endl;
		}

		
	};
	

};







#endif
