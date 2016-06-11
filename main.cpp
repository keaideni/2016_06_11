#include "Parameter.h"
#include "OP.h"
#include "SuperEnergy.h"
#include "DMRGP.h"



int OP::Max;
int main()
{
	Parameter para;
	para.read();
	
	OP::Max = para.ParticleNo + 1;

	
	/*Sub a(para, 1), b(para, 2);

	Super sup(para, a, b, a, b, 1);
	SuperEnergy supp(para, sup);

	std::cout << para.Energy << std::endl;

	sup.Wave.show();*/
	
	







	DMRGP DMRG(para);

	
	system("pause");
}