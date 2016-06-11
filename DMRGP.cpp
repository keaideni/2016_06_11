#include "DMRGP.h"
#include "Sub.h"
#include "QWave.h"
#include "Super.h"


DMRGP::DMRGP()
{

}
DMRGP::~DMRGP()
{

}

DMRGP::DMRGP(Parameter& para)
{
	clock_t Abegin;
	double allT;
	int OS(1);
	int OE;
	int dir(1);


	//===========this label for the way of growth the blocks======
	OrbitalM = 4;
	OrbitalN = 4;
	Gdir = 1;
	//====================================================



	OE = para.LatticeSize;


	SaveAll.open("./result/SaveAll");
	if (!SaveAll.is_open())
	{
		std::cout << "the file doesn't exit!" << std::endl;
	}

	saveT = 0;

	SaveAll << "==============build up================" << std::endl;
	Abegin = clock();





	BuildUpP(para, OS, OE, dir);

	allT = difftime(clock(), Abegin) / CLOCKS_PER_SEC;
	SaveAll << "===========build up finished============" << std::endl;

	SaveAll << "The build up process takes " << allT << " seconds!" << std::endl;
	SaveAll << "The eigenstate calculation takes " << saveT << " seconds!" << std::endl;

	SaveAll << "===============Sweep================" << std::endl;
	saveT = 0;
	Abegin = clock();





	SweepP(para, OS, OE, dir);

	allT = difftime(clock(), Abegin) / CLOCKS_PER_SEC;
	SaveAll << "===========Sweep finished==============" << std::endl;
	SaveAll << "The build up process takes " << allT << " seconds!" << std::endl;
	SaveAll << "The eigenstate calculation takes " << saveT << " seconds!" << std::endl;
	SaveAll.close();

	Fdata.open("./result/data");
	if (!Fdata.is_open())
	{
		std::cout << "the file doesn't exit!" << std::endl;
	}


	Fdata << "Q = " << para.ParticleNo << "    LatticeSize = " << std::setw(4) << para.LatticeSize << ",      gr = " << std::setw(4) << para.gr << ",    gl = " << std::setw(4) << para.gl
		<< ",        MiuP = " << std::setprecision(15) << MiuP << ",    MiuN = " << std::setprecision(15) << MiuN << ",    trace = " << std::setprecision(15) << FTrace
		<< ",    truncerr = " << std::setprecision(15) << FTruncerr << std::endl;

	Fdata.close();
}










void DMRGP::BuildUpP(Parameter& para, int& OS, int& OE, int& dir)
{

	befortruncateP(para, OS, OE);
	//int judge = (OE-OS ==1)? 1:0;
	int i(1);

	while (true)
	{
		SaveAll << i << std::endl;

		Sys.read(OS);//Sys.show();
		Env.read(OE);//Env.show();

		m.Initial(para, OrbitalM);
		n.Initial(para, OrbitalN);

		getEnergyP(para, dir);






		if ((OE - OS) == 3)break;


		truncUpdateP(para, OS, OE, dir);

		if (Gdir == 1)
		{
			OrbitalM += 1;
		}
		else
		{
			OrbitalN += 1;
		}
		Gdir *= -1;


		i++;
	}

}


void DMRGP::befortruncateP(const Parameter& para, const int& OS, const int& OE) const
{
	Sub Sys1(para, OS);
	Sys1.save();
	//std::cout<<Sys.SubSys.Max<<std::endl;
	Sub Env1(para, 1);
	Env1.Orbital = OE;
	Env1.save();
}



void DMRGP::getEnergyP(Parameter& para, int dir)
{



	int qtot = Sys.Orbital + 1;
	//std::cout<<qtot<<std::endl;
	Super Sup(para, Sys, m, n, Env, qtot);
	//std::cout<<"hehe"<<std::endl;
	begin = clock();
	SuperEnergy Supp(para, Sup);
	saveT += difftime(clock(), begin) / CLOCKS_PER_SEC;
	//std::cout<<"haha"<<std::endl;
	double trace;
	double truncerr;


	OP temp;
	OP dentemp;


	Supp.wave.Wave2OP(temp, Sys.SubSysEye, m.SubSysEye, n.SubSysEye, Env.SubSysEye, Gdir);

	dentemp.getDenS(temp);
	truncU.DengetTruncU(para, dentemp, trace, truncerr);//temp.show();

	temp.clear();
	Supp.wave.Wave2OP(temp, Sys.SubSysEye, m.SubSysEye, n.SubSysEye, Env.SubSysEye, -1 * Gdir);




	//Sup.Wave.Wave2OP(temp, Sys.SubSysEye, m.SubSysEye, n.SubSysEye, Env.SubSysEye, Gdir);
	dentemp.getDenE(temp);
	truncUR.DengetTruncU(para, dentemp, trace, truncerr);

	//truncU.show();
	SaveAll << "Q = " << qtot << "    WaveD = " << std::setw(4) << Sup.Dim
		<< "      OS =" << std::setw(2) << Sys.Orbital << ",  OE =" << std::setw(2) << Env.Orbital
		<< ",    E = " << std::setw(10) << std::setprecision(15) << para.Energy << ",    trace = " << std::setprecision(15) << trace
		<< ",    truncerr = " << std::setprecision(15) << truncerr << std::endl;


	std::cout << "Q = " << qtot << "    WaveD = " <<std::setw(4)<< Sup.Dim
	<< "      OS ="  <<std::setw(2)<<Sys.Orbital << ",  OE =" <<std::setw(2)<< Env.Orbital
	<< ",    E = " <<std::setw(10)<< std::setprecision(15)<<para.Energy <<",    trace = "<< std::setprecision(15)<<trace
	<<",    truncerr = "<< std::setprecision(15)<<truncerr<<std::endl;
	//FEnergy = para.Energy;
	//FTrace = trace;
	//FTruncerr = truncerr;



}


void DMRGP::truncUpdateP(const Parameter& para, int& OS, int& OE, int dir)
{

	OS += dir;
	OE -= dir;

	if (Gdir == 1)
	{
		if (m.QorRl == 0)
		{
			newS.update(para, OS, Sys, m, para.gr);
			newE.update(para, OE, Env, m, para.gl);

		}
		else
		{
			newS.update(para, OS, Sys, m, para.gl);
			newE.update(para, OE, Env, m, para.gr);
		}

	}
	else
	{
		if (n.QorRr == 0)
		{
			newS.update(para, OS, n, Sys, para.gl);
			newE.update(para, OE, n, Env, para.gr);
		}
		else
		{
			newS.update(para, OS, n, Sys, para.gr);
			newE.update(para, OE, n, Env, para.gl);
		}

	}





	//if(OS != 1)
	//{
	newS.trunc(truncU);
	newE.trunc(truncUR);
	//}

	newS.save();





	newE.save();

}





//=============sweep================
void DMRGP::SweepP(Parameter& para, int& OS, int& OE, int& dir)
{

	int flag(0);

	FEnergy = 10000000000;
	while (flag<para.SweepNo)
	{

		SaveAll << "the " << (flag + 1) << "th Sweep" << std::endl;
		//std::cout<<"the "<<(flag+1)<<"th Sweep"<<std::endl;
		//dir*=(-1);//local here for the first left direction sweep


		//FEnergy = 1000000000;
		while (true)
		{
			m.Initial(para, OrbitalM);
			n.Initial(para, OrbitalN);
			//===============consider the ways in the xishoudian d fangshi ==================================================================================
			//==========saomiao guocheng zhong d zhangdian, zai zhangdao zhongdian d shihou you yige zhangdian fangxiang d fanzhuang, yuanben yinggai ============
			//==========youbian zhangdian d shihou huancheng l zuobian zhangdian.=========================================================================
			//=======zhuyao kan Wave2OP d fangshi====================


			//==================this aprt is for the first right sweep, if first left sweep, it should absent==========
			if (OS == (para.LatticeSize - 2) / 2)
				Gdir *= -1;
			//============================present with line 356=================================

			Sys.read(OS);//Sys.show();
			Env.read(OE);//Env.show();

			OP truncU;

			//exit(1);

			getEnergySweepP(para, dir);

			if (dir == 1)
			{

				if (para.LatticeSize - OS == 3)
				{

					break;
				}
			}
			else
			{
				if (OE == 4)
				{

					break;
				}
			}


			truncUpdateSweepP(para, OS, OE, dir);

			if (dir == 1)
			{
				if (Gdir == 1)
				{
					OrbitalM += 1;
				}
				else
				{
					OrbitalN += 1;
				}
				Gdir *= -1;
			}
			else
			{

				if (Gdir == 1)
				{
					OrbitalN -= 1;
				}
				else
				{
					OrbitalM -= 1;
				}
				Gdir *= -1;
			}


			//==================this aprt is for the first left sweep, if first left sweep, it should absent==========
			/*if(OS == (para.LatticeSize -2)/2)
			Gdir *= -1;*/
			//============================present with line 277========================================


		}
		dir *= (-1);    //local the for the first right sweep
		flag++;
	}

}






void DMRGP::getEnergySweepP(Parameter& para, int dir)
{


	int qtot = para.ParticleNo;
	double trace;
	double truncerr;













	Super Sup(para, Sys, m, n, Env, qtot);






	begin = clock();
	SuperEnergy Supp(para, Sup);
	saveT += difftime(clock(), begin) / CLOCKS_PER_SEC;




	OP temp;
	OP Dentemp;
	Supp.wave.Wave2OP(temp, Sys.SubSysEye, m.SubSysEye, n.SubSysEye, Env.SubSysEye, Gdir);
	if (dir == 1)
	{
		Dentemp.getDenS(temp);
	}
	else
	{
		Dentemp.getDenE(temp);
	}
	DenOPWave.time(Dentemp, 1.0 / 3);

	Energy = para.Energy;
	para.Energy = 0;
	Super SupL(para, Sys, m, n, Env, qtot - 1);
	begin = clock();
	SuperEnergy SuppL(para, SupL);
	saveT += difftime(clock(), begin) / CLOCKS_PER_SEC;



	temp.clear();
	Dentemp.clear();
	SuppL.wave.Wave2OP(temp, Sys.SubSysEye, m.SubSysEye, n.SubSysEye, Env.SubSysEye, Gdir);
	if (dir == 1)
	{
		Dentemp.getDenS(temp);
	}
	else
	{
		Dentemp.getDenE(temp);
	}
	Dentemp.time(1.0 / 3);
	DenOPWave.addWave(Dentemp);


	LEnergy = para.Energy;


	para.Energy = 0;
	Super SupR(para, Sys, m, n, Env, qtot + 1);
	begin = clock();
	SuperEnergy SuppR(para, SupR);
	saveT += difftime(clock(), begin) / CLOCKS_PER_SEC;



	temp.clear();
	Dentemp.clear();
	SuppR.wave.Wave2OP(temp, Sys.SubSysEye, m.SubSysEye, n.SubSysEye, Env.SubSysEye, Gdir);
	if (dir == 1)
	{
		Dentemp.getDenS(temp);
	}
	else
	{
		Dentemp.getDenE(temp);
	}
	Dentemp.time(1.0 / 3);
	DenOPWave.addWave(Dentemp);


	REnergy = para.Energy;





	truncU.DengetTruncU(para, DenOPWave, trace, truncerr);//temp.show();






	//truncU.show();
	SaveAll << "Q = " << qtot << ",    E = " << std::setprecision(15) << Energy << ",    LE = " << std::setprecision(15) << LEnergy << ",    RE = " << std::setprecision(15) << REnergy << std::endl
		<< "      OS =" << std::setw(2) << Sys.Orbital << ",  OE =" << std::setw(2) << Env.Orbital
		<< ",    trace = " << std::setprecision(15) << trace
		<< ",    truncerr = " << std::setprecision(15) << truncerr << std::endl << std::endl;




	std::cout << "Q = " << qtot << ",    E = " << std::setprecision(15)<<Energy << ",    LE = " << std::setprecision(15)<<LEnergy  << ",    RE = " << std::setprecision(15)<<REnergy<<std::endl
	<< "      OS ="  <<std::setw(2)<<Sys.Orbital << ",  OE =" <<std::setw(2)<< Env.Orbital
	<<",    trace = "<< std::setprecision(15)<<trace
	<<",    truncerr = "<< std::setprecision(15)<<truncerr << std::endl<<std::endl;


	if (Sys.Orbital == (para.LatticeSize - 2) / 2)
	{
		FEnergy = Energy;
		MiuP = REnergy - Energy;
		MiuN = Energy - LEnergy;
		FTrace = trace;
		FTruncerr = truncerr;
	}



}



void DMRGP::truncUpdateSweepP(const Parameter& para, int& OS, int& OE, int dir)
{


	//std::cout<<dir<<std::endl;
	if (dir == 1)
	{

		if (Gdir == 1)
		{
			if (m.QorRl == 0)
			{
				newS.update(para, OS + 1, Sys, m, para.gr);
			}
			else
			{
				newS.update(para, OS + 1, Sys, m, para.gl);
			}
		}
		else
		{
			if (n.QorRr == 0)
			{
				newS.update(para, OS + 1, n, Sys, para.gl);
			}
			else
			{
				newS.update(para, OS + 1, n, Sys, para.gr);
			}
		}
		//newS.show();
		//truncU.show();

		newS.trunc(truncU);
		//newS.show();

		newS.save();
	}
	else
	{

		if (Gdir == 1)
		{
			if (n.QorRr == 0)
			{
				newS.update(para, OE - 1, n, Env, para.gr);
			}
			else
			{
				newS.update(para, OE - 1, n, Env, para.gl);
			}
		}
		else
		{
			if (m.QorRl == 0)
			{
				newS.update(para, OE - 1, Env, m, para.gl);
			}
			else
			{
				newS.update(para, OE - 1, Env, m, para.gr);
			}
		}




		newS.trunc(truncU);
		//truncUR.show();

		newS.save();
	}
	OE += dir;
	OS += dir;
}
