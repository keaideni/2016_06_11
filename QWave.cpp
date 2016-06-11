#include "QWave.h"


QWave::QWave(){}



QWave::~QWave(){}




QWave::QWave(const QWave& wave)
{
	WavePart = wave.WavePart;
}






QWave::QWave(const OP& Sys, const OP& m, const OP& n, const OP& Env, int QTot)
{
	for (auto itm = m.QDim.begin(); itm != m.QDim.end(); itm++)
	{
		for (auto itn = n.QDim.begin(); itn != n.QDim.end(); itn++)
		{
			OP tempOP;
			for (auto its = Sys.QDim.begin(); its != Sys.QDim.end(); its++)
			{

				for (auto ite = Env.QDim.begin(); ite != Env.QDim.end(); ite++)
				{
					int q1, q2, qtot_;
					q1=(itm->first+ itn->first);
					q2=(its->first+ ite->first);
					qtot_=(q1+ q2);
					if (qtot_ == QTot)
					{
						tempOP.RLQ[ite->first] = its->first;
						int dimL = its->second;
						int dimR = ite->second;
						MatrixXd tempM(MatrixXd::Zero(dimL, dimR));
						
						tempOP.QMat[ite->first] = tempM;

					}
				}
			}
			if (tempOP.QMat.size()>0)
			{
				WavePart[std::pair<int, int>(itm->first, itn->first)] = tempOP;
			}

		}
	}
}






void QWave::Initial(const OP& Sys, const OP& m, const OP& n, const OP& Env, int QTot)
{
	for (auto itm = m.QDim.begin(); itm != m.QDim.end(); itm++)
	{
		for (auto itn = n.QDim.begin(); itn != n.QDim.end(); itn++)
		{
			OP tempOP;
			for (auto its = Sys.QDim.begin(); its != Sys.QDim.end(); its++)
			{

				for (auto ite = Env.QDim.begin(); ite != Env.QDim.end(); ite++)
				{
					int q1, q2, qtot_;
					q1=(itm->first + itn->first);
					q2=(its->first + ite->first);
					qtot_=(q1 + q2);
					if (qtot_ == QTot)
					{
						tempOP.RLQ[ite->first] = its->first;
						int dimL = its->second;
						int dimR = ite->second;
						MatrixXd tempM(MatrixXd::Zero(dimL, dimR));
						tempOP.QMat[ite->first] = tempM;

					}
				}
			}
			if (tempOP.QMat.size()>0)
			{
				WavePart[std::pair<int, int>(itm->first, itn->first)] = tempOP;
			}

		}
	}
}







void QWave::setZero()
{
	for (auto it = WavePart.begin(); it != WavePart.end(); it++)
	{
		for (auto itt = (it->second).QMat.begin(); itt != (it->second).QMat.end(); itt++)
		{
			int diml = (itt->second).rows();
			int dimr = (itt->second).cols();
			itt->second = (MatrixXd::Zero(diml, dimr));
		}
	}
}






void QWave::normalize()
{
	double x(0);
	for (auto it = WavePart.begin(); it != WavePart.end(); it++)
	{
		for (auto tempit = it->second.QMat.begin(); tempit != it->second.QMat.end(); tempit++)
		{
			int diml = tempit->second.rows();
			int dimr = tempit->second.cols();
			for (int i = 0; i<diml; i++)
			{
				for (int j = 0; j<dimr; j++)
				{
					x +=  ((tempit->second)(i, j)) * ((tempit->second)(i, j));
				}
			}
		}
	}
	double y = sqrt(x);


	for (auto it = WavePart.begin(); it != WavePart.end(); it++)
	{
		for (auto tempit = it->second.QMat.begin(); tempit != it->second.QMat.end(); tempit++)
		{
			tempit->second /= y;
		}
	}
}






int QWave::getDim() const
{
	int n(0);

	for (auto it = WavePart.begin(); it != WavePart.end(); it++)
	{
		for (auto tempit = it->second.QMat.begin(); tempit != it->second.QMat.end(); tempit++)
		{
			n += (tempit->second.rows())*(tempit->second.cols());
		}
	}

	return n;
}





void QWave::add(const QWave&a, const QWave& b)
{
	WavePart = a.WavePart;
	for (auto it = b.WavePart.begin(); it != b.WavePart.end(); it++)
	{
		auto tempit = a.WavePart.find(it->first);
		if (tempit == a.WavePart.end())
		{
			WavePart.insert(std::pair<std::pair<int, int>, OP>(it->first, it->second));
		}
		else
		{
			OP tempOP;
			tempOP.addWave(it->second, tempit->second);
			WavePart.at(it->first) = tempOP;
		}
	}
}







void QWave::Wave2f(std::vector<double>& f) const
{
	f.clear();
	for (auto it = WavePart.begin(); it != WavePart.end(); it++)
	{
		for (auto tempit = it->second.QMat.begin(); tempit != it->second.QMat.end(); tempit++)
		{
			int diml = tempit->second.rows();
			int dimr = tempit->second.cols();
			for (int i = 0; i<diml; i++)
			{
				for (int j = 0; j < dimr; j++)
				{
					double tempv = (tempit->second)(i, j);
					f.push_back(tempv);
				}
			}
		}
	}
}







void QWave::f2Wave(const std::vector<double>& f)
{
	int n(0);
	for (auto it = WavePart.begin(); it != WavePart.end(); it++)
	{
		for (auto tempit = it->second.QMat.begin(); tempit != it->second.QMat.end(); tempit++)
		{
			int diml = tempit->second.rows();
			int dimr = tempit->second.cols();
			for (int i = 0; i<diml; i++)
			{
				for (int j = 0; j < dimr; j++)
				{
					(tempit->second)(i, j) = f.at(n);
					n++;
				}
			}
		}
	}
}



void QWave::f2Wave(const VectorXd& f)
{
	int n(0);
	for (auto it = WavePart.begin(); it != WavePart.end(); it++)
	{
		for (auto tempit = it->second.QMat.begin(); tempit != it->second.QMat.end(); tempit++)
		{
			int diml = tempit->second.rows();
			int dimr = tempit->second.cols();
			for (int i = 0; i<diml; i++)
			{
				for (int j = 0; j < dimr; j++)
				{
					(tempit->second)(i, j) = f(n);
					n++;
				}
			}
		}
	}
}







//translate the QWave to OP
void QWave::Wave2OP(OP& O, const OP& sys, const OP& m, const OP& n, const OP& env) const
{
	OP tempS, tempE;
	std::map<std::pair<int, int>, int, classcom> startDimS, startDimE;
	tempS.findDim(sys, m, tempS.QDim, startDimS);
	tempE.findDim(env, n, tempE.QDim, startDimE);
	//tempS.show();

	for (auto it = WavePart.begin(); it != WavePart.end(); it++)
	{
		for (auto itt = it->second.RLQ.begin(); itt != it->second.RLQ.end(); itt++)
		{
			int tempSQ, tempEQ;
			tempSQ = (it->first.first + itt->second);
			tempEQ = (it->first.second + itt->first);

			int dimS = tempS.QDim.at(tempSQ);
			int dimE = tempE.QDim.at(tempEQ);

			MatrixXd tempmat(MatrixXd::Zero(dimS, dimE));
			
			

			int startL = startDimS.at(std::pair<int, int>(itt->second, it->first.first));
			int startR = startDimE.at(std::pair<int, int>(itt->first, it->first.second));

			int size_row = it->second.QMat.at(itt->first).rows();
			int size_col = it->second.QMat.at(itt->first).cols();

			




			auto temppp = O.QMat.find(tempEQ);
			if (temppp != O.QMat.end())
			{
				tempmat = temppp->second;
			}
			
			//tempmat.block(startL, startR, size_row, size_col) = it->second.QMat.at(itt->first);

			for (int i = 0; i < size_row; ++i)
			{
				for (int j = 0; j < size_col; ++j)
				{
					tempmat(startL + i, startR + j) = it->second.QMat.at(itt->first)(i, j);
				}
			}

			O.QMat[tempEQ] = tempmat;


			auto tempp = O.RLQ.find(tempEQ);
			if (tempp == O.RLQ.end())
			{
				O.RLQ[tempEQ] = tempSQ;
			}

		}
	}

}






//translate the QWave to OP
void QWave::Wave2OP(OP& O, const OP& sys, const OP& m, const OP& n, const OP& env, const int& way) const
{
	O.clear();
	OP tempS, tempE;
	std::map<std::pair<int, int>, int, classcom> startDimS, startDimE;
	switch (way)
	{
	case 1://==============add m on the right edge of system=============
	{
		       tempS.findDim(sys, m, tempS.QDim, startDimS);
		       tempE.findDim(n, env, tempE.QDim, startDimE);
		       break;
	}
	case -1://==========add n on the left edge of system=============
	{
			tempS.findDim(n, sys, tempS.QDim, startDimS);
			tempE.findDim(env, m, tempE.QDim, startDimE);

			break;
	}
	}
	for (auto it = WavePart.begin(); it != WavePart.end(); it++)
	{
		for (auto itt = it->second.RLQ.begin(); itt != it->second.RLQ.end(); itt++)
		{
			int tempSQ, tempEQ;


			int startL, startR;
			switch (way)
			{
			case 1://==============add m on the right edge of system=============
			{

				       tempSQ = (it->first.first + itt->second);
				       tempEQ = (it->first.second + itt->first);

				       startL = startDimS.at(std::pair<int, int>(itt->second, it->first.first));
				       startR = startDimE.at(std::pair<int, int>(it->first.second, itt->first));

				       break;


			}
			case -1://==========add n on the left edge of system=============
			{
					tempSQ = (itt->second + it->first.second);
					tempEQ = (itt->first + it->first.first);

					startL = startDimS.at(std::pair<int, int>(it->first.second, itt->second));
					startR = startDimE.at(std::pair<int, int>(itt->first, it->first.first));


					break;

			}

			}


			int dimS = tempS.QDim.at(tempSQ);
			int dimE = tempE.QDim.at(tempEQ);

			MatrixXd tempmat(MatrixXd::Zero(dimS, dimE));
			
			

			

			int size_row = it->second.QMat.at(itt->first).rows();
			int size_col = it->second.QMat.at(itt->first).cols();

			




			auto temppp = O.QMat.find(tempEQ);
			if (temppp != O.QMat.end())
			{
				tempmat = temppp->second;
			}
			
			tempmat.block(startL, startR, size_row, size_col) = it->second.QMat.at(itt->first);
			O.QMat[tempEQ] = tempmat;


			auto tempp = O.RLQ.find(tempEQ);
			if (tempp == O.RLQ.end())
			{
				O.RLQ[tempEQ] = tempSQ;
			}

		}
	}

}




//the reroad of some operator.
void QWave::operator=(const QWave& wave)
{
	WavePart.clear();
	WavePart = wave.WavePart;

}



QWave QWave::operator+(const QWave& wave)
{
	QWave sumwave;
	sumwave.WavePart = WavePart;
	for (auto it = wave.WavePart.begin(); it != wave.WavePart.end(); it++)
	{
		auto tempit = WavePart.find(it->first);
		if (tempit == WavePart.end())
		{
			sumwave.WavePart.insert(std::pair<std::pair<int, int>, OP>(it->first, it->second));
		}
		else
		{
			OP tempOP;
			tempOP.addWave(it->second, tempit->second);
			sumwave.WavePart.at(it->first) = tempOP;
		}
	}
	return sumwave;

}






void QWave::show()
{
	for (auto it = WavePart.begin(); it != WavePart.end(); it++)
	{
		std::cout << "< " << it->first.first << ", " << it->first.second << "> :" << std::endl;
		it->second.show();
	}
}





void QWave::clear()
{
	WavePart.clear();

}






//==============operator function on the QWave 1============================
//|WavePart> = O |wave>
void QWave::OPWave2New(const QWave& wave, const OP& O, int flag)
{
	switch (flag)
	{
	case 1:
	{
		      OSWave2New(O, wave);
		      break;
	}
	case 2:
	{
		      OMWave2New(O, wave);
		      break;
	}
	case 3:
	{
		      ONWave2New(O, wave);
		      break;
	}
	case 4:
	{
		      OEWave2New(O, wave);
		      break;
	}
	}
}






void QWave::OSWave2New(const OP& O, const QWave& wave)
{
	if (O.RLQ.size() == 0)
	{
		exit(1);
	}

	for (auto it = wave.WavePart.begin(); it != wave.WavePart.end(); it++)
	{
		OP temp;
		if (temp.ltime(O, it->second) > 0)
		{
			WavePart.insert(std::pair<std::pair<int, int>, OP>(it->first, temp));
		}
	}
}






void QWave::OEWave2New(const OP& O, const QWave& wave)
{
	if (O.RLQ.size() == 0)
	{
		exit(1);
	}

	for (auto it = wave.WavePart.begin(); it != wave.WavePart.end(); it++)
	{
		OP temp;
		if (temp.rtime(O, it->second) > 0)
		{
			WavePart.insert(std::pair<std::pair<int, int>, OP>(it->first, temp));
		}
	}
}







void QWave::OMWave2New(const OP& O, const QWave& wave)
{
	if (O.RLQ.size() == 0)
	{
		exit(1);
	}

	auto it = O.RLQ.begin();
	int dQ = it->second - it->first;

	for (auto tempit = wave.WavePart.begin(); tempit != wave.WavePart.end(); tempit++)
	{
		auto itt = O.QMat.find(tempit->first.first);
		if (itt != O.QMat.end())
		{
			double x = (itt->second)(0, 0);
			int tempQ(tempit->first.first + dQ);
			OP temp;
			temp.time(x, tempit->second);
			WavePart.insert(std::pair<std::pair<int, int>, OP>(std::pair<int, int>(tempQ, tempit->first.second), temp));
		}
	}
}







void QWave::ONWave2New(const OP& O, const QWave& wave)
{
	if (O.RLQ.size() == 0)
	{
		exit(1);
	}

	auto it = O.RLQ.begin();
	int dQ(it->second - it->first);

	for (auto tempit = wave.WavePart.begin(); tempit != wave.WavePart.end(); tempit++)
	{
		auto itt = O.QMat.find(tempit->first.second);
		if (itt != O.QMat.end())
		{
			double x = (itt->second)(0, 0);
			int tempQ(tempit->first.second + dQ);
			OP temp;
			temp.time(x, tempit->second);
			WavePart.insert(std::pair<std::pair<int, int>, OP>(std::pair<int, int>(tempit->first.first, tempQ), temp));
		}
	}
}







//==================operator function on the QWave 2====================
//|storewave> = O|wavepart> + |storewave>         ==============why it has the second term +|storewave>?==================
void QWave::OPWave(QWave& wave, const OP& O, int flag) const
{
	switch (flag)
	{
	case 1:
	{
		      OSWave(O, wave);
		      break;
	}
	case 2:
	{
		      OMWave(O, wave);
		      break;
	}
	case 3:
	{
		      ONWave(O, wave);
		      break;
	}
	case 4:
	{
		      OEWave(O, wave);
		      break;
	}
	}
}









void QWave::OSWave(const OP& O, QWave& storewave) const
{
	if (O.RLQ.size() == 0)
	{
		exit(1);
	}

	for (auto it = WavePart.begin(); it != WavePart.end(); it++)
	{
		auto itt = storewave.WavePart.find(it->first);
		if (itt != storewave.WavePart.end())
		{
			OP temp;
			if (temp.ltime(O, it->second) > 0)
			{
				itt->second.addWave(temp);
			}
		}
		else
		{
			OP temp;
			if (temp.ltime(O, it->second)>0)
			{
				storewave.WavePart.insert(std::pair<std::pair<int, int>, OP>(it->first, temp));
			}
		}
	}
}






void QWave::OEWave(const OP& O, QWave& storewave) const
{
	if (O.RLQ.size() == 0)
	{
		exit(1);
	}

	for (auto it = WavePart.begin(); it != WavePart.end(); it++)
	{
		auto itt = storewave.WavePart.find(it->first);
		if (itt != storewave.WavePart.end())
		{
			OP temp;
			if (temp.rtime(O, it->second) > 0)
			{
				itt->second.addWave(temp);
			}
		}
		else
		{
			OP temp;
			if (temp.rtime(O, it->second)>0)
			{
				storewave.WavePart.insert(std::pair<std::pair<int, int>, OP>(it->first, temp));
			}
		}
	}
}








void QWave::OMWave(const OP& O, QWave& storewave) const
{
	if (O.RLQ.size() == 0)
	{
		exit(1);
	}

	auto it = O.RLQ.begin();
	int dQ(it->second - it->first);

	QWave tempQW;
	for (auto tempit = WavePart.begin(); tempit != WavePart.end(); tempit++)
	{
		auto itt = O.QMat.find(tempit->first.first);
		if (itt != O.QMat.end())
		{
			double x = (itt->second)(0, 0);
			int tempQ(tempit->first.first + dQ);
			OP temp;
			temp.time(x, tempit->second);
			tempQW.WavePart.insert(std::pair<std::pair<int, int>, OP>(std::pair<int, int>(tempQ, tempit->first.second), temp));
		}
	}
	storewave = tempQW + storewave;
}








void QWave::ONWave(const OP& O, QWave& storewave) const
{
	if (O.RLQ.size() == 0)
	{
		exit(1);
	}

	auto it = O.RLQ.begin();
	int dQ(it->second - it->first);

	QWave tempQW;
	for (auto tempit = WavePart.begin(); tempit != WavePart.end(); tempit++)
	{
		auto itt = O.QMat.find(tempit->first.second);
		if (itt != O.QMat.end())
		{
			double x = (itt->second)(0, 0);
			int tempQ(tempit->first.second + dQ);
			OP temp;
			temp.time(x, tempit->second);
			tempQW.WavePart.insert(std::pair<std::pair<int, int>, OP>(std::pair<int, int>(tempit->first.first, tempQ), temp));
		}
	}
	storewave = tempQW + storewave;
}









