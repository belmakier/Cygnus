#include "NucleusReader.h"

NucleusReader::NucleusReader() : fNucleus(NULL)
{
}

NucleusReader::NucleusReader(const char *filename, bool GOSIA) : fNucleus(NULL)
{
	if(GOSIA)
		ReadGOSIANucleus(filename);
	else
		ReadNucleusFile(filename);	
}

NucleusReader::NucleusReader(const NucleusReader &n) : fNucleus(n.fNucleus)
{
}

NucleusReader& NucleusReader::operator = (const NucleusReader &n)
{

	fNucleus	= n.fNucleus;

	return *this;

}

void NucleusReader::ReadNucleusFile(const char *filename){
	
	std::ifstream infile(filename);
	
	int counter = 0;
	int tmpA, tmpZ, tmpS, tmpL;
	int tmpI, tmpP; 
	double tmpE, tmpJ;

	int tmpI1, tmpI2, tmp_Lambda;
	std::string tmp_Lambda_s;
	double tmp_ME;

	std::string line;
	while(std::getline(infile,line)){
		std::size_t found = line.find("!");
		if(found == std::string::npos){
			std::istringstream ss(line);

			if(counter == 0){
				ss >> tmpA >> tmpZ >> tmpS >> tmpL;
				fNucleus = new Nucleus(tmpA,tmpZ,tmpS,tmpL);
			}
			else if(counter <= tmpS){
				ss >> tmpI >> tmpE >> tmpJ >> tmpP;
				fNucleus->SetState(tmpI,tmpE,tmpJ,tmpP);	
			}
			else{
				std::size_t foundE = line.find("E");
				std::size_t foundM = line.find("M");
				if(foundE != std::string::npos || foundM != std::string::npos){ // Found a transition specified E or M
					ss >> tmpI1 >> tmpI2 >> tmp_ME >> tmp_Lambda_s;
					if(foundE != std::string::npos && foundM !=std::string::npos){
						std::cout << "Both E and M found, not possible. Skipping line!" << std::endl;
						continue;
					}
					tmp_Lambda = tmp_Lambda_s[1] - '0'; // Since we know the second digit is the multipolarity
					if(foundM != std::string::npos)
						tmp_Lambda += 6;	
				}
				else{
					ss >> tmpI1 >> tmpI2 >> tmp_ME >> tmp_Lambda;	
				}
				fNucleus->SetMatrixElement(tmp_Lambda-1,tmpI1,tmpI2,tmp_ME);
			}

			counter++;
		}

	}

	//fNucleus->PrintNucleus();

}

void NucleusReader::ReadGOSIANucleus(const char* filename){

	std::ifstream infile(filename);

	std::string line;

	std::vector<int>	tI;
	std::vector<int>	tP;
	std::vector<double>	tJ;
	std::vector<double>	tE;

	int	nExpt;
	int	tmpZ;
	int	tmpA;

	bool	exptFlag = false;
	while(std::getline(infile,line)){
		if(exptFlag){
			std::istringstream ss(line);
			ss	>> nExpt >> tmpZ >> tmpA;	
			break;	// That's all we need from here
		}
		if(line == "EXPT")
			exptFlag = true;
	}

	infile.clear();
	infile.seekg(0);

	bool	levelFlag = false;
	while(std::getline(infile,line)){

		if(levelFlag && (line == "0,0,0,0" || line == "0\t0\t0\t0")){
			levelFlag = false;
		}

		if(levelFlag){
			std::cout	<< line
					<< std::endl;
			std::istringstream ss(line);
			int	tmpI,tmpP;
			double	tmpJ,tmpE;
			ss	>> tmpI >> tmpP >> tmpJ >> tmpE;
			tI.push_back(tmpI-1);
			tP.push_back(tmpP);
			tJ.push_back(tmpJ);
			tE.push_back(tmpE);		
		}

		if(line == "LEVE"){
			std::cout	<< line
					<< std::endl;
			levelFlag = true;
		}
	}

	fNucleus = new Nucleus(tmpZ,tmpA,tI.size());
	for(size_t i=0;i<tI.size();i++)
		fNucleus->SetState(tI.at(i),tE.at(i),tJ.at(i),tP.at(i));

	infile.clear();
	infile.seekg(0);

	bst_i.clear();
	bst_f.clear();
	bst_l.clear();

	int	tmpLambda = -1;
	bool	MEFlag = false;
	while(std::getline(infile,line)){
		if(MEFlag && (line == "0,0,0,0,0" || line == "0\t0\t0\t0\t0"))
			MEFlag = false;

		if(MEFlag && (line == "1,0,0,0,0" || line == "1\t0\t0\t0\t0")){
			tmpLambda = 0;
		}
		else if(MEFlag && (line == "2,0,0,0,0" || line == "2\t0\t0\t0\t0")){
			tmpLambda = 1;
		}
		else if(MEFlag && (line == "3,0,0,0,0" || line == "3\t0\t0\t0\t0")){
			tmpLambda = 2;
		}
		else if(MEFlag && (line == "4,0,0,0,0" || line == "4\t0\t0\t0\t0")){
			tmpLambda = 3;
		}
		else if(MEFlag && (line == "5,0,0,0,0" || line == "5\t0\t0\t0\t0")){
			tmpLambda = 4;
		}
		else if(MEFlag && (line == "6,0,0,0,0" || line == "6\t0\t0\t0\t0")){
			tmpLambda = 5;
		}
		else if(MEFlag && (line == "7,0,0,0,0" || line == "7\t0\t0\t0\t0")){
			tmpLambda = 6;
		}
		else if(MEFlag){
			std::istringstream ss(line);
			int	tmpI,tmpF;
			float	tmpME,tmpLL,tmpUL;
			ss	>> tmpI >> tmpF >> tmpME >> tmpLL >> tmpUL;
			fNucleus->SetMatrixElement(tmpLambda,tmpI-1,tmpF-1,tmpME);
			bst_i.push_back(tmpI-1);
			bst_f.push_back(tmpF-1);
			bst_l.push_back(tmpLambda);
		}
			
		if(line == "ME")
			MEFlag = true;
	}


}
