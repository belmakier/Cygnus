#include "DataReader.h"

DataReader::DataReader(Nucleus* nucl, const char *datafilename, bool GOSIA){
	exptData.clear();
	fNucleus = *nucl;
	if(GOSIA)
		ReadGOSIAData(datafilename);
	else
		ReadDataFile(datafilename);
}

void	DataReader::ReadDataFile(const char* datafilename){

	std::ifstream infile(datafilename);

	int nExpt = 0;

	int initial_state;
	int final_state;
	double counts;
	double uncertainty;

	ExperimentData tmpExpt;

	std::string line;
	while(std::getline(infile,line)){
		std::size_t found = line.find("!"); // Comment
		if(found == std::string::npos){
			std::istringstream ss(line);
			found = line.find("EXPT");  // Experiment definition
			if(found != std::string::npos && nExpt !=0){
				exptData.push_back(tmpExpt);	
				tmpExpt.ClearData();
				nExpt++;
				continue;
			}
			else if(found != std::string::npos){
				nExpt++;
				continue;
			}
			found = line.find("END");
			if(found != std::string::npos){
				exptData.push_back(tmpExpt);
				tmpExpt.ClearData();
				break;
			}
			ss >> initial_state >> final_state >> counts >> uncertainty;
			if(initial_state > 100 && final_state > 100){	// Doublet
				int	i1	= initial_state / 100;
				int	f1	= final_state / 100;
				int 	i2	= initial_state - i1 * 100;			
				int 	f2	= final_state - f1 * 100;
				std::cout	<< "Doublet:"
						<< std::endl;
				std::cout	<< i1	<< "\t" << f1 << "\t" << i2 << "\t" << f2	
						<< std::endl;
				tmpExpt.AddDoublet(i1,f1,i2,f2,counts,uncertainty);
			}
			else if(initial_state < fNucleus.GetNstates() && final_state < fNucleus.GetNstates() && initial_state >= 0 && final_state >= 0){
				tmpExpt.AddData(initial_state,final_state,counts,uncertainty);
			}
			else{
				if(initial_state < fNucleus.GetNstates())
					std::cout << "State " << initial_state << " out of range" << std::endl;
				if(final_state < fNucleus.GetNstates())
					std::cout << "State " << final_state << " out of range" << std::endl;
				fNucleus.PrintNucleus();
			}
		}
	}	

}


void	DataReader::ReadGOSIAData(const char* datafilename){

	std::ifstream infile(datafilename);

	int nExpt = 0;

	int initial_state;
	int final_state;
	double counts;
	double uncertainty;

	ExperimentData tmpExpt;

	int	lineCounter = 0;
	int	exptCounter = 0;
	int	dataCounter = 0;
	int	nData;
	bool	Flag = false;

	std::string line;
	while(std::getline(infile,line)){
		std::istringstream ss(line);
		if(lineCounter == 0){ // First line
			int	tmp1,tmp2,tmp3,tmp4,tmp5,tmp6;
			float	tmp7;
			ss	>> tmp1 >> tmp2 >> tmp3 >> tmp4 >> tmp5 >> tmp6 >> tmp7;
			nData	= tmp6;
		}
		else if(dataCounter < nData){ // Data
			int	tmp1,tmp2;
			float	tmp3,tmp4;
			ss	>> tmp1 >> tmp2 >> tmp3 >> tmp4;
			dataCounter++;
			tmpExpt.AddData(tmp1-1,tmp2-1,tmp3,tmp4);
		}
		else{	// New experiment
			exptData.push_back(tmpExpt);
			tmpExpt.ClearData();
			exptCounter++;
			int	tmp1,tmp2,tmp3,tmp4,tmp5,tmp6;
			float	tmp7;
			ss	>> tmp1 >> tmp2 >> tmp3 >> tmp4 >> tmp5 >> tmp6 >> tmp7;
			nData	= tmp6;
			dataCounter = 0;
		}
		lineCounter++;
	}	

	if(tmpExpt.GetData().size() > 0){
		exptData.push_back(tmpExpt);
		tmpExpt.ClearData();
	}
}



