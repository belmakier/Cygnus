#include "ExperimentalInput.h"

ExptDoublet::ExptDoublet(const ExptDoublet& e){

 	StateIndex_I1	= e.StateIndex_I1;
  	StateIndex_F1	= e.StateIndex_F1;
 	StateIndex_I2	= e.StateIndex_I2;
  	StateIndex_F2	= e.StateIndex_F2;
	Counts		= e.Counts;
	UpUncertainty 	= e.UpUncertainty;
	DnUncertainty	= e.DnUncertainty;
	DetEff		= e.DetEff;	

}
ExptDoublet& ExptDoublet::operator = (const ExptDoublet& e){

 	StateIndex_I1	= e.StateIndex_I1;
  	StateIndex_F1	= e.StateIndex_F1;
 	StateIndex_I2	= e.StateIndex_I2;
  	StateIndex_F2	= e.StateIndex_F2;
	Counts		= e.Counts;
	UpUncertainty 	= e.UpUncertainty;
	DnUncertainty	= e.DnUncertainty;
	DetEff		= e.DetEff;

	return *this;

}

ExptData::ExptData(const ExptData& e){

 	StateIndex_I	= e.StateIndex_I;
  	StateIndex_F	= e.StateIndex_F;
	Counts		= e.Counts;
	UpUncertainty 	= e.UpUncertainty;
	DnUncertainty	= e.DnUncertainty;
	DetEff		= e.DetEff;	

}
ExptData& ExptData::operator = (const ExptData& e){

 	StateIndex_I	= e.StateIndex_I;
  	StateIndex_F	= e.StateIndex_F;
	Counts		= e.Counts;
	UpUncertainty 	= e.UpUncertainty;
	DnUncertainty	= e.DnUncertainty;
	DetEff		= e.DetEff;

	return *this;

}
ExperimentData::ExperimentData(const ExperimentData& e){

	Data.resize(e.Data.size());
	for(size_t d=0;d<e.Data.size();d++)
		Data.at(d) = e.Data.at(d);		
	Doublet.resize(e.Doublet.size());
	for(size_t d=0;d<e.Doublet.size();d++)
		Doublet.at(d) = e.Doublet.at(d);		
	thetaCM		= e.thetaCM;
	

}
ExperimentData& ExperimentData::operator = (const ExperimentData& e){

	Data.resize(e.Data.size());
	for(size_t d=0;d<e.Data.size();d++)
		Data.at(d) = e.Data.at(d);		
	Doublet.resize(e.Doublet.size());
	for(size_t d=0;d<e.Doublet.size();d++)
		Doublet.at(d) = e.Doublet.at(d);		
	thetaCM		= e.thetaCM;

	return *this;

}

void ExperimentData::Print() const {

	if(thetaCM < 0)
		std::cout << "Theta undefined" << std::endl;
	else
		std::cout << "Theta [CM]: " << thetaCM << " degrees" << std::endl;

	std::cout 	<< std::setw(9) << std::left << "Initial"
			<< std::setw(7) << std::left << "Final"
			<< std::setw(10) << std::left << "Counts"
			<< std::setw(10) << std::left << "UpUncert."
			<< std::setw(10) << std::left << "DnUncert."
			<< std::endl;
	for(unsigned int i=0;i<Data.size();i++)
		Data.at(i).Print();

}

void ExptData::Print() const {
	std::cout 	<< std::setw(9) << std::left << StateIndex_I 
			<< std::setw(7) << std::left << StateIndex_F 
			<< std::setw(10) << std::left << Counts 
			<< std::setw(10) << UpUncertainty 
			<< std::setw(10) << DnUncertainty 
			<< std::endl; 
}
