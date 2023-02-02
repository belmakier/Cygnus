#include "Integral.h"

void ThreadTaskIntegral(PointCoulEx **pm, const std::vector<int> &mEs, const std::vector<int> &mTs, const std::vector<std::vector<double> > &tm){
  for (int i = 0; i<mEs.size(); ++i) {
    pm[mEs[i]*tm.at(0).size() + mTs[i]]->CalculatePointProbabilities(tm[mEs[i]][mTs[i]]);
  }

	return;
	
}

Integral::Integral() : fNucleus(NULL), fReaction(NULL){

	fAccuracy	= 1e-5;

	fProjectileExcitation	= true;
	targetDetection = false;
	verbose		= false;

	theta_meshpoints.clear();
	energy_meshpoints.clear();
		
	point_calculations.clear();	

	nThreads 	= 1;

	fComplete	= false;

	fUseFixedStep	= false;
	fUseSymmetry	= true;

}

Integral::Integral(Nucleus *nucl, Reaction *reac) : fNucleus(nucl), fReaction(reac) {

	fProjectileExcitation	= true;
	targetDetection = false;
	verbose		= false;

	theta_meshpoints.clear();
	energy_meshpoints.clear();
		
	point_calculations.clear();

	nThreads	= 1;
	fAccuracy	= 1e-5;

	fComplete	= false;

	fUseFixedStep	= false;
	fUseSymmetry	= true;

}

Integral::Integral(const Integral& x) : 
	fNucleus(x.fNucleus), fReaction(x.fReaction), 
	point_calculations(x.point_calculations.size()), energymeshpoint_reaction(x.energymeshpoint_reaction.size()), 
	meshpointCrossSections(x.meshpointCrossSections.size()), meshpointProbabilities(x.meshpointProbabilities.size())	
{

	fUseFixedStep		= x.fUseFixedStep;
	fUseSymmetry		= x.fUseSymmetry;

	fProjectileExcitation	= x.fProjectileExcitation;

	nThreads	= x.nThreads;
	fAccuracy	= x.fAccuracy;

	verbose		= x.verbose;
	targetDetection	= x.targetDetection;

	fComplete	= x.fComplete;

	for(size_t i = 0; i<x.point_calculations.size(); i++)
		point_calculations.at(i) = x.point_calculations.at(i);
	for(size_t i = 0; i<x.energymeshpoint_reaction.size(); i++)
		energymeshpoint_reaction.at(i) = x.energymeshpoint_reaction.at(i);
	for(size_t i = 0; i<x.meshpointCrossSections.size(); i++){
		meshpointCrossSections.at(i).resize(x.meshpointCrossSections.at(i).size());
		for(size_t j=0;j<x.meshpointCrossSections.at(i).size();j++){
			meshpointCrossSections.at(i).at(j).ResizeTo(x.meshpointCrossSections.at(i).at(j));
			meshpointCrossSections.at(i).at(j) = x.meshpointCrossSections.at(i).at(j);
		}
	}
	for(size_t i = 0; i<x.meshpointProbabilities.size(); i++){
		meshpointProbabilities.at(i).resize(x.meshpointProbabilities.at(i).size());
		for(size_t j=0;j<x.meshpointProbabilities.at(i).size();j++){
			meshpointProbabilities.at(i).at(j).ResizeTo(x.meshpointProbabilities.at(i).at(j));
			meshpointProbabilities.at(i).at(j) = x.meshpointProbabilities.at(i).at(j);
		}
	}

}
Integral& Integral::operator = (const Integral& x){

	fUseFixedStep		= x.fUseFixedStep;
	fUseSymmetry		= x.fUseSymmetry;

	fProjectileExcitation	= x.fProjectileExcitation;

	nThreads	= x.nThreads;
	fAccuracy	= x.fAccuracy;

	verbose		= x.verbose;
	targetDetection	= x.targetDetection;

	fNucleus 	= x.fNucleus;
	fReaction 	= x.fReaction;

	fComplete	= x.fComplete;
	
	point_calculations.resize(x.point_calculations.size());
	energymeshpoint_reaction.resize(x.energymeshpoint_reaction.size());
	meshpointCrossSections.resize(x.meshpointCrossSections.size());
	meshpointProbabilities.resize(x.meshpointProbabilities.size());

	for(size_t i = 0; i<x.point_calculations.size(); i++)
		point_calculations.at(i) = x.point_calculations.at(i);
	for(size_t i = 0; i<x.energymeshpoint_reaction.size(); i++)
		energymeshpoint_reaction.at(i) = x.energymeshpoint_reaction.at(i);
	for(size_t i = 0; i<x.meshpointCrossSections.size(); i++){
		meshpointCrossSections.at(i).resize(x.meshpointCrossSections.at(i).size());
		for(size_t j=0;j<x.meshpointCrossSections.at(i).size();j++){
			meshpointCrossSections.at(i).at(j).ResizeTo(x.meshpointCrossSections.at(i).at(j));
			meshpointCrossSections.at(i).at(j) = x.meshpointCrossSections.at(i).at(j);
		}
	}
	for(size_t i = 0; i<x.meshpointProbabilities.size(); i++){
		meshpointProbabilities.at(i).resize(x.meshpointProbabilities.at(i).size());
		for(size_t j=0;j<x.meshpointProbabilities.at(i).size();j++){
			meshpointProbabilities.at(i).at(j).ResizeTo(x.meshpointProbabilities.at(i).at(j));
			meshpointProbabilities.at(i).at(j) = x.meshpointProbabilities.at(i).at(j);
		}
	}

	return *this;

}

void Integral::SetThetaMeshpoint(int i, double t){

	if(i<(int)theta_meshpoints.size()){
		theta_meshpoints[i] = t;
	}
	else{
		std::cout << "Energy meshpoint index " << i << " too large, appending to list" << std::endl;
		AddThetaMeshpoint(t);
	}

}
void Integral::SetEnergyMeshpoint(int i, double e){

	if(i<(int)energy_meshpoints.size()){
		energy_meshpoints[i] = e;
	}
	else{
		std::cout << "Theta meshpoint index " << i << " too large, appending to list" << std::endl;
		AddEnergyMeshpoint(e);
	}

}

void Integral::CalculateIntegral(){  

	cmTheta.clear();	
	meshpointCrossSections.clear();
	meshpointProbabilities.clear();

	TVectorD out;
	out.ResizeTo(fNucleus->GetNstates());

	int part = 2;
	if(targetDetection)
		part = 3;

	std::vector<int> index;

	//	To simplify things for multithreading, here we pre-create
	//	all of the vectors and matrices containing the point-
	//	calculations.

	energymeshpoint_reaction.clear();
  PointCoulEx *PoinMat[energy_meshpoints.size()][theta_meshpoints.size()];
	std::vector< std::vector <double> >		ThetaMat;

	for(unsigned int mE = 0; mE < energy_meshpoints.size(); mE++){

		Reaction tmpReac = *fReaction;
		energymeshpoint_reaction.push_back(tmpReac);
		energymeshpoint_reaction[mE].SetLabEnergy(energy_meshpoints.at(mE));

		std::vector <PointCoulEx>	tmpPoinVec;
		std::vector <double>		tmpThetaVec;

		for(unsigned int mT = 0; mT < theta_meshpoints.size(); mT++){
			Nucleus nucl = *fNucleus;
      //create in place to avoid copying
      PoinMat[mE][mT] = new PointCoulEx(&nucl,&energymeshpoint_reaction[mE]);
      PoinMat[mE][mT]->SetMaxMatrix();
		  PoinMat[mE][mT]->SetUseSymmetry(fUseSymmetry);
		  PoinMat[mE][mT]->FixStep(fUseFixedStep);
		  PoinMat[mE][mT]->SetAccuracy(fAccuracy);
		  PoinMat[mE][mT]->SetProjectileExcitation(fProjectileExcitation);
			double tmpTheta = energymeshpoint_reaction[mE].ConvertThetaLabToCm(TMath::DegToRad() * theta_meshpoints.at(mT),part) * TMath::RadToDeg();
			tmpThetaVec.push_back(tmpTheta);
		}

		ThetaMat.push_back(tmpThetaVec);

	}


	std::vector<TVectorD>	tmpVector_P;
	std::vector<TVectorD>	tmpVector_CS;	
	std::vector<double>	tmpVector_cmTheta;
	

	//	If we've specified that multiple threads will be used:
	if(nThreads > 1){

		std::vector<std::thread> Threads;
		Threads.resize(nThreads);

    std::vector<std::vector<int> > energyIndices(nThreads);
    std::vector<std::vector<int> > thetaIndices(nThreads);

    int count(0);

    // prepare list of points for each thread
    // rather than creating and destroying a thread for each meshpoint, we pre-calculate which
    // points will be handled by what threads, and pass them all at once
    for (unsigned int mE = 0; mE < energy_meshpoints.size(); ++mE) {
      for (unsigned int mT = 0; mT < theta_meshpoints.size(); ++mT) {
        PoinMat[mE][mT]->SetThread(count%nThreads);
        energyIndices.at(count%nThreads).push_back(mE);
        thetaIndices.at(count%nThreads).push_back(mT);
        ++count;
      }
    }

    //assign each thread its calculations
    for (size_t t = 0; t<nThreads; ++t) {
      Threads[t] = std::thread(ThreadTaskIntegral, &PoinMat[0][0], energyIndices.at(t), thetaIndices.at(t), ThetaMat);
    }

    //wait until they are complete
    for(size_t t=0; t < nThreads; t++) {
      Threads[t].join();
    }

    //get all the results
		for(unsigned int mE = 0; mE < energy_meshpoints.size(); mE++){

			tmpVector_cmTheta.clear();
			tmpVector_CS.clear();
			tmpVector_P.clear();

			for(unsigned int mT = 0; mT < theta_meshpoints.size(); mT++){

				TVectorD tmpVec_CS;
				TVectorD tmpVec_P;
				tmpVec_P.ResizeTo(PoinMat[mE][mT]->GetProbabilitiesVector().GetNrows());
				tmpVec_CS.ResizeTo(PoinMat[mE][mT]->GetProbabilitiesVector().GetNrows());

				tmpVec_P = PoinMat[mE][mT]->GetProbabilitiesVector();

				tmpVec_CS = PoinMat[mE][mT]->GetProbabilitiesVector();
				//tmpVec_CS *= energymeshpoint_reaction.at(mE).RutherfordCM(ThetaMat.at(mE).at(mT));
		
				tmpVector_P.push_back(tmpVec_P);
				tmpVector_CS.push_back(tmpVec_CS);
				tmpVector_cmTheta.push_back(ThetaMat.at(mE).at(mT));

				index.push_back(mE);
				track_theta.push_back(theta_meshpoints.at(mT));
        point_calculations.push_back(PoinMat[mE][mT]);   //store pointer instead of value to avoid copying

			}

			cmTheta.push_back(tmpVector_cmTheta);	
			meshpointCrossSections.push_back(tmpVector_CS);	
			meshpointProbabilities.push_back(tmpVector_P);
	

		}
  }
	else{
		for(unsigned int mE = 0; mE < energy_meshpoints.size(); mE++){

			tmpVector_cmTheta.clear();
			tmpVector_CS.clear();
			tmpVector_P.clear();

			for(unsigned int mT = 0; mT < theta_meshpoints.size(); mT++){        
        PoinMat[mE][mT]->SetThread(0);
			  PoinMat[mE][mT]->CalculatePointProbabilities(ThetaMat.at(mE).at(mT));

				TVectorD tmpVec_CS;
				TVectorD tmpVec_P;
				tmpVec_P.ResizeTo( PoinMat[mE][mT]->GetProbabilitiesVector().GetNrows());
				tmpVec_CS.ResizeTo(PoinMat[mE][mT]->GetProbabilitiesVector().GetNrows());

				tmpVec_P = PoinMat[mE][mT]->GetProbabilitiesVector();

				tmpVec_CS = PoinMat[mE][mT]->GetProbabilitiesVector();
				tmpVec_CS *= energymeshpoint_reaction.at(mE).RutherfordCM(ThetaMat.at(mE).at(mT));
		
				tmpVector_P.push_back(tmpVec_P);
				tmpVector_CS.push_back(tmpVec_CS);
				tmpVector_cmTheta.push_back(ThetaMat.at(mE).at(mT));

				index.push_back(mE);
				track_theta.push_back(theta_meshpoints.at(mT));
        point_calculations.push_back(PoinMat[mE][mT]);  //store pointer instead of value to avoid copying

			}

			cmTheta.push_back(tmpVector_cmTheta);	
			meshpointCrossSections.push_back(tmpVector_CS);	
			meshpointProbabilities.push_back(tmpVector_P);
	

		}
	}

	fTensors.resize(point_calculations.size());		

	fComplete = true;
  
}

std::vector<TVectorD> Integral::GetProbabilities(){
	
	std::vector<TVectorD>	tmpProb;
	for(unsigned int i=0;i < point_calculations.size(); i++)
		tmpProb.push_back(point_calculations.at(i)->GetProbabilitiesVector());

	return tmpProb;

}

