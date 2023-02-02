#include "Substate.h"

Substate::Substate(double M, int ss, int s){

	fM		= M;
	fSubstateIndex	= ss;
	fStateIndex	= s;
	fMirrorIndex	= ss;

}


Substate::Substate(const Substate& ss){

	fM		= ss.fM;
	fSubstateIndex	= ss.fSubstateIndex;
	fStateIndex	= ss.fStateIndex;
	fMirrorIndex	= ss.fMirrorIndex;
	fConnections.clear();
	fConnections.resize(ss.fConnections.size());
	for(size_t c = 0; c < fConnections.size(); c++)
		fConnections.at(c) = ss.fConnections.at(c);

}

Substate& Substate::operator = (const Substate& ss){

	fM		= ss.fM;
	fSubstateIndex	= ss.fSubstateIndex;
	fStateIndex	= ss.fStateIndex;
	fMirrorIndex	= ss.fMirrorIndex;
	fConnections.clear();
	fConnections.resize(ss.fConnections.size());
	for(size_t c = 0; c < fConnections.size(); c++)
		fConnections.at(c) = ss.fConnections.at(c);

	return *this;

}

void Substate::SetPar(Nucleus &nuc) {
  if(GetM() < 0){
    int stateindex = GetStateIndex();
    fPar = pow(-1,nuc.GetLevelP()[0]-nuc.GetLevelP()[stateindex]-nuc.GetLevelJ()[stateindex]);
  }
}

void Substate::AddLambda(const int &i, const int &l, const double &xi, const double &psi, const double &zeta) {
  fConnections.at(i).AddLambda(l, xi, psi, zeta);
}

void Substate::AddLambdaLast(const int &l, const double &xi, const double &psi, const double &zeta) {
  fConnections.at(fConnections.size()-1).AddLambda(l, xi, psi, zeta);
}
