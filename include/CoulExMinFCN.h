#ifndef CoulExMinFCN_h
#define CoulExMinFCN_h

#include "ScalingParameter.h"
#include "Nucleus.h"
#include "Minuit2/FCNBase.h"
#include "Literature.h"
#include "ExperimentalInput.h"
#include "MatrixElement.h"
#include "TransitionRates.h"
#include "PointCoulEx.h"
#include "Reaction.h"
#include "GammaYield.h"

#include <ctime>
#include <thread>
#include <iomanip>
#include <vector>
#include <chrono>

class MatrixElements;
class ExperimentData;
class LitLifetime;
class LitBranchingRatio;
class LitMixingRatio;

///
///	\class CoulExMinFCN
///
///	\brief Contains the definition of the chi-squared function used in the minimization
///	routines
///


class CoulExMinFCN { // : public ROOT::Minuit2::FCNBase{

	public : 

		CoulExMinFCN(std::vector<ExperimentData> d) : exptData(d)		{ 
												verbose = false;  
												iter = 0;	
												nThreads = 1;
											}	/*!< Construct object with vector of experimental data to be fit */
		virtual ~CoulExMinFCN()							{;				}

		void	SetupCalculation();	/*!< Prepare the calculation */

		virtual void ClearAll();	/*!< Clear all vectors */

		double up() const	 						{ return theErrorDef;			}	/*!< Required by ROOT::Minimizer */
		double operator()(const double*);											/*!< Required by ROOT::Minimizer */
		void	setErrorDef(double def)						{ theErrorDef = def;			}	/*!< Required by ROOT::Minimizer */

		void	SetMatrixElements(std::vector<MatrixElement> m)			{ ME = m;				}	/*!< Define the vector of MatrixElement objects to be fitted */
		std::vector<MatrixElement>	GetMatrixElements() 		const	{ return ME;				}	/*!< Return the vector of MatrixElement objects to be fitted */

		std::vector<ScalingParameter>	GetScalingParameters()			{ return scalingParameters;		}	/*!< Return the vector of ScalingParameter objects for fitting */
		void	SetScalingParameters(std::vector<ScalingParameter> s)		{ scalingParameters = s;		}	/*!< Define the vector of ScalingParameter objects for fitting */
		void	AddScalingParameter(ScalingParameter s)				{ scalingParameters.push_back(s);	}	/*!< Append a new ScalingParameter to the vector */
		void	ClearScalingParameters()					{ scalingParameters.clear();		}	/*!< Clear the vector of ScalingParameter objects */

		void	SetPointCalcs(std::vector<PointCoulEx> p)			{ pointCalcs = p;			}	/*!< Define the vector of PointCoulEx objects (one for each experiment) */
		std::vector<PointCoulEx>	GetPointCalcs() 		const	{ return pointCalcs;			}	/*!< Return the vector of PointCoulEx objects (one for each experiment) */

		void	SetData(std::vector<ExperimentData> d)				{ exptData = d;				}	/*!< Define the vector of ExperimentData objects (one for each experiment) */
		std::vector<ExperimentData>	GetData() 			const	{ return exptData;			}	/*!< Return the vector of ExperimentData objects (one for each experiment) */

		void	SetLitLifetimes(std::vector<LitLifetime> l)			{ litLifetimes = l;			}	/*!< Define the vector of LitLifetime objects defining the literature lifetime data for fitting */
		std::vector<LitLifetime>	GetLitLifetimes() 		const	{ return litLifetimes;			}	/*!< Return the vector of LitLifetime objects defining the literature lifetime data for fitting */

		void	SetLitBranching(std::vector<LitBranchingRatio> b) 		{ litBranchingRatios = b;		}	/*!< Define the vector of LitBranchingRatio objects defining the literature branching ratio data for fitting */
		std::vector<LitBranchingRatio>	GetLitBranching() 		const	{ return litBranchingRatios;		}	/*!< Return the vector of LitBranchingRatio objects defining the literature branching ratio data for fitting */
	
		void	SetLitMixing(std::vector<LitMixingRatio> m)			{ litMixingRatios = m;			}	/*!< Define the vector of LitMixingRatio objects defining the literature mixing ratio data for fitting */
		std::vector<LitMixingRatio>	GetLitMixing() 			const	{ return litMixingRatios;		}	/*!< Return the vector of LitMixingRatio objects defining the literature mixing ratio data for fitting */

		std::vector<TMatrixD>		GetEffectiveCrossSection() 	const	{ return EffectiveCrossSection;		}	/*!< Return the "effective cross section" = direct population + feeding */

		void	SetBaseNucleus(Nucleus* nucl)					{ fNucleus_Base = *nucl;		}	/*!< Define the base nucleus (not to be varied in fitting) */
		Nucleus				GetBaseNucleus() 		const	{ return fNucleus_Base;			}	/*!< Return the base nucleus (not to be varied in fitting) */
	
		void	SetNucleus(Nucleus *nucl)					{ fNucleus = *nucl;			}	/*!< Define the fitting nucleus (varied in fitting) */
		Nucleus				GetNucleus() 			const	{ return fNucleus;			}	/*!< Return the fitting nucleus (varied in fitting) */

		void	SetNpar(int n)							{ nPar = n;				}	/*!< Define the number of fitting parameters (fitting matrix elements + scaling parameters) */
		int	GetNpar()						const	{ return nPar;				}	/*!< Return the number of fitting parameters (fitting matrix elements + scaling parameters) */

		void	SetCorrectionFactors(std::vector<TVectorD> v)			{ correctionFactors = v;		}	/*!< Define the vector of correction factors between point and integrated cross sections */
		std::vector<TVectorD> GetCorrectionFactors()			const	{ return correctionFactors;		}	/*!< Return the vector of correciton factors between point and integrated cross sections */

		void	SetVerbose(bool b = true)					{ verbose = b;				}	/*!< Define the verbocity of the minimization */	
		bool	GetVerbose() 						const	{ return verbose;			}	/*!< Return the verbocity of the minimization */

		void	SetIter(int i)							{ nIterations = i;			}	/*!< Define the number of iterations (MINUIT2) */
		void 	SetCalls(int i)							{ nCalls = i;				}	/*!< Define the number of function calls (GSL) */

		void	SetNthreads(int n)						{ nThreads = n;				}	/*!< Define the number of cores the function is allowed to use */
		int	GetNthreads()						const	{ return nThreads;			}	/*!< Return the number of cores the function is allowed to use */

		double	GetParameter(int i)					const	{ return parameters.at(i);		}	/*!< Return the fitting parameter indexed by i */
		std::vector<double> GetParameters()					{ return parameters;			}	/*!< Return the vector of fitting parameters */

		void	ResetIter()							{ iter = 0;				}	/*!< Reset the iteration number */

	private :

		std::vector<double>		parameters;			// Matrix elements + scaling factors
		std::vector<MatrixElement>	ME;				// Preset to related parameters to matrix elements
		std::vector<ScalingParameter>	scalingParameters;		// Scaling parameters

		std::vector<TVectorD>		correctionFactors;

		std::vector<PointCoulEx>	pointCalcs;			// Point calculations
		std::vector<ExperimentData>	exptData;			// Experimental data (one vector entry for each data subset)
		std::vector<LitLifetime>	litLifetimes;			// Literature data, lifetimes
		std::vector<LitBranchingRatio>	litBranchingRatios;		// Literature data, branching ratios
		std::vector<LitMixingRatio>	litMixingRatios;		// Literature data, mixing ratios
		double				theErrorDef;

		std::vector<TMatrixD>		EffectiveCrossSection;

		Nucleus				fNucleus;
		Nucleus				fNucleus_Base;

		int				nPar;

		TransitionRates			fRates;

		bool				verbose;

		int				nCalls;
		int				nIterations;
		int				iter;

		int				nThreads;

		std::vector<int>		exptIndex;

};
#endif