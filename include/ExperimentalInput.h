#ifndef ExperimentalInput_h
#define ExperimentalInput_h

#include <iostream>
#include <iomanip>
#include <vector>

class ExptDoublet  {

	public:
		ExptDoublet() 							{ DetEff = 1;													}
		///
		///	Constructor containing data with symmetric uncertainties:
		///	Initial state = i
		///	Final state = f
		///	Yield = c
		///	Uncertainty = e
		///
		ExptDoublet(int i1, int f1, int i2, int f2, double c, double e)			
		{
			StateIndex_I1 = i1; 
			StateIndex_F1 = f1; 
			StateIndex_I2 = i2; 
			StateIndex_F2 = f2;
 			Counts = c; 	
			UpUncertainty =  e; 
			DnUncertainty =  e;	
			DetEff = 1;	
		}		
		~ExptDoublet()	{;}
		ExptDoublet(const ExptDoublet& e);				/*!< Copy constructor */
		ExptDoublet& operator = (const ExptDoublet& e);		/*!< Assignment operator */

		///
		///	Set data with symmetric uncertainties:
		///	Initial state = i
		///	Final state = f
		///	Yield = c
		///	Uncertainty = e
		///
		void	Set(int i1, int f1, int i2, int f2, double c, double e)			
		{ 
			StateIndex_I1 = i1; 
			StateIndex_F1 = f1; 
			StateIndex_I2 = i2; 
			StateIndex_F2 = f2;
 			Counts = c; 	
			UpUncertainty =  e; 
			DnUncertainty =  e;
		}
		void	SetCounts(double c)					{ Counts = c;											}	/*!< Set experimental yield, c */

		void	SetEfficiency(double eff)				{ DetEff = eff;											}	/*!< Define the detection efficiency. Default = 1 (i.e. efficiency corrected)  */

		int	GetInitialIndex1() 	const	{ return StateIndex_I1;	}	/*!< Return initial state index */
		int	GetFinalIndex1() 	const	{ return StateIndex_F1;	}	/*!< Return final state index */
		int	GetInitialIndex2() 	const	{ return StateIndex_I2;	}	/*!< Return initial state index */
		int	GetFinalIndex2() 	const	{ return StateIndex_F2;	}	/*!< Return final state index */
		double	GetCounts() 		const	{ return Counts;	}	/*!< Return yield */
		double	GetUpUnc() 		const	{ return UpUncertainty;	}	/*!< Return positive uncertainty */
		double	GetDnUnc() 		const	{ return DnUncertainty;	}	/*!< Return negative uncertainty */
		double	GetEfficiency()		const	{ return DetEff;	}	/*!< Return detection efficiency */

		void	Print()			const;					/*!< Print yield information */

	private:
		int 	StateIndex_I1;		/*!< Initial state index */
		int 	StateIndex_F1;           /*!< Final state index */
		int 	StateIndex_I2;		/*!< Initial state index */
		int 	StateIndex_F2;           /*!< Final state index */
		double	Counts;                 /*!< Yield */
		double	UpUncertainty;          /*!< Positive uncertainty */
		double	DnUncertainty;          /*!< Negative uncertainty */
		double	DetEff;			/*!< Detection efficiency */



};

///
///	\class ExptData
///
///	\brief Single transition information, stored within ExperimentData
///
///	The ExptData class is contains the yield information for a single 
///	transition. Each ExperimentData class object can then contain a vector 
///	of ExptData objects which includes all experimental data.	
///
///	Compatible with asymmetric uncertainties.
///

class ExptData {

	public:
		ExptData() 							{ DetEff = 1;													}
		///
		///	Constructor containing data with symmetric uncertainties:
		///	Initial state = i
		///	Final state = f
		///	Yield = c
		///	Uncertainty = e
		///
		ExptData(int i, int f, double c, double e)			{ StateIndex_I = i; StateIndex_F = f; Counts = c; UpUncertainty =  e; DnUncertainty =  e;	DetEff = 1;	}		
		///
		///	Constructor containing data with asymmetric uncertainties:
		///	Initial state = i
		///	Final state = f
		///	Yield = c
		///	Positive uncertainty = ue
		///	Negative uncertainty = de
		///
		ExptData(int i, int f, double c, double ue, double de)		{ StateIndex_I = i; StateIndex_F = f; Counts = c; UpUncertainty = ue; DnUncertainty = de;	DetEff = 1;	}
		~ExptData()	{;}
		ExptData(const ExptData& e);				/*!< Copy constructor */
		ExptData& operator = (const ExptData& e);		/*!< Assignment operator */

		///
		///	Set data with symmetric uncertainties:
		///	Initial state = i
		///	Final state = f
		///	Yield = c
		///	Uncertainty = e
		///
		void	Set(int i, int f, double c, double e)			{ Set(i,f,c,e,e);										}
		///
		///	Set data with asymmetric uncertainties:
		///	Initial state = i
		///	Final state = f
		///	Yield = c
		///	Positive uncertainty = ue
		///	Negative uncertainty = de
		///
		void	Set(int i, int f, double c, double ue, double de)	{ StateIndex_I = i; StateIndex_F = f; Counts = c; UpUncertainty = ue; DnUncertainty = de;	}
		void	SetCounts(double c)					{ Counts = c;											}	/*!< Set experimental yield, c */

		void	SetEfficiency(double eff)				{ DetEff = eff;											}	/*!< Define the detection efficiency. Default = 1 (i.e. efficiency corrected)  */

		int	GetInitialIndex() 	const	{ return StateIndex_I;	}	/*!< Return initial state index */
		int	GetFinalIndex() 	const	{ return StateIndex_F;	}	/*!< Return final state index */
		double	GetCounts() 		const	{ return Counts;	}	/*!< Return yield */
		double	GetUpUnc() 		const	{ return UpUncertainty;	}	/*!< Return positive uncertainty */
		double	GetDnUnc() 		const	{ return DnUncertainty;	}	/*!< Return negative uncertainty */
		double	GetEfficiency()		const	{ return DetEff;	}	/*!< Return detection efficiency */

		void	Print()			const;					/*!< Print yield information */

	private:
		int 	StateIndex_I;		/*!< Initial state index */
		int 	StateIndex_F;           /*!< Final state index */
		double	Counts;                 /*!< Yield */
		double	UpUncertainty;          /*!< Positive uncertainty */
		double	DnUncertainty;          /*!< Negative uncertainty */
		double	DetEff;			/*!< Detection efficiency */
	
};

///
///	\class ExperimentData
///
///	\brief Holder class containing experimental data, stored in a vector
///	of ExptData objects.
///
///	The mean theta in the center-of-mass frame is also held in this class
///	for convenience later
///
class ExperimentData{

	public:	
		ExperimentData()	{ Data.clear(); thetaCM = -1; }
		~ExperimentData()	{;}
		ExperimentData(const ExperimentData& e);	       	/*!< Copy constructor */
		ExperimentData& operator = (const ExperimentData& e);	/*!< Assignment operator */

		///
		///	Add a new ExptData value with:
		///	initial state = i
		///	final state = f
		///	yield = c
		/// 	uncertainty = e
		///
		void				AddData(int i, int f, double c, double e)	{ ExptData tmp(i,f,c,e);	Data.push_back(tmp);	}
		void				SetData(int i, ExptData exp)			{ Data.at(i) = exp;					}	/*!< Set experimental yield, index i to ExptData object exp */
		void				ClearData()					{ Data.clear();	Doublet.clear();			}	/*!< Clear experimental data */
	
		void				AddDoublet(int i1, int f1, int i2, int f2, double c, double e){
			ExptDoublet	tmp(i1,f1,i2,f2,c,e);
			Doublet.push_back(tmp);
		}

		void				SetDataEfficiency(int i, double eff)		{ Data.at(i).SetEfficiency(eff);			}	/*!<	*/
	
		ExptData			GetDataPoint(int i)	const			{ return Data.at(i);					}	/*!< Return ExptData object at index i */
		std::vector<ExptData>		GetData()		const			{ return Data;						}	/*!< Return vector of ExptData */
	
		ExptDoublet			GetDoubletPoint(int i)	const			{ return Doublet.at(i);					}	/*!< Return ExptData object at index i */
		std::vector<ExptDoublet>	GetDoublet()		const			{ return Doublet;					}	/*!< Return vector of ExptData */

		void				SetThetaCM(double t)				{ thetaCM = t;						}	/*!< Define theta CM */
		double				GetThetaCM()		const			{ return thetaCM;					}	/*!< Return theta CM */

		void				Print()			const;	/*!< Print all experimental data */

	private:
		std::vector<ExptData>		Data;		/*!< Vector of experimental yields in ExptData objects */
		std::vector<ExptDoublet>	Doublet;
		double				thetaCM;	/*!< Theta CM (mean) for this experiment */

};
#endif
