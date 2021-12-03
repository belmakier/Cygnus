#ifndef NucleusReader_h
#define NucleusReader_h

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <regex>

#include "Nucleus.h"

///
///	\class NucleusReader
///
///	\brief Grabs nucleus data from formatted input file
///
///	This simple class is used to read nucleus information
///	from a text file and convert it into a Nucleus class. 
///

class Nucleus;

class NucleusReader
{

	public:
		NucleusReader();
		NucleusReader(const char*, bool GOSIA=false);		/*!< Construct with a filename */
		NucleusReader(const NucleusReader &n);	/*!< Copy constructor */
		NucleusReader& operator = (const NucleusReader &n); /*!< Assignment operator */
		~NucleusReader() {;}

		void		ReadNucleusFile(const char*);	/*!< Read formatted input file and create Nucleus */
		Nucleus*	GetNucleus()	{ return fNucleus; }	/*!< Return Nucleus created from input file */

		void		ReadGOSIANucleus(const char*);	/*!< Read GOSIA input file and create Nucleus */

		std::vector<int>	GetBSTInit()	const	{ return bst_i;	}
		std::vector<int>	GetBSTFinal()	const 	{ return bst_f;	}
		std::vector<int>	GetBSTLambda()	const	{ return bst_l;	}

	private:
		Nucleus 	*fNucleus;

		std::vector<int>	bst_i;
		std::vector<int>	bst_f;
		std::vector<int>	bst_l;

};
#endif
