#pragma  once
#include "Hj1Optimize.h"

namespace Hj1Optimize{
namespace innerOptimize{
	NEWMAT::ColumnVector getResidue(const vector< HJ1OptStruct >& hj1OptStructList);
	void buildNormalEquation(const vector< HJ1OptStruct >& hj1OptStructList, 
		const vector<int>& parameterList,
		NEWMAT::SymmetricMatrix& A,
		NEWMAT::ColumnVector& residue,
		NEWMAT::ColumnVector& projResidue,
		double pstep_scale = 1e-6);
	void adjustment(const vector< HJ1OptStruct >& hj1OptStructList, const vector<int>& parameterList);
}
}