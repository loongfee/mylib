#pragma  once
#include "Cbers04Optimize.h"

namespace Cbers04Optimize{
namespace innerOptimize{
	NEWMAT::ColumnVector getResidue(const vector< Cbers04OptStruct >& optStructList);
	void buildNormalEquation(const vector< Cbers04OptStruct >& optStructList,
		const vector<int>& parameterList,
		NEWMAT::SymmetricMatrix& A,
		NEWMAT::ColumnVector& residue,
		NEWMAT::ColumnVector& projResidue,
		double pstep_scale = 1e-6);
	void adjustment(const vector< Cbers04OptStruct >& optStructList, const vector<int>& parameterList);
}
}