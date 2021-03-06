#include "mySpectralMatching.h"
#include <arlsmat.h>
#include <arlssym.h>

using namespace std;

namespace mylib{

#define sigma_d_theta 5.0
#define sigma_d_rho 0.2
//the resolution scale of theta histogram, used when compute angle histogram of lines in the image
#define ResolutionScale  20  //10 degree
/*The following two thresholds are used to decide whether the estimated global rotation is acceptable.
*Some image pairs don't have a stable global rotation angle, e.g. the image pair of the wide baseline
*non planar scene. */
#define AcceptableAngleHistogramDifference 0.49
#define AcceptableLengthVectorDifference   0.4


/*The following four thresholds are used to decide whether a line in the left and a line in the right
*are the possible matched line pair. If they are, then their differences should be smaller than these
*thresholds.*/
#define LengthDifThreshold                 4
#define AngleDifferenceThreshold           0.7854//45degree
#define DescriptorDifThreshold             0.35//0.35, or o.5 are good values for LBD
	/*The following four threshold are used to decide whether the two possible matched line pair agree with each other.
	*They reflect the similarity of pairwise geometric information.
	*/
#define RelativeAngleDifferenceThreshold   0.7854//45degree
#define IntersectionRationDifThreshold     1
#define ProjectionRationDifThreshold       1

#define DifThreshold             10.35//0.35, or o.5 are good values for LBD
//this is used when get matching result from principal eigen vector
#define WeightOfMeanEigenVec               0.1


void mySpectralMatching::Matching(const std::vector<cvline_polar> &dataSource,const std::vector<cvline_polar> &dataReference,
		std::vector<unsigned int> &matchResult)
{
	globalRotationAngle_ = GlobalRotationOfImagePair_(dataSource, dataReference);
	BuildAdjacencyMatrix_(dataSource,dataReference);
	MatchingResultFromPrincipalEigenvector_(dataSource,dataReference,matchResult);
}

double lineLen(const cvline_polar& line)
{
	double delta_x = line.pt1.x - line.pt2.x;
	double delta_y = line.pt1.y - line.pt2.y;
	return sqrt(delta_x*delta_x + delta_y*delta_y);
}

double mySpectralMatching::GlobalRotationOfImagePair_(const std::vector<cvline_polar> &dataSource, const std::vector<cvline_polar> &dataReference)
{
	double TwoPI = 2 * M_PI;
	double rotationAngle = TwoPI;

	//step 1: compute the angle histogram of lines in the left and right images
	unsigned int dim = 360 / ResolutionScale; //number of the bins of histogram
	unsigned int index;//index in the histogram
	double direction;
	double scalar = 180 / (ResolutionScale*3.1415927);//used when compute the index
	double angleShift = (ResolutionScale*M_PI) / 360;//make sure zero is the middle of the interval

	Eigen::VectorXd angleHistLeft(dim);
	Eigen::VectorXd angleHistRight(dim);
	Eigen::VectorXd lengthLeft(dim);//lengthLeft[i] store the total line length of all the lines in the ith angle bin.
	Eigen::VectorXd lengthRight(dim);
	angleHistLeft.setZero();//.SetZero();
	angleHistRight.setZero();//.SetZero();
	lengthLeft.setZero();//.SetZero();
	lengthRight.setZero();//.SetZero();

	for (unsigned int linenum = 0; linenum<dataSource.size(); linenum++){
		direction = dataSource[linenum].theta + M_PI + angleShift;
		direction = direction<TwoPI ? direction : (direction - TwoPI);
		index = floor(direction*scalar);
		angleHistLeft[index] ++;
		lengthLeft[index] += lineLen(dataSource[linenum]);
	}
	for (unsigned int linenum = 0; linenum<dataReference.size(); linenum++){
		direction = dataReference[linenum].theta + M_PI + angleShift;
		direction = direction<TwoPI ? direction : (direction - TwoPI);
		index = floor(direction*scalar);
		angleHistRight[index] ++;
		lengthRight[index] += lineLen(dataReference[linenum]);
	}

	angleHistLeft = (1 / angleHistLeft.norm())*angleHistLeft;
	angleHistRight = (1 / angleHistRight.norm())*angleHistRight;
	lengthLeft = (1 / lengthLeft.norm())*lengthLeft;
	lengthRight = (1 / lengthRight.norm())*lengthRight;

	//  angleHistLeft.Save("histLeft.txt");
	//  angleHistRight.Save("histRight.txt");

	//step 2: find shift to decide the approximate global rotation
	Eigen::VectorXd difVec(dim);//the difference vector between left histogram and shifted right histogram
	double minDif = 10;//the minimal angle histogram difference
	double secondMinDif = 10;//the second minimal histogram difference
	unsigned int minShift;//the shift of right angle histogram when minimal difference achieved
	unsigned int secondMinShift;//the shift of right angle histogram when second minimal difference achieved

	Eigen::VectorXd lengthDifVec(dim);//the length difference vector between left and right
	double minLenDif = 10;//the minimal length difference
	double secondMinLenDif = 10;//the second minimal length difference
	unsigned int minLenShift;//the shift of right length vector when minimal length difference achieved
	unsigned int secondMinLenShift;//the shift of right length vector when the second minimal length difference achieved

	double normOfVec;
	for (unsigned int shift = 0; shift<dim; shift++){
		for (unsigned int j = 0; j<dim; j++){
			index = j + shift;
			index = index<dim ? index : (index - dim);
			difVec[j] = angleHistLeft[j] - angleHistRight[index];
			lengthDifVec[j] = lengthLeft[j] - lengthRight[index];
		}
		//find the minShift and secondMinShift for angle histogram
		normOfVec = difVec.norm();//.NormL2();
		if (normOfVec<secondMinDif){
			if (normOfVec<minDif){
				secondMinDif = minDif;
				secondMinShift = minShift;
				minDif = normOfVec;
				minShift = shift;
			}
			else{
				secondMinDif = normOfVec;
				secondMinShift = shift;
			}
		}
		//find the minLenShift and secondMinLenShift of length vector
		normOfVec = lengthDifVec.norm();//.NormL2();
		if (normOfVec<secondMinLenDif){
			if (normOfVec<minLenDif){
				secondMinLenDif = minLenDif;
				secondMinLenShift = minLenShift;
				minLenDif = normOfVec;
				minLenShift = shift;
			}
			else{
				secondMinLenDif = normOfVec;
				secondMinLenShift = shift;
			}
		}
	}

	//first check whether there exist an approximate global rotation angle between image pair
	if (minDif<AcceptableAngleHistogramDifference && minLenDif<AcceptableLengthVectorDifference){
		rotationAngle = minShift*ResolutionScale;
		if (rotationAngle>90 && 360 - rotationAngle>90){
			//In most case we believe the rotation angle between two image pairs should belong to [-Pi/2, Pi/2]
			rotationAngle = rotationAngle - 180;
		}
		rotationAngle = rotationAngle*M_PI / 180;
	}

	cout << "minimal histgram distance = " << minDif << ", Approximate global rotation angle = " << rotationAngle << endl;
	return rotationAngle;
}

fPoint mySpectralMatching::dist2(const cvline_polar& l1, const cvline_polar& l2)
{
	fPoint dist;
	// dist1
	double d1 = -l1.pt1.x * sin(l2.theta) + l1.pt1.y * cos(l2.theta) - l2.rho;
	double d2 = -l1.pt2.x * sin(l2.theta) + l1.pt2.y * cos(l2.theta) - l2.rho;
	dist.x = sqrt(d1*d1+d2*d2);
	// dist2
	d1 = -l2.pt1.x * sin(l1.theta) + l2.pt1.y * cos(l1.theta) - l1.rho;
	d2 = -l2.pt2.x * sin(l1.theta) + l2.pt2.y * cos(l1.theta) - l1.rho;
	dist.y = sqrt(d1*d1+d2*d2);
	return dist;
}

double mySpectralMatching::dist(const cvline_polar& l1, const cvline_polar& l2)
{
	return l1.theta - l2.theta;
	double d1 = -l1.pt1.x * sin(l2.theta) + l1.pt1.y * cos(l2.theta) - l2.rho;
	double d2 = -l1.pt2.x * sin(l2.theta) + l1.pt2.y * cos(l2.theta) - l2.rho;
	return sqrt(d1*d1+d2*d2);
}

void mySpectralMatching::setParameters(Eigen::VectorXd parameters)
{
	m_parameters = parameters;
}

cvline_polar mySpectralMatching::transform(const cvline_polar& l)
{
	cvline_polar outLine;
	fPoint line[2];
	line[0].x = m_parameters[0] + m_parameters[1] * l.pt1.x + m_parameters[2] * l.pt1.y;
	line[0].y = m_parameters[3] + m_parameters[4] * l.pt1.x + m_parameters[5] * l.pt1.y;
	line[1].x = m_parameters[0] + m_parameters[1] * l.pt2.x + m_parameters[2] * l.pt2.y;
	line[1].y = m_parameters[3] + m_parameters[4] * l.pt2.x + m_parameters[5] * l.pt2.y;
	outLine = line2polar(line);
	return outLine;
}

void mySpectralMatching::BuildAdjacencyMatrix_(const std::vector<cvline_polar> &dataSource,const std::vector<cvline_polar> &dataReference)
{
	double TwoPI = 2*M_PI;
	unsigned int numSource  = dataSource.size();
	unsigned int numReference = dataReference.size();
	/*first step, find nodes which are possible correspondent lines in the left and right images according to
	 *their direction, gray value  and gradient magnitude.
	 */
	nodesList_.clear();
	double angleDif;
	double lengthDif;

	Eigen::MatrixXd disMat(numSource, numReference);//store the descriptor distance of lines in left and right images.
	//std::vector<float> desLeft;
	//std::vector<float> desRight;

// //store descriptor for debug
//	Matrix<double> desCripLeft(numLineLeft,dimOfDes);
//	Matrix<double> desCripRight(numLineRight,dimOfDes);
//	for(unsigned int i=0; i<numLineLeft; i++){
//		for(unsigned int j=0; j<dimOfDes; j++){
//			desCripLeft[i][j] = linesInLeft[i].decriptor[j];
//		}
//	}
//	for(unsigned int i=0; i<numLineRight; i++){
//		for(unsigned int j=0; j<dimOfDes; j++){
//			desCripRight[i][j] = linesInRight[i].decriptor[j];
//		}
//	}
//	desCripLeft.Save("DescriptorLeft.txt");
//	desCripRight.Save("DescriptorRight.txt");

	//first compute descriptor distances

	float *desL, *desR, *desMax, *desOld;
	//std::vector<cvline_polar> transLines(numSource);
	//for (int i = 0; i < numSource; i++)
	//{
	//	transLines.push_back(transform(dataSource[i]));
	//}

	float minDis,dis,temp;
	for(int idS=0; idS<numSource; idS++){
		for(int idR=0; idR<numReference; idR++){
			fPoint d2 = dist2(dataSource[idS], dataReference[idR]);
			double d = sqrt(d2.x*d2.x+d2.y*d2.y);
			//double d = fabs(dataSource[idS].theta - dataReference[idR].theta);
			disMat(idS,idR) = d;
		}//end for(int idR=0; idR<numReference; idR++)
	}// end for(int idS=0; idS<numSource; idS++)


	for(unsigned int i=0; i<numSource; i++){
		for(unsigned int j=0; j<numReference; j++){
			if (disMat(i, j) > DifThreshold){
				continue;//the descriptor difference is too large;
			}

			Node node;//line i in left image and line j in right image pass the test, (i,j) is a possible matched line pair.
			node.srcID = i;
			node.refID = j;
			nodesList_.push_back(node);
		}//end inner loop
	}
	std::cout<<"the number of possible matched line pair = "<<nodesList_.size()<<endl;
//	desDisMat.Save("DescriptorDis.txt");

	/*Second step, build the adjacency matrix which reflect the geometric constraints between nodes.
	 *The matrix is stored in the Compressed Sparse Column(CSC) format.
	 */
	unsigned int dim = nodesList_.size();// Dimension of the problem.
	int nnz = 0;// Number of nonzero elements in adjacenceMat.
	/*adjacenceVec only store the lower part of the adjacency matrix which is a symmetric matrix.
	 *								| 0  1  0  2  0 |
	 *								| 1  0  3  0  1 |
	 *eg:  adjMatrix =		| 0  3  0  2  0 |
	 *								| 2  0  2  0  3 |
	 *								| 0  1  0  3  0 |
	 *     adjacenceVec = [0,1,0,2,0,0,3,0,1,0,2,0,0,3,0]
	 */
	//	Matrix<double> testMat(dim,dim);
	//	testMat.SetZero();
	Eigen::VectorXd adjacenceVec(dim*(dim+1)/2);
	adjacenceVec.setZero();//.SetZero();

	unsigned int idSource1, idSource2;//the id of lines in the left pair
	unsigned int idReference1, idReference2;//the id of lines in the right pair

	double similarity;
	double sigma_d_theta_3 = 3 * sigma_d_theta;
	double sigma_d_theta_square = sigma_d_theta*sigma_d_theta;
	double sigma_d_rho_3 = 3 * sigma_d_rho;
	double sigma_d_rho_square = sigma_d_rho*sigma_d_rho;
	double lammda = 0.2;
	for(unsigned int j=0; j<dim; j++){//column
		idSource1  = nodesList_[j].srcID;
		idReference1 = nodesList_[j].refID;
		for(unsigned int i=j+1; i<dim; i++){//row
			idSource2  = nodesList_[i].srcID;
			idReference2 = nodesList_[i].refID;
			if((idSource1==idSource2)||(idReference1==idReference2)){
				continue;//not satisfy the one to one match condition
			}

			//now compute the similarity score between two line pairs.
			double d_src = dist(dataSource[idSource1], dataSource[idSource2]);
			double d_ref = dist(dataReference[idReference1], dataReference[idReference2]);
			double d = fabs(d_src - d_ref);
			if (d < sigma_d_rho_3)
			{
				similarity = d*d/(2*sigma_d_rho_square);
				adjacenceVec[(2*dim-j-1)*j/2+i] = similarity;
				nnz++;
			}
			//double theta_src = dataSource[idSource1].theta - dataSource[idSource2].theta;
			//double theta_ref = dataReference[idReference1].theta - dataReference[idReference2].theta;
			//double d_theta = fabs(theta_src - theta_ref)*360.0/3.14159265;
			//double d_rho_src = dataSource[idSource1].rho - dataSource[idSource2].rho;
			//double d_rho_ref = dataReference[idReference1].rho - dataReference[idReference2].rho;
			//double d_rho = fabs(d_rho_src - d_rho_ref);
			//if (d_theta < sigma_d_theta_3)
			//{
			//	d_theta = 4.5 - d_theta*d_theta/(2*sigma_d_theta_square);
			//	if (d_rho < sigma_d_rho_3)
			//	{
			//		d_rho = 4.5 - d_rho*d_rho/(2*sigma_d_rho_square);
			//	}
			//	else
			//	{
			//		d_rho = 0.0;
			//	}
			//}
			//else
			//{
			//	d_theta = 0.0;
			//	d_rho = 0.0;
			//}
			//if (d_theta > 0.0 || d_rho > 0.0)
			//{
			//	similarity = lammda * d_theta + (1.0-lammda) * d_rho;
			//	adjacenceVec[(2*dim-j-1)*j/2+i] = similarity;
			//	nnz++;
			//}
		}
	}

	// pointer to an array that stores the nonzero elements of Adjacency matrix.
	double* adjacenceMat = new double[nnz];
	// pointer to an array that stores the row indices of the non-zeros in adjacenceMat.
	int*    irow = new int[nnz];
	// pointer to an array of pointers to the beginning of each column of adjacenceMat.
	int*    pcol = new int[dim+1];
	int idOfNNZ = 0;//the order of none zero element
	pcol[0] = 0;
	unsigned int tempValue;
	for(unsigned int j=0; j<dim; j++){//column
		for(unsigned int i=j; i<dim; i++){//row
			tempValue = (2*dim-j-1)*j/2+i;
			if(adjacenceVec[tempValue]!=0){
				adjacenceMat[idOfNNZ] = adjacenceVec[tempValue];
				irow[idOfNNZ] = i;
				idOfNNZ++;
			}
		}
		pcol[j+1] = idOfNNZ;
	}

	/*Third step, solve the principal eigenvector of the adjacency matrix using Arpack lib.
	 */
	ARluSymMatrix<double> arMatrix(dim, nnz, adjacenceMat, irow, pcol);
	ARluSymStdEig<double> dprob(2, arMatrix, "LM");// Defining what we need: the first eigenvector of arMatrix with largest magnitude.
	// Finding eigenvalues and eigenvectors.
	dprob.FindEigenvectors();
	std::cout<<"Number of 'converged' eigenvalues  : " << dprob.ConvergedEigenvalues() <<endl;
	//  cout<<"eigenvalue is = "<<dprob.Eigenvalue(0)<<", and "<<dprob.Eigenvalue(1)<<endl;
	//  if(dprob.EigenvectorsFound()){
	//  	for(unsigned int j=0; j<dim; j++){
	//  		cout<< dprob.Eigenvector(1,j) <<", ";
	//  	}
	//  	cout<<endl;
	//  }
	eigenMap_.clear();

  double meanEigenVec = 0;
  if(dprob.EigenvectorsFound()){
  	double value;
  	for(unsigned int j=0; j<dim; j++){
  		value = fabs(dprob.Eigenvector(1,j));
  		meanEigenVec += value;
  		eigenMap_.insert(std::make_pair(value,j));
  	}
  }
  minOfEigenVec_ = WeightOfMeanEigenVec*meanEigenVec/dim;
	delete[] adjacenceMat;
	delete[] irow;
	delete[] pcol;
}

void mySpectralMatching::MatchingResultFromPrincipalEigenvector_(const std::vector<cvline_polar> &dataSource,const std::vector<cvline_polar> &dataReference,
		std::vector<unsigned int > &matchResult)
{
	double TwoPI = 2*M_PI;
	std::vector<unsigned int > matchRet1;
	std::vector<unsigned int > matchRet2;
	double matchScore1 = 0;
	double matchScore2 = 0;
	EigenMAP mapCopy = eigenMap_;
	unsigned int dim = nodesList_.size();
	if (dim < 1)
	{
		return;
	}
	EigenMAP::iterator iter;
	unsigned int id,idSource2,idReference2;
	double sideValueS, sideValueR;
	double pointX,pointY;


  ////store eigenMap for debug
  //std::fstream resMap;
  //ostringstream fileNameMap;
  //fileNameMap<<"eigenVec.txt";
  //resMap.open(fileNameMap.str().c_str(), std::ios::out);

  //Eigen::MatrixXd mat(dataSource.size(),dataReference.size());
  //mat.setZero();//.SetZero();
  //for(iter = eigenMap_.begin();iter!=eigenMap_.end(); iter++){
  //	id = iter->second;
  //	resMap<<nodesList_[id].srcID<<"    "<<nodesList_[id].refID<<"   "<<iter->first<<endl;
  //	mat(nodesList_[id].srcID,nodesList_[id].refID) = iter->first;
  //}
  ////mat.Save("eigenMap.txt");
  //resMap.flush();
  //resMap.close();


	/*first try, start from the top element in eigenmap */
	while(1){
		iter = eigenMap_.begin();
		//if the top element in the map has small value, then there is no need to continue find more matching line pairs;
		if(iter->first < minOfEigenVec_){
			break;
		}
		id = iter->second;
		unsigned int idSource1 = nodesList_[id].srcID;
		unsigned int idReference1= nodesList_[id].refID;
		matchRet1.push_back(idSource1);
		matchRet1.push_back(idReference1);
		matchScore1 += iter->first;
		eigenMap_.erase(iter++);
		//remove all potential assignments in conflict with top matched line pair
		for( ; iter->first >= minOfEigenVec_; ){
			id = iter->second;
			idSource2 = nodesList_[id].srcID;
			idReference2= nodesList_[id].refID;
			//check one to one match condition
			if((idSource1==idSource2)||(idReference1==idReference2)){
				eigenMap_.erase(iter++);
				continue;//not satisfy the one to one match condition
			}
			iter++;
		}
	}//end while(stillLoop)
	matchResult = matchRet1;
	std::cout<<"matchRet1.size"<<matchRet1.size()<<", minOfEigenVec_= "<<minOfEigenVec_<<endl;
}
}