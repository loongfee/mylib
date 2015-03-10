#include "ecmlr.h"

ecmlr::ecmlr( const vector<cvline_polar>& source_lines_polar, const vector<cvline_polar>& reference_lines_polar)
{

	m_source_lines_polar = source_lines_polar;
	//for (int i = 0;i < (int)m_model_lines_polar.size();++i)
	//{
	//	normalize_segment(m_model_lines_polar[i]);
	//}
	m_reference_lines_polar = reference_lines_polar;
	//for (int i = 0;i < (int)m_observe_lines_polar.size();++i)
	//{
	//	normalize_segment(m_observe_lines_polar[i]);
	//}

	m_nX = (int)m_source_lines_polar.size();
	m_nY = (int)m_reference_lines_polar.size();
	m_ndim = 2;
	//initAdjustableParameters();
}

bool ecmlr::solve()
{
	int nX = m_nX;
	int nY = m_nY;
	int ndim = m_ndim;
	initialization();

	double convergence_epsilon = 2.0e-4;
	double delta;
	int max_iteration = 200;
	int iTime = 0;

	string foldName = m_sourceImageFile.substr(0, m_sourceImageFile.find_last_of('\\'));

	Eigen::VectorXd oldParameters = getParameters();
	do 
	{
		// compute transformed straight lines
		updateTransLines();
		
		double tmp = m_delta_2;
		updateDeltaSquare();
		cout<<"iteration:\t"<<iTime+1<<endl;
		cout<<"rms:\t"<<m_delta_2<<endl;
		cout<<oldParameters<<endl;
		if (iTime > 0 && m_delta_2 > tmp)
		{
			cout<<"Warning: the model is getting worse!"<<endl;
			break;
		}
		if (m_delta_2 < 1.0)
		{
			cout<<"rms epsilon reached!"<<endl;
			break;
		}
		if (fabs(m_delta_2 - tmp) < 0.01)
		{
			cout<<"delta rms epsilon reached!"<<endl;
			break;
		}

		// E-step
		E_Step();

		// CM-step A
		// estimate the trans parameters
		parameters_optimization();
		Eigen::VectorXd tmpParameters = getParameters();
		delta = (oldParameters - tmpParameters).norm();
		oldParameters = tmpParameters;

		// CM-step B
		// estimate the covariance matrices
		//updateCovariances();
		//m_delta_2 *= 0.9;
	} while (++iTime < max_iteration && fabs(delta) > convergence_epsilon);
	if(iTime == max_iteration)
	{
		cout<<"Warning: the max number of iteration reached, and the results may be not accurate."<<endl;
	}
	classification();
	int nRemoved = 0;

	double delta_ = sqrt(m_delta_2);
	for(int j = 1; j < nY;j++)
	{
		if(m_z[j] == nX) continue;
		// remove outliers
		fPoint dist = segment2polarline(m_trans_lines_polar[m_z[j]], m_reference_lines_polar[j]);
		double dist2 = sqrt(dist.x*dist.x + dist.y*dist.y);
		cvline_polar l_obs = m_reference_lines_polar[j];
		cvline_polar l_trans = m_trans_lines_polar[m_z[j]];
		fPoint center_obs((l_obs.pt1.x+l_obs.pt2.x)*0.5, (l_obs.pt1.y+l_obs.pt2.y)*0.5);
		fPoint center_trans((l_trans.pt1.x+l_trans.pt2.x)*0.5, (l_trans.pt1.y+l_trans.pt2.y)*0.5);
		double cen_dis = (center_obs.x-center_trans.x)*(center_obs.x-center_trans.x) + (center_obs.y-center_trans.y)*(center_obs.y-center_trans.y);
		if (cen_dis > 1000.0)
		{
			nRemoved++;
			continue;
		}
		if (dist2 > 2*delta_)
		{	
			nRemoved++;
			continue;
		}
		m_inliers.push_back(j);
		m_inliers.push_back(m_z[j]);
	}
	int nTotal = (int)m_inliers.size() / 2;
	vector<int> inliers(nTotal);
	//mat dataPoints(nTotal, 6);
	for(int i = 0;i < nTotal;++i)
	{
		inliers[i] = i;
	}

	////output precision report
	//FILE* pf;
	//pf = fopen((foldName+"\\report.txt").c_str(), "w+");
	//fprintf(pf, "%4s%15s%15s%15s%15s%15s%15s%15s\n", "ID", "cen_dis", "e1", "e2", "x1", "y1", "x2", "y2");
	////fReport.open("pic\\report.txt", std::ios_base::out);
	////fReport<<"ID\t e1\t e1\t x1\t y1\t x2\t y2\n";
	////int n = (int)m_inliers.size() / 2;
	int n = (int)inliers.size();
	//for (int p = 0;p < n;++p)
	//{
	//	int j = inliers[p];
	//	fPoint dist = segment2polarline(m_trans_lines_polar[m_inliers[2*j+1]], m_reference_lines_polar[m_inliers[2*j]]);
	//	//fReport<<j+1
	//	//	<<"\t"<<dist.x
	//	//	<<"\t"<<dist.y
	//	//	<<"\t"<<m_observe_lines_polar[m_inliers[j]].pt1.x
	//	//	<<"\t"<<m_observe_lines_polar[m_inliers[j]].pt1.y
	//	//	<<"\t"<<m_observe_lines_polar[m_inliers[j]].pt2.x
	//	//	<<"\t"<<m_observe_lines_polar[m_inliers[j]].pt2.y
	//	//	<<endl;
	//	cvline_polar l_obs = m_reference_lines_polar[m_inliers[2*j]];
	//	cvline_polar l_trans = m_trans_lines_polar[m_inliers[2*j+1]];
	//	fPoint center_obs((l_obs.pt1.x+l_obs.pt2.x)*0.5, (l_obs.pt1.y+l_obs.pt2.y)*0.5);
	//	fPoint center_trans((l_trans.pt1.x+l_trans.pt2.x)*0.5, (l_trans.pt1.y+l_trans.pt2.y)*0.5);
	//	double cen_dis = (center_obs.x-center_trans.x)*(center_obs.x-center_trans.x) + (center_obs.y-center_trans.y)*(center_obs.y-center_trans.y);
	//	fprintf(pf, "%4d%15f%15f%15f%15f%15f%15f%15f\n", j+1, cen_dis, dist.x, dist.y, m_reference_lines_polar[m_inliers[2*j]].pt1.x,
	//		m_reference_lines_polar[m_inliers[2*j]].pt1.y, m_reference_lines_polar[m_inliers[2*j]].pt2.x, m_reference_lines_polar[m_inliers[2*j]].pt2.y);
	//}
	//fclose(pf);

	//// control line segments
	//pf = fopen((foldName+"\\lines.txt").c_str(), "w+");
	//for (int p = 0;p < n;++p)
	//{
	//	int j = inliers[p];
	//	fprintf(pf, "%4d%15f%15f%15f%15f%15f\n", j+1, m_source_lines_polar[m_inliers[2*j+1]].pt1.x, m_source_lines_polar[m_inliers[2*j+1]].pt1.y, 
	//		m_reference_lines_polar[m_inliers[2*j]].pt1.x, m_reference_lines_polar[m_inliers[2*j]].pt1.y);
	//	fprintf(pf, "%4d%15f%15f%15f%15f%15f\n", j+1, m_source_lines_polar[m_inliers[2*j+1]].pt2.x, m_source_lines_polar[m_inliers[2*j+1]].pt2.y,
	//		m_reference_lines_polar[m_inliers[2*j]].pt2.x, m_reference_lines_polar[m_inliers[2*j]].pt2.y);
	//}
	//fclose(pf);
	////fReport.close();


	// draw mathes
	//vector<cvline_polar> model_ls;
	//vector<cvline_polar> observe_ls;
	m_matched_src_ls.clear();
	m_matched_ref_ls.clear();
	//int n = (int)m_inliers.size() / 2;
	for (int p = 0;p < n;++p)
	{
		int j = inliers[p];
		m_matched_src_ls.push_back(m_source_lines_polar[m_inliers[2*j+1]]);
		m_matched_ref_ls.push_back(m_reference_lines_polar[m_inliers[2*j]]);
	}

	//cv::Mat img_observe = cv::imread(m_referenceImageFile);
	//cv::Mat img_model = cv::imread(m_sourceImageFile);
	//cv::Mat img_matches;
	//drawLineMatches(img_model, m_matched_ref_ls, img_observe, m_matched_src_ls,
	//	img_matches, cv::Scalar::all(-1), cv::Scalar::all(-1),
	//	vector<vector<char> >(), cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS, 2);

	//cv::imwrite((foldName+"\\matches.png").c_str(), img_matches);

	vector<double> param(get_PARAMETER_NUM());
	for (int i = 0;i < get_PARAMETER_NUM();++i) param[i] = oldParameters[i];
	//drawLastResult(m_observeImageFile, m_modelImageFile, param, m_observe_lines_polar, m_trans_lines_polar, m_inliers, "pic\\result1.jpg", "pic\\result2.jpg");
	//drawLastResult(m_sourceImageFile, m_referenceImageFile, m_reference_lines_polar, m_source_lines_polar, m_inliers, "pic\\result1.jpg", "pic\\result2.jpg");
	//drawLastResult(m_observe_lines_polar, m_trans_lines_polar, m_inliers, "pic\\result1.jpg", "pic\\result2.jpg");
	cout<<"EM based point set registration finished after "<<iTime<<" times of iterations."<<endl;
	cout<<"********************************************************************************"<<endl<<endl;
	return true;
}


bool ecmlr::classification()
{
	int nX = m_nX;
	int nY = m_nY;
	int ndim = m_ndim;

	// classification
	for (int j = 0;j < nY;++j)
	{
		int maxIndex = 0;
		double maxValue = 0.0;
		for (int i = 0;i < nX + 1;++i)
		{
			if(m_alpha(j, i) > maxValue)
			{
				maxValue = m_alpha(j, i);
				maxIndex = i;
			}
		}
		m_z[j] = maxIndex;
	}
	return true;
}

double ecmlr::acceptance_ratio(const correspondence_struct& c, const correspondence_struct& c_prime)
{
	double result = 1.0;
	for (int j = 0;j < (int)c.forward_map.size();++j)
	{
		double test = m_alpha(34, 91);
		double a_prime = m_alpha(j, c_prime.forward_map[j]);
		double a = m_alpha(j, c.forward_map[j]);
		result *= (m_alpha(j, c_prime.forward_map[j]))/(m_alpha(j, c.forward_map[j]) + DBL_EPSILON);
	}
	return result;
}

correspondence_struct ecmlr::getValidCorrespondence()
{
	int nX = m_nX;
	int nY = m_nY;
	int ndim = m_ndim;

	// valid correspondence vector
	// indices start from 0
	// nX indicates invisible or outlier
	correspondence_struct correspond;
	correspond.forward_map.resize(nY, nX);
	correspond.inverse_map.resize(nX, -1);

	// randomly choose the start index
	// initialize random seed:
	srand ( time(NULL) );
	// generate start index:
	int start_pos = rand() % nY;

	vector<int> candidate_list(nX+1);
	for (int i = 0;i < nX+1;++i) candidate_list[i] = i;

	for (int j = 0;j < nY;++j)
	{
		int real_pose = (j + start_pos) % nY;
		int maxIndex = 0;
		double maxValue = 0.0;
		if(candidate_list.size() < 2) break;
		for (int i = 0;i < (int)candidate_list.size();++i)
		{
			if(m_alpha(real_pose, candidate_list[i]) > maxValue)
			{
				maxValue = m_alpha(real_pose, candidate_list[i]);
				maxIndex = i;
			}
		}
		correspond.forward_map[real_pose] = candidate_list[maxIndex];

		// if it is not outlier
		if(candidate_list[maxIndex] != nX)
		{
			correspond.inverse_map[candidate_list[maxIndex]] = real_pose;
			candidate_list.erase(candidate_list.begin()+maxIndex);
		}
	}
	return correspond;
}

correspondence_struct ecmlr::assignment_multi()
{
	int nX = m_nX;
	int nY = m_nY;
	int ndim = m_ndim;

	// valid correspondence vector
	// indices start from 0
	// nX indicates invisible or outlier
	correspondence_struct correspond;
	correspond.forward_map.resize(nY, nX);
	correspond.inverse_map.resize(nX, -1);

	vector<int> openList(nY);
	for (int i = 0;i < nY;++i) openList[i] = i;

	// randomly choose the start index
	// initialize random seed:
	srand ( time(NULL) );
	// generate start index:
	//int j = rand() % openList.size();
	int j = 0;
	vector<int> closeList;
	while((int)openList.size() > 0)
	{
		int real_pose = j % (int)openList.size();
		//if (27 == openList[real_pose])
		//{
		//	cout<<"debug.."<<endl;
		//}

		int maxIndex = 0;
		double maxValue = 0.0;
		for (int i = 0;i < nX+1;++i)
		{
			if(m_alpha(openList[real_pose], i) > maxValue)
			{
				maxValue = m_alpha(openList[real_pose], i);
				maxIndex = i;
			}
		}

		if(maxIndex == nX)
		{
			closeList.push_back(openList[real_pose]);
			openList.erase(openList.begin()+real_pose);
			continue;
		}

		correspond.forward_map[openList[real_pose]] = maxIndex;
		correspond.inverse_map[maxIndex] = openList[real_pose];
		closeList.push_back(openList[real_pose]);
		openList.erase(openList.begin()+real_pose);
	}

	return correspond;
}

correspondence_struct ecmlr::assignment_one()
{
	int nX = m_nX;
	int nY = m_nY;
	int ndim = m_ndim;

	// valid correspondence vector
	// indices start from 0
	// nX indicates invisible or outlier
	correspondence_struct correspond;
	correspond.forward_map.resize(nY, nX);
	correspond.inverse_map.resize(nX, -1);

	vector<int> openList(nY);
	for (int i = 0;i < nY;++i) openList[i] = i;

	// randomly choose the start index
	// initialize random seed:
	srand ( time(NULL) );
	// generate start index:
	//int j = rand() % openList.size();
	int j = 0;
	vector<int> closeList;
	while((int)openList.size() > 0)
	{
		int real_pose = j % (int)openList.size();
		//if (27 == openList[real_pose])
		//{
		//	cout<<"debug.."<<endl;
		//}


		int maxIndex = nX;
		double maxValue = 0.0;
		for (int i = 0;i < nX+1;++i)
		{
			if(m_alpha(openList[real_pose], i) > maxValue)
			{
				maxValue = m_alpha(openList[real_pose], i);
				maxIndex = i;
			}
		}
		//if (46 == maxIndex)
		//{
		//	cout<<"debug.."<<endl;
		//}

		if(maxIndex == nX)
		{
			closeList.push_back(openList[real_pose]);
			openList.erase(openList.begin()+real_pose);
			continue;
		}
		// 
		if(correspond.inverse_map[maxIndex] != -1)
		{
			// assigned
			int old_observation = correspond.inverse_map[maxIndex];
			cvline_polar old_line = m_reference_lines_polar[old_observation];
			cvline_polar new_line = m_reference_lines_polar[openList[real_pose]];
			cvline_polar trans_line = m_trans_lines_polar[maxIndex];
			fPoint center_old((old_line.pt1.x+old_line.pt2.x)*0.5, (old_line.pt1.y+old_line.pt2.y)*0.5);
			fPoint center_new((new_line.pt1.x+new_line.pt2.x)*0.5, (new_line.pt1.y+new_line.pt2.y)*0.5);
			fPoint center_trans((trans_line.pt1.x+trans_line.pt2.x)*0.5, (trans_line.pt1.y+trans_line.pt2.y)*0.5);
			double dis_old = (center_old.x-center_trans.x)*(center_old.x-center_trans.x) + (center_old.y-center_trans.y)*(center_old.y-center_trans.y);
			double dis_new = (center_new.x-center_trans.x)*(center_new.x-center_trans.x) + (center_new.y-center_trans.y)*(center_new.y-center_trans.y);
			if(dis_new < dis_old)
			{
				correspond.forward_map[openList[real_pose]] = maxIndex;
				correspond.inverse_map[maxIndex] = openList[real_pose];
				closeList.push_back(openList[real_pose]);
				openList.erase(openList.begin()+real_pose);

				correspond.forward_map[old_observation] = nX;
				openList.push_back(old_observation);
				//j++;
			}
			else
			{
				int maxIndex = nX;
				double maxValue = 0.0;
				for (int i = 0;i < nX+1;++i)
				{
					if(correspond.inverse_map[i] == -1 && m_alpha(openList[real_pose], i) > maxValue)
					{
						maxValue = m_alpha(openList[real_pose], i);
						maxIndex = i;
					}
				}
				if(maxIndex == nX)
				{
					closeList.push_back(openList[real_pose]);
					openList.erase(openList.begin()+real_pose);
					continue;
				}

				correspond.forward_map[openList[real_pose]] = maxIndex;
				correspond.inverse_map[maxIndex] = openList[real_pose];
				closeList.push_back(openList[real_pose]);
				openList.erase(openList.begin()+real_pose);
			}
		}
		else
		{
			// not assigned
			correspond.forward_map[openList[real_pose]] = maxIndex;
			correspond.inverse_map[maxIndex] = openList[real_pose];
			closeList.push_back(openList[real_pose]);
			openList.erase(openList.begin()+real_pose);
		}
	}

	return correspond;
}

bool ecmlr::classification1()
{
	int nX = m_nX;
	int nY = m_nY;
	int ndim = m_ndim;

	int nSample = 1000;
	int nDiscard = 50;
	correspondence_struct current_correspondence = getValidCorrespondence();

	Eigen::MatrixXi accumulateMatrix(nY, nX + 1);
	accumulateMatrix.fill(0);
	int iCount = 0;
	do 
	{
		//double f1 = getCorrepondenceProbility(current_correspondence);
		correspondence_struct new_correspondence = MarkovNext(current_correspondence);
		//double f2 = getCorrepondenceProbility(new_correspondence);
		double acceptance = acceptance_ratio(current_correspondence, new_correspondence);
		if(acceptance < 0.999)
		{
			// reject
			continue;
		}
		else
		{
			// accept
			current_correspondence = new_correspondence;
			iCount++;
		}
	} while (iCount < nSample);

	m_z = current_correspondence.forward_map;

	return true;
}

bool ecmlr::MetropolisHastings()
{
	int nX = m_nX;
	int nY = m_nY;
	int ndim = m_ndim;

	int nSample = 10000;
	int nDiscard = 50;
	correspondence_struct current_correspondence = getValidCorrespondence();

	Eigen::MatrixXi accumulateMatrix(nY, nX + 1);
	accumulateMatrix.fill(0);
	int iCount = 0;
	do 
	{
		//double f1 = getCorrepondenceProbility(current_correspondence);
		correspondence_struct new_correspondence = MarkovNext(current_correspondence);
		//double f2 = getCorrepondenceProbility(new_correspondence);
		double acceptance = acceptance_ratio(current_correspondence, new_correspondence);
		if(acceptance < 0.999)
		{
			// reject
			continue;
		}
		else
		{
			// accept
			current_correspondence = new_correspondence;
			if (iCount >= nDiscard)
			{
				for (int k = 0;k < (int)current_correspondence.forward_map.size();++k)
				{
					accumulateMatrix(k, current_correspondence.forward_map[k]) += 1;
				}
			}
			iCount++;
		}
	} while (iCount < nSample);

	// update the corresponding weights
	int num = nSample - nDiscard;
	for (int j = 0;j < nY;++j)
	{
		for(int i = 0;i < nX;++i)
		{
			m_alpha(j, i) = accumulateMatrix(j, i) / (double)num;
		}
	}

	//m_alpha.fill(0.0);
	//for (int j = 0;j < nY;++j)
	//{
	//	m_alpha(j, current_correspondence.forward_map[j]) = 1.0;
	//}

	return true;
}

double ecmlr::getCorrepondenceProbility(const correspondence_struct& correspond)
{
	double result = 1.0;
	for (int j = 0;j < (int)correspond.forward_map.size();++j)
	{
		result *= m_alpha(j, correspond.forward_map[j]);
	}
	return result;
}

correspondence_struct ecmlr::MarkovNext(const correspondence_struct& correspond)
{
	int nX = m_nX;
	int nY = m_nY;
	int ndim = m_ndim;

	correspondence_struct newCorrespond = correspond;

	// randomly choose the start index
	// initialize random seed:
	srand ( time(NULL) );
	// generate start index:
	int oldSource = rand() % nY;

	vector<int> openList(nX + 1);
	for (int i = 0;i < nX+1;++i) openList[i] = i;
	do 
	{
		int maxIndex = 0;
		double maxValue = 0.0;
		for (int i = 0;i < (int)openList.size();++i)
		{
			if(m_alpha(oldSource, openList[i]) > maxValue)
			{
				maxValue = m_alpha(oldSource, openList[i]);
				maxIndex = i;
			}
		}

		int oldTarget = newCorrespond.forward_map[oldSource];
		int newTarget = openList[maxIndex];
		if(openList[maxIndex] == nX || newCorrespond.inverse_map[newTarget] == -1)
		{
			// outlier
			if (oldTarget != nX) newCorrespond.inverse_map[oldTarget] = -1;
			newCorrespond.forward_map[oldSource] = newTarget;
			break;
		}
		else
		{
			int newSource = newCorrespond.inverse_map[newTarget];
			if (oldTarget != nX) newCorrespond.inverse_map[oldTarget] = -1;
			newCorrespond.inverse_map[newTarget] = oldSource;
			newCorrespond.forward_map[oldSource] = newTarget;
			oldSource = newSource;

			openList.erase(openList.begin() + maxIndex);
		}
	} while ((int)openList.size() > 0);

	return newCorrespond;
}


void ecmlr::updateDeltaSquare()
{
	int nX = m_nX;
	int nY = m_nY;
	int ndim = m_ndim;

	double sum = 0.0;
	for(int i = 0;i < nX;++i)
	{
		for (int j = 0;j < nY;++j)
		{
			//fPoint dist = distance2(m_trans_lines_polar[i], m_observe_lines_polar[j]);
			fPoint dist = segment2polarline(m_trans_lines_polar[i], m_reference_lines_polar[j]);
			//fPoint dist;
			//dist.x = point2polarline(m_trans_lines_polar[i].pt1, m_source_lines_polar[j]);
			//dist.y = point2polarline(m_trans_lines_polar[i].pt2, m_source_lines_polar[j]);
			sum += (dist.x*dist.x + dist.y*dist.y)*m_alpha(j,i);
		}
	}

	//m_delta_2 = sum / (nY * ndim);
	//m_delta_2 = sum / nY;
	m_delta_2 = sum / (m_lammda_sum * ndim);
	if(m_delta_2 < 1.0) m_delta_2 = 1.0;
}

bool ecmlr::initialization()
{
	//set the size of posterior matrix
	m_alpha.resize(m_nY, m_nX + 1);
	double w = 1.0 / (m_nX + 1);
	m_alpha.fill(w);

	// set the size of virtual observation matrix
	m_w.resize(m_nX);

	// set the weight vector of virtual observations
	m_lambda.resize(m_nX);

	m_lammda_sum = 0.0;
	for (int i = 0;i < m_nX;++i)
	{
		m_lambda[i] = 0;
		double theta = 0.0;
		double rho = 0.0;
		for (int j = 0;j < m_nY;++j)
		{
			m_lambda[i] += m_alpha(j, i);
			theta += m_alpha(j, i) * m_reference_lines_polar[j].theta;
			rho += m_alpha(j, i) * m_reference_lines_polar[j].rho;
		}
		m_lammda_sum += m_lambda[i];
		m_w[i].theta = theta;
		m_w[i].rho = rho;
	}

	//
	m_z.resize(m_nY);
	//m_z.set_size(nY);

	//
	m_trans_lines_polar.resize(m_nX);
	//m_transData.set_size(nX, ndim);

	return true;
}

bool ecmlr::setMatched()
{
	initAdjustableParameters();
	if (m_nY != m_nX)
	{
		return false;
	}

	//set the size of posterior matrix
	m_alpha.resize(m_nY, m_nX + 1);
	double w = 1.0 / (m_nX + 1);
	m_alpha.fill(0.0);
	for (int i = 0;i < m_nY;++i)
	{
		m_alpha(i,i) = 1.0;
	}

	// set the size of virtual observation matrix
	m_w.resize(m_nX);

	// set the weight vector of virtual observations
	m_lambda.resize(m_nX);

	m_lammda_sum = 0.0;
	for (int i = 0;i < m_nX;++i)
	{
		m_lambda[i] = 0;
		double theta = 0.0;
		double rho = 0.0;
		for (int j = 0;j < m_nY;++j)
		{
			m_lambda[i] += m_alpha(j, i);
			theta += m_alpha(j, i) * m_reference_lines_polar[j].theta;
			rho += m_alpha(j, i) * m_reference_lines_polar[j].rho;
		}
		m_lammda_sum += m_lambda[i];
		m_w[i].theta = theta;
		m_w[i].rho = rho;
	}

	//
	m_z.resize(m_nY);
	//m_z.set_size(nY);

	//
	m_trans_lines_polar.resize(m_nX);
	//m_transData.set_size(nX, ndim);

	return true;
}

Eigen::VectorXd ecmlr::getParameters() const
{
	int num = get_PARAMETER_NUM();
	Eigen::VectorXd parameters(num);
	for (int i = 0;i < num;++i)
	{
		parameters[i] = m_parameters[i];
	}
	return parameters;
}

bool ecmlr::E_Step()
{
	// E-step evaluate the posteriors

	//double outlier = 2.0 * m_delta_2 / (r * r + FLOAT_EPSILON);
	//double outlier = pow(2.0*PI*m_delta_2, ndim / 2) * nX / nY;
	//double outlier_dist = 3.0 * sqrt(m_delta_2);
	//double r = 5.0;
	//double r = outlier_dist;
	//double r = 2.0 / 9.0;
	//double outlier = 3.0 * m_delta_2 / (r * r + DBL_EPSILON);
	double outlier = 2.0 / 9.0;
	//double outlier = 2.0 / 4.0;
	for(int j = 0;j < m_nY;++j)
	{
		double sum = 0.0;
		vector<double> tmp(m_nX);
		double n_ = 1.0 / (double)m_nX;
		for (int i = 0;i < m_nX;++i)
		{
			//fPoint dist = distance2(m_observe_lines_polar[j], m_trans_lines_polar[i]);
			fPoint dist = segment2polarline(m_trans_lines_polar[i], m_reference_lines_polar[j]);
			//fPoint dist;
			//dist.x = point2polarline(m_trans_lines_polar[i].pt1, m_observe_lines_polar[j]);
			//dist.y = point2polarline(m_trans_lines_polar[i].pt2, m_observe_lines_polar[j]);
			double dist_2 = dist.x*dist.x + dist.y*dist.y;

			//////////////////////////////////////////////////////////////////////////
			double centerx1 = (m_trans_lines_polar[i].pt1.x + m_trans_lines_polar[i].pt2.x)*0.5;
			double centery1 = (m_trans_lines_polar[i].pt1.y + m_trans_lines_polar[i].pt2.y)*0.5;
			double centerx2 = (m_reference_lines_polar[j].pt1.x + m_reference_lines_polar[j].pt2.x)*0.5;
			double centery2 = (m_reference_lines_polar[j].pt1.y + m_reference_lines_polar[j].pt2.y)*0.5;
			//fPoint shift = fPoint(centerx1-centerx2, centery1-centery2);
			double shift = sqrt((centerx1-centerx2)*(centerx1-centerx2) + (centery1-centery2)*(centery1-centery2));
			//shift = max(dist.x,dist.y)*exp(-1.0/(shift*shift+DBL_EPSILON));
			shift = 0.0*exp(-m_delta_2/(shift*shift+DBL_EPSILON));
			//////////////////////////////////////////////////////////////////////////
			tmp[i] = exp(-(dist_2 + shift*shift) / (2.0 * m_delta_2));
			sum += tmp[i];
		}
		
		//cout<<sum<<"\t"<<outlier<<endl;
		sum += outlier;
		for(int i = 0;i < m_nX;++i)
		{
			m_alpha(j, i) = tmp[i] / (sum + DBL_EPSILON);
		}
		m_alpha(j, m_nX) = outlier / (sum + DBL_EPSILON);
	}

	//MetropolisHastings();
	// calculate the weights of virtual observations
	// calculate the virtual observation matrix
	//m_w.fill(0.0);
	m_lammda_sum = 0.0;
	for (int i = 0;i < m_nX;++i)
	{
		m_lambda[i] = 0;
		double theta = 0.0;
		double rho = 0.0;
		double x1 = 0.0;
		double y1 = 0.0;
		double x2 = 0.0;
		double y2 = 0.0;

		for (int j = 0;j < m_nY;++j)
		{
			m_lambda[i] += m_alpha(j, i);
			theta += m_alpha(j, i) * m_reference_lines_polar[j].theta;
			rho += m_alpha(j, i) * m_reference_lines_polar[j].rho;
			x1 += m_alpha(j, i) * m_reference_lines_polar[j].pt1.x;
			y1 += m_alpha(j, i) * m_reference_lines_polar[j].pt1.y;
			x2 += m_alpha(j, i) * m_reference_lines_polar[j].pt2.x;
			y2 += m_alpha(j, i) * m_reference_lines_polar[j].pt2.y;
		}
		m_lammda_sum += m_lambda[i];
		fPoint line[2];
		line[0] = fPoint(x1, y1);
		line[1] = fPoint(x2, y2);
		m_w[i] = line2polar(line);
		//m_w[i].theta = theta;
		//m_w[i].rho = rho;
	}
	return true;
}

void ecmlr::funcErrorEquation(double *param, double *hx, int nparameter, int nequation, void *adata)
{
	ecmlr *pThis = (ecmlr*)adata;
	int nX = (int)pThis->m_nX;
	int nY = (int)pThis->m_nY;
	int ndim = (int)pThis->m_ndim;

	int pos = 0;
	int i;

	for(i=0;i<nparameter;++i)
	{
		pThis->m_parameters[i] = param[i];
	}

	for(i = 0;i < nequation;++i)
	{
		cvline_polar outPt = pThis->forward(pThis->m_source_lines_polar[i]);
		double res = pThis->m_lambda[i] * outPt.rho * outPt.rho;
		res += scale_angle * pThis->m_lambda[i] * outPt.theta * outPt.theta;
		res -= 2.0 * outPt.rho * pThis->m_w[i].rho;
		res -= scale_angle * 2.0 * outPt.theta * pThis->m_w[i].theta;
		//double dist = segment_distance(pThis->m_w[i], outPt);
		hx[pos++] = res;
	}
}

void ecmlr::jacErrorEquation(double *param, double *jac, int nparameter, int nequation, void *adata)
{
	ecmlr *pThis = (ecmlr*)adata;
	int nX = (int)pThis->m_nX;
	int nY = (int)pThis->m_nY;
	int ndim = (int)pThis->m_ndim;

	int pos = 0;
	int i;
	for(i=0;i<nparameter;++i)
	{
		pThis->m_parameters[i] = param[i];
	}

	double pstep_scale = 1e-4;
	static double den = 0.5 / pstep_scale;
	int c = 0;

	for(i = 0;i < nequation;++i)
	{
		cvline_polar outPt = pThis->forward(pThis->m_source_lines_polar[i]);
		double temp_theta = 2.0 * (pThis->m_lambda[i] * outPt.theta - pThis->m_w[i].theta);
		double temp_rho = 2.0 * (pThis->m_lambda[i] * outPt.rho - pThis->m_w[i].rho);
		for(int p=0;p<nparameter;++p)
		{
			double middle = pThis->m_parameters[p];
			pThis->m_parameters[p] = middle + pstep_scale;
			cvline_polar outLine1 = pThis->forward(pThis->m_source_lines_polar[i]);

			pThis->m_parameters[p] = middle - pstep_scale;
			cvline_polar outLine2 = pThis->forward(pThis->m_source_lines_polar[i]);

			pThis->m_parameters[p] = middle;

			double derivative_theta = (outLine1.theta - outLine2.theta) * den;
			double derivative_rho = (outLine1.rho - outLine2.rho) * den;

			jac[c++] = derivative_rho * temp_rho + scale_angle * derivative_theta * temp_theta;
		}
	}
}

void ecmlr::setParameters(const vector<double> & parameters)
{
	m_parameters = parameters;
}

fPoint ecmlr::distance2(const cvline_polar& l1, const cvline_polar& l2)
{
	double centerx1 = (l1.pt1.x + l1.pt2.x)*0.5;
	double centery1 = (l1.pt1.y + l1.pt2.y)*0.5;
	double centerx2 = (l2.pt1.x + l2.pt2.x)*0.5;
	double centery2 = (l2.pt1.y + l2.pt2.y)*0.5;
	//fPoint shift = fPoint(centerx1-centerx2, centery1-centery2);
	double shift = sqrt((centerx1-centerx2)*(centerx1-centerx2) + (centery1-centery2)*(centery1-centery2));
	//double shift_weight = 1.0E-2;
	//double shift_weight = 0.08;

	fPoint dist;
	// dist1
	dist.x = distance(l1, l2);
	// dist2
	dist.y = distance(l2, l1);

	//shift = max(dist.x,dist.y)*exp(-1.0/(shift*shift+DBL_EPSILON));
	//shift = 3.0*exp(-m_delta_2/(shift*shift+DBL_EPSILON));
	shift = 0.0*exp(-1.0/(shift+DBL_EPSILON));
	//if((dist.x*dist.x+dist.y*dist.y) < 4.0)
	//{
	dist.x += shift;
	dist.y += shift;
	//}

	//dist.x = fabs(l1.rho - l2.rho);
	//dist.y = scale_angle * fabs(l1.theta - l2.theta);
	////dist.y = dist.y > lammda_rho*dist.x ? lammda_rho*dist.x : dist.y;
	////dist.y = fabs(tan(l1.theta) - tan(l2.theta));
	return dist;
}

double ecmlr::point2polarline(const fPoint& p, const cvline_polar& l)
{
	double dist = -p.x * sin(l.theta) + p.y * cos(l.theta) - l.rho;
	return dist;
}

fPoint ecmlr::segment2polarline(const cvline_polar& segment, const cvline_polar& line)
{
	return distance2(line, segment);
	fPoint dist(point2polarline(segment.pt1, line), point2polarline(segment.pt2, line));
	return dist;
}

double ecmlr::distance(const cvline_polar& l1, const cvline_polar& l2)
{
	//return (point2polarline(l1.pt1, l2) + point2polarline(l1.pt2, l2))*0.5;
	double d1 = point2polarline(l1.pt1, l2);
	double d2 = point2polarline(l1.pt2, l2);
	return sqrt(d1*d1+d2*d2);
	//return (fabs(d1)+fabs(d2))*0.5;
	// dist1
	fPoint gpt1 = l2.pt1;
	fPoint gpt2 = l2.pt2;

	double hemline = sqrt((gpt1.x - gpt2.x)*(gpt1.x - gpt2.x) + (gpt1.y - gpt2.y)*(gpt1.y - gpt2.y));
	double k1 = (gpt1.y - gpt2.y) / (hemline + DBL_EPSILON);
	double k2 = (gpt1.x - gpt2.x) / (hemline + DBL_EPSILON);

	//for the first endpoint of a line
	double dist1 = fabs((gpt2.x * gpt1.y - gpt1.x * gpt2.y) / (hemline + DBL_EPSILON) - k1 * l1.pt1.x + k2 * l1.pt1.y);

	//for the second endpoint of a line
	double dist2 = fabs((gpt2.x * gpt1.y - gpt1.x * gpt2.y) / (hemline + DBL_EPSILON) - k1 * l1.pt2.x + k2 * l1.pt2.y);

	// dist2
	return (dist1+dist2)/2.0;
}

void ecmlr::updateTransLines()
{
	for (int i = 0;i < m_nX;++i)
	{
		cvline_polar transedPt = forward(m_source_lines_polar[i]);
		m_trans_lines_polar[i] = transedPt;
	}
}
void ecmlr::levmar_function_fvec(double *param, double *hx, int nparameter, int nequation, void *adata)
{
	ecmlr *pThis = (ecmlr*)adata;
	int nX = (int)pThis->m_nX;
	int nY = (int)pThis->m_nY;
	int ndim = (int)pThis->m_ndim;

	int pos = 0;
	int i;
	for(i=0;i<nparameter;++i)
	{
		pThis->m_parameters[i] = param[i];
		//pThis->m_parameters[i] = param[i];
	}

	int c = 0;

	for (int i = 0;i < nX;++i)
	{
		cvline_polar outPt = pThis->forward(pThis->m_source_lines_polar[i]);
		for (int j = 0;j < nY;++j)
		{
			//fPoint dist = pThis->distance2(pThis->m_observe_lines_polar[j], outPt);
			fPoint dist = pThis->segment2polarline(outPt, pThis->m_reference_lines_polar[j]);

			hx[c++] = pThis->m_alpha(j,i) * dist.x;
			hx[c++] = pThis->m_alpha(j,i) * dist.y;
		}
	}
}

void ecmlr::levmar_function_jac(double *param, double *jac, int nparameter, int nequation, void *adata)
{
	ecmlr *pThis = (ecmlr*)adata;
	int nX = (int)pThis->m_nX;
	int nY = (int)pThis->m_nY;
	int ndim = (int)pThis->m_ndim;

	int pos = 0;
	int i;
	for(i=0;i<nparameter;++i)
	{
		pThis->m_parameters[i] = param[i];
		//pThis->m_parameters[i] = param[i];
	}

	double pstep_scale = 1e-4;
	static double den = 0.5 / pstep_scale;
	int c = 0;

	for (int i = 0;i < nX;++i)
	{
		cvline_polar outPt = pThis->forward(pThis->m_source_lines_polar[i]);
		for (int j = 0;j < nY;++j)
		{
			for(int p=0;p<nparameter;++p)
			{
				double middle = pThis->m_parameters[p];
				pThis->m_parameters[p] = middle + pstep_scale;
				cvline_polar outLine1 = pThis->forward(pThis->m_source_lines_polar[i]);
				//fPoint dist1 = pThis->distance2(pThis->m_observe_lines_polar[j], outLine1);
				fPoint dist1 = pThis->segment2polarline(outLine1, pThis->m_reference_lines_polar[j]);

				pThis->m_parameters[p] = middle - pstep_scale;
				cvline_polar outLine2 = pThis->forward(pThis->m_source_lines_polar[i]);
				//fPoint dist2 = pThis->distance2(pThis->m_observe_lines_polar[j], outLine2);
				fPoint dist2 = pThis->segment2polarline(outLine2, pThis->m_reference_lines_polar[j]);

				pThis->m_parameters[p] = middle;

				double derivative_x = (dist1.x - dist2.x) * den;
				double derivative_y = (dist1.y - dist2.y) * den;

				jac[c*nparameter+p] = pThis->m_alpha(j,i) * derivative_x;
				jac[(c+1)*nparameter+p] = pThis->m_alpha(j,i) * derivative_y;
			}
			c += 2;
		}
	}
}

bool ecmlr::parameters_optimization()
{
	// levmar
	return levmar_optimization();

	// spase LM
	//return spaseLM_optimization();

	// alglib improved LM
	//return alglib_optimization();

	// Eigen
	return eigen_levmar_optimization();


	int nX = m_nX;
	int nY = m_nY;
	int ndim = m_ndim;

	int nparam = get_PARAMETER_NUM();
	std::vector<double> cparm(nparam);

	for(int i=0;i<nparam;++i) cparm[i] = m_parameters[i];

	//Eigen::VectorXd outParameters(nparam);
	//// lm
	//double *p = &cparm[0];
	////double *x = new double[2*nX];
	////vector<double> x(2*nX, 0.0);
	////for(int i=0; i<nX*nY; i++) x[i]=0.0;

	//double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
	//opts[0]=LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-20;
	//opts[4]= LM_DIFF_DELTA; // relevant only if the Jacobian is approximated using finite differences; specifies forward differencing 

	////int ret = dlevmar_dif(levmar_function_fvec, &cparm[0], &x[0], nparam, 2*nX, 1000, opts, info, NULL, NULL, this);  // no Jacobian
	//int ret = dlevmar_der(levmar_function_fvec, levmar_function_jac, p, NULL, nparam, 2*nX*nY, 1000, opts, info, NULL, NULL, this); // with analytic Jacobian
	////delete []x;

	
	return true;
}

bool ecmlr::levmar_optimization()
{
	int nX = m_nX;
	int nY = m_nY;
	int ndim = m_ndim;

	int nparam = get_PARAMETER_NUM();
	std::vector<double> cparm(nparam);

	for(int i=0;i<nparam;++i) cparm[i] = m_parameters[i];

	// lm
	double *p = &cparm[0];
	//double *x = new double[2*nX];
	//vector<double> x(2*nX, 0.0);
	//for(int i=0; i<nX*nY; i++) x[i]=0.0;

	double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
	opts[0]=LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-20;
	opts[4]= LM_DIFF_DELTA; // relevant only if the Jacobian is approximated using finite differences; specifies forward differencing 

	//int ret = dlevmar_dif(levmar_function_fvec, &cparm[0], &x[0], nparam, 2*nX, 1000, opts, info, NULL, NULL, this);  // no Jacobian
	int ret = dlevmar_der(levmar_function_fvec, levmar_function_jac, p, NULL, nparam, 2*nX*nY, 1000, opts, info, NULL, NULL, this); // with analytic Jacobian
	//delete []x;

	for(int i=0;i<nparam;++i) m_parameters[i] = p[i];
	return true;
}

bool ecmlr::eigen_levmar_optimization()
{
	//int nX = m_nX;
	//int nY = m_nY;
	//int ndim = m_ndim;

	//int nparam = get_PARAMETER_NUM();
	//
	//int info;
	////double fnorm, covfac;
	//Eigen::VectorXd x(nparam);
	//for(int i=0;i<nparam;++i) x[i] = m_parameters[i];
	//
	//// do the computation
	//lmder_functor functor(nparam, 2*nX*nY);
	//functor.pData = this;
	//Eigen::LevenbergMarquardt<lmder_functor> lm(functor);
	//info = lm.minimize(x);
	//
	//for(int i=0;i<nparam;++i) m_parameters[i] = x[i];
	return true;
}

//
//int lmder_functor::operator()(const VectorXd &x, VectorXd &fvec) const
//{
//	ecmlr *pThis = (ecmlr*)pData;
//	int nX = (int)pThis->m_nX;
//	int nY = (int)pThis->m_nY;
//	int ndim = (int)pThis->m_ndim;
//
//	int nparameter = pThis->get_PARAMETER_NUM();
//
//	for(int i=0;i<nparameter;++i)
//	{
//		pThis->m_parameters[i] = x[i];
//	}
//
//	double pstep_scale = 1e-4;
//	static double den = 0.5 / pstep_scale;
//
//	int c = 0;
//	for (int i = 0;i < nX;++i)
//	{
//		cvline_polar outPt = pThis->forward(pThis->m_source_lines_polar[i]);
//		for (int j = 0;j < nY;++j)
//		{
//			fPoint dist = pThis->segment2polarline(outPt, pThis->m_reference_lines_polar[j]);
//
//			fvec[c++] = pThis->m_alpha(j,i) * dist.x;
//			fvec[c++] = pThis->m_alpha(j,i) * dist.y;
//		}
//	}
//	return 0;
//}
//
//int lmder_functor::df(const VectorXd &x, MatrixXd &fjac) const
//{
//	ecmlr *pThis = (ecmlr*)pData;
//	int nX = (int)pThis->m_nX;
//	int nY = (int)pThis->m_nY;
//	int ndim = (int)pThis->m_ndim;
//
//	int nparameter = pThis->get_PARAMETER_NUM();
//
//	for(int i=0;i<nparameter;++i)
//	{
//		pThis->m_parameters[i] = x[i];
//	}
//
//	double pstep_scale = 1e-4;
//	static double den = 0.5 / pstep_scale;
//	int c = 0;
//
//
//	for (int i = 0;i < nX;++i)
//	{
//		cvline_polar outPt = pThis->forward(pThis->m_source_lines_polar[i]);
//		for (int j = 0;j < nY;++j)
//		{
//			//fPoint dist = pThis->distance2(pThis->m_observe_lines_polar[j], outPt);
//			fPoint dist = pThis->segment2polarline(outPt, pThis->m_reference_lines_polar[j]);
//
//			for(int p=0;p<nparameter;++p)
//			{
//				double middle = pThis->m_parameters[p];
//				pThis->m_parameters[p] = middle + pstep_scale;
//				cvline_polar outLine1 = pThis->forward(pThis->m_source_lines_polar[i]);
//				//fPoint dist1 = pThis->distance2(pThis->m_observe_lines_polar[j], outLine1);
//				fPoint dist1 = pThis->segment2polarline(outLine1, pThis->m_reference_lines_polar[j]);
//
//				pThis->m_parameters[p] = middle - pstep_scale;
//				cvline_polar outLine2 = pThis->forward(pThis->m_source_lines_polar[i]);
//				//fPoint dist2 = pThis->distance2(pThis->m_observe_lines_polar[j], outLine2);
//				fPoint dist2 = pThis->segment2polarline(outLine2, pThis->m_reference_lines_polar[j]);
//
//				pThis->m_parameters[p] = middle;
//
//				double derivative_x = (dist1.x - dist2.x) * den;
//				double derivative_y = (dist1.y - dist2.y) * den;
//
//				fjac(c,p) = pThis->m_alpha(j,i) * derivative_x;
//				fjac(c+1,p) = pThis->m_alpha(j,i) * derivative_y;
//			}
//			c += 2;
//		}
//	}
//	return 0;
//}
