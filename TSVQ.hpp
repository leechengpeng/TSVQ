#pragma once
#include <vector>

template <typename T>
class CTSVQ
{
#if defined(DEBUG) | defined(_DEBUG)
	template <typename T> friend void isClusterRight(const CTSVQ<T> &vTSVQ);
#endif

public:
	CTSVQ(void);
	~CTSVQ(void);
	
	void quantizeVectors(std::vector<T*> vInputDataSet, unsigned int vDimension);

	inline void quantizeVectors(std::vector<T*> vInputDataSet, unsigned int vDimension, unsigned int vCodeVectorsNums);
	inline void quantizeVectors(std::vector<T*> vInputDataSet, unsigned int vDimension, unsigned int vCodeVectorsNums, unsigned int vMaxInterations);
	inline void quantizeVectors(std::vector<T*> vInputDataSet, unsigned int vDimension, unsigned int vCodeVectorsNums, unsigned int vMaxInterations, double vEpsilon);

	inline void clear();

	const T* findMostSimilarVectors(const T *vVectors);
	const T* findCodeVectors(const T *vVectors);

private:
	struct SNode
	{
		SNode() : CodeVectors(NULL), pLeft(NULL), pRight(NULL) {}
		~SNode() 
		{
			if (pLeft != NULL)
			{
				delete pLeft;
				pLeft = NULL;
			}

			if (pRight != NULL)
			{
				delete pRight;
				pRight = NULL;
			}

			if (CodeVectors != NULL)
			{
				delete[] CodeVectors;
				CodeVectors = NULL;
			}
		}

		SNode* pLeft;
		SNode* pRight;

		T* CodeVectors;
		std::vector<T*> VectorsSet;
	};

private:
	void __initRootNode();
	void __split(std::vector<SNode*> &vSplitNodeSet);
	void __clusterInputVectorsSet(std::vector<SNode*> &vSplitNodeSet);
	void __calCentroid(const std::vector<T*> &vVectorsSet, T *voCentroid);
	void __iterate(std::vector<SNode*> &vSplitNodeSet, const double vDistortionMeasure);

	double __calDistortionMeasure(const std::vector<SNode*> &vNodeSet);
	double __calEuclideanDistance(const T *vVectorsA, const T *vVectorsB);

private:
	double m_Epsilon;
	SNode* m_RootNode;

	unsigned int m_Dimension;
	unsigned int m_MaxInterations; 
	unsigned int m_CodeVectorsNums;

	std::vector<T*> m_InputVectorsSet;
};

template <typename T>
CTSVQ<T>::CTSVQ(void) : m_Dimension(0), m_MaxInterations(10), m_CodeVectorsNums(1), m_Epsilon(0.001), m_RootNode(NULL)
{

}

template <typename T>
CTSVQ<T>::~CTSVQ(void)
{
	clear();
}

//******************************************************************************
//FUNCTION:
template <typename T>
void CTSVQ<T>::clear()
{
	m_Epsilon		  = 0.001;
	m_Dimension		  = 0;
	m_MaxInterations  = 10;
	m_CodeVectorsNums = 1;

	if (m_RootNode != NULL)
	{
		delete m_RootNode;
		m_RootNode = NULL;
	}
}

//******************************************************************************
//FUNCTION:
template <typename T>
void CTSVQ<T>::__initRootNode()
{
	m_RootNode = new SNode();

	m_RootNode->VectorsSet  = m_InputVectorsSet;
	m_RootNode->CodeVectors = new T[m_Dimension]();

	__calCentroid(m_InputVectorsSet, m_RootNode->CodeVectors);
}

//******************************************************************************
//FUNCTION:
template <typename T>
const T* CTSVQ<T>::findMostSimilarVectors(const T *vVectors)
{
	_ASSERT(vVectors);

	auto Node = m_RootNode;

	while (Node->pLeft != NULL && Node->pRight != NULL)
	{
		auto LeftDistance  = __calEuclideanDistance(vVectors, Node->pLeft->CodeVectors);
		auto RightDistance = __calEuclideanDistance(vVectors, Node->pRight->CodeVectors);

		LeftDistance < RightDistance ? Node = Node->pLeft : Node = Node->pRight;
	}

	auto MostSimilarVectors = Node->VectorsSet.at(0);
	auto MostSimilarity		= __calEuclideanDistance(vVectors, MostSimilarVectors);
	
	for (const auto NodeVectors : Node->VectorsSet)
	{
		auto Similarity = __calEuclideanDistance(vVectors, NodeVectors);
		if (Similarity < MostSimilarity)
		{
			MostSimilarity	   = Similarity;
			MostSimilarVectors = NodeVectors;
		}
	}

	return MostSimilarVectors;
}

//******************************************************************************
//FUNCTION:
template <typename T>
const T* CTSVQ<T>::findCodeVectors(const T *vVectors)
{
	_ASSERT(vVectors);

	auto Node = m_RootNode;

	// PIXME: if there is no LeafNode?
	while (Node->pLeft != NULL && Node->pRight != NULL)
	{
		auto LeftDistance  = __calEuclideanDistance(vVectors, Node->pLeft->CodeVectors);
		auto RightDistance = __calEuclideanDistance(vVectors, Node->pRight->CodeVectors);

		LeftDistance < RightDistance ? Node = Node->pLeft : Node = Node->pRight;
	}

	return Node->CodeVectors;
}

//******************************************************************************
//FUNCTION:
template <typename T>
void CTSVQ<T>::quantizeVectors(std::vector<T*> vInputDataSet, unsigned int vDimension)
{
	_ASSERT(vInputDataSet.size() > 0);

	m_Dimension		  = vDimension;
	m_InputVectorsSet = vInputDataSet;

	__initRootNode();

	std::vector<SNode*> NodeSet;
	NodeSet.push_back(m_RootNode);

	auto DistortionMeasure = __calDistortionMeasure(NodeSet);

	int MaxLevel = static_cast<int>(log(m_CodeVectorsNums)/log(2.0));
	for (int TreeLevel=0; TreeLevel<MaxLevel; ++TreeLevel)
	{
		__split(NodeSet);

		std::vector<SNode*> SplitNodeSet;
		for (auto &Node : NodeSet)
		{
			SplitNodeSet.push_back(Node->pLeft);
			SplitNodeSet.push_back(Node->pRight);
		}

		__iterate(SplitNodeSet, DistortionMeasure);

		swap(SplitNodeSet, NodeSet);
	}
}

//******************************************************************************
//FUNCTION:
template <typename T>
void CTSVQ<T>::quantizeVectors(std::vector<T*> vInputDataSet, unsigned int vDimension, unsigned int vCodeVectorsNums)
{
	_ASSERT(vCodeVectorsNums != 0);
	m_CodeVectorsNums = vCodeVectorsNums;

	quantizeVectors(vInputDataSet, vDimension);
}

//******************************************************************************
//FUNCTION:
template <typename T>
void CTSVQ<T>::quantizeVectors(std::vector<T*> vInputDataSet, unsigned int vDimension, unsigned int vCodeVectorsNums, unsigned int vMaxInterations)
{
	
	m_MaxInterations = vMaxInterations;

	quantizeVectors(vInputDataSet, vDimension, vCodeVectorsNums);
}

//******************************************************************************
//FUNCTION:
template <typename T>
void CTSVQ<T>::quantizeVectors(std::vector<T*> vInputDataSet, unsigned int vDimension, unsigned int vCodeVectorsNums, unsigned int vMaxInterations, double vEpsilon)
{
	_ASSERT(vEpsilon != 0);
	m_Epsilon = vEpsilon;

	quantizeVectors(vInputDataSet, vDimension, vCodeVectorsNums, vMaxInterations);
}

//******************************************************************************
//FUNCTION:
template <typename T>
void CTSVQ<T>::__split(std::vector<SNode*> &vSplitNodeSet)
{
	for (auto &Node : vSplitNodeSet)
	{
		Node->pLeft  = new SNode();
		Node->pRight = new SNode();

		Node->pLeft->CodeVectors  = new T[m_Dimension]();
		Node->pRight->CodeVectors = new T[m_Dimension]();

		auto CodeVectors = Node->CodeVectors;

		for (unsigned int Index=0; Index<m_Dimension; ++Index)
		{
			Node->pLeft->CodeVectors[Index]  = static_cast<T>(CodeVectors[Index] * (1 + m_Epsilon));
			Node->pRight->CodeVectors[Index] = static_cast<T>(CodeVectors[Index] * (1 - m_Epsilon));
		}
	}
}

//******************************************************************************
//FUNCTION:
template <typename T>
void CTSVQ<T>::__clusterInputVectorsSet(std::vector<SNode*> &vSplitNodeSet)
{
	for (const auto &Vectors : m_InputVectorsSet)
	{
		unsigned int Index = 0;

		auto MinIndex	 = Index;
		auto MinDistance = __calEuclideanDistance(Vectors, vSplitNodeSet.at(Index)->CodeVectors);

		for (++Index; Index<vSplitNodeSet.size(); ++Index)
		{
			auto Distance = __calEuclideanDistance(Vectors, vSplitNodeSet.at(Index)->CodeVectors);

			if (Distance < MinDistance)
			{
				MinIndex	= Index;
				MinDistance = Distance;
			}
		}

		vSplitNodeSet.at(MinIndex)->VectorsSet.push_back(Vectors);
	}
}

//******************************************************************************
//FUNCTION:
template <typename T>
void CTSVQ<T>::__iterate(std::vector<SNode*> &vSplitNodeSet, const double vDistortionMeasure)
{
	double CurrentDistortionMeasure  = 0.0;
	double PreviousDistortionMeasure = vDistortionMeasure;

	for (unsigned int IterCounter=0; IterCounter<m_MaxInterations; ++IterCounter)
	{
		__clusterInputVectorsSet(vSplitNodeSet);

		for (auto &Node : vSplitNodeSet)
		{
			__calCentroid(Node->VectorsSet, Node->CodeVectors);
		}

		CurrentDistortionMeasure = __calDistortionMeasure(vSplitNodeSet);

		if (PreviousDistortionMeasure - CurrentDistortionMeasure > m_Epsilon * PreviousDistortionMeasure)
		{
			for (auto &Node : vSplitNodeSet) Node->VectorsSet.clear();

			PreviousDistortionMeasure = CurrentDistortionMeasure;
		}
		else 
		{
			break;
		}
	}
}

//******************************************************************************
//FUNCTION:
template <typename T>
void CTSVQ<T>::__calCentroid(const std::vector<T*> &vVectorsSet, T *voCentroid)
{
	_ASSERT(voCentroid);

	if (vVectorsSet.size() == 0) return;

	double *Centroid = new double[m_Dimension]();

	for (const auto &Vectors : vVectorsSet)
	{
		for (unsigned int Index=0; Index<m_Dimension; ++Index)
		{
			Centroid[Index] += Vectors[Index];
		}
	}

	for (unsigned int Index=0; Index<m_Dimension; ++Index)
	{
		voCentroid[Index] = static_cast<T>(Centroid[Index]/vVectorsSet.size());
	}

	delete[] Centroid;
}

//******************************************************************************
//FUNCTION:
template <typename T>
double CTSVQ<T>::__calDistortionMeasure(const std::vector<SNode*> &vNodeSet)
{
	double EuclideanDistanceSum = 0.0;

	for (const auto &Node : vNodeSet)
	{
		for (const auto &Vectors : Node->VectorsSet)
		{
			EuclideanDistanceSum += __calEuclideanDistance(Vectors, Node->CodeVectors);
		}
	}

	_ASSERT(m_InputVectorsSet.size() != 0 && m_Dimension != 0);

	return EuclideanDistanceSum/(m_InputVectorsSet.size()*m_Dimension);
}

//******************************************************************************
//FUNCTION:
template <typename T>
double CTSVQ<T>::__calEuclideanDistance(const T *vVectorsA, const T *vVectorsB)
{
	_ASSERT(vVectorsA && vVectorsB);

	double Sum = 0.0;

	for(unsigned int Index=0; Index<m_Dimension; ++Index)
	{
		Sum += (vVectorsA[Index] - vVectorsB[Index])*(vVectorsA[Index] - vVectorsB[Index]);
	}

	return sqrt(Sum);
}

//******************************************************************************
//FUNCTION:
template <typename T>
double calEuclideanDistance(const T *vVectorsA, const T *vVectorsB, unsigned int vDim)
{
	_ASSERT(vVectorsA && vVectorsB);

	double Sum = 0.0;

	for(unsigned int Index=0; Index<vDim; ++Index)
	{
		Sum += (vVectorsA[Index] - vVectorsB[Index])*(vVectorsA[Index] - vVectorsB[Index]);
	}

	return sqrt(Sum);
}

//******************************************************************************
//FUNCTION:
template <typename T> void isClusterRight(const CTSVQ<T> &vTSVQ)
{
	auto RootNode = vTSVQ.m_RootNode;

	std::vector<decltype(RootNode)> NodeSet;
	NodeSet.push_back(RootNode);

	while (NodeSet.at(0)->pLeft != NULL && NodeSet.at(0)->pRight != NULL)
	{
		std::vector<decltype(RootNode)> NextNodeSet;
		for (auto &Node : NodeSet)
		{
			NextNodeSet.push_back(Node->pLeft);
			NextNodeSet.push_back(Node->pRight);
		}

		NodeSet.swap(NextNodeSet);
	}

	for (auto &Node : NodeSet)
	{
		for (auto &Vectors : Node->VectorsSet)
		{
			std::cout << "Vectors: ( ";
			for (unsigned int i=0; i<vTSVQ.m_Dimension; ++i)
			{
				std::cout << Vectors[i] << " ";
			}
			std::cout << ")";

			std::cout << "\nIts CodeVectors: ( ";
			for (unsigned int i=0; i<vTSVQ.m_Dimension; ++i)
			{
				std::cout << Node->CodeVectors[i] << " ";
			}
			std::cout << ")" << std::endl << std::endl;

			bool IsRight = true;
			auto MinDistance = calEuclideanDistance(Vectors, Node->CodeVectors, vTSVQ.m_Dimension);
			for (auto &TempNode : NodeSet)
			{
				auto Distance = calEuclideanDistance(Vectors, TempNode->CodeVectors, vTSVQ.m_Dimension);
				if (Distance < MinDistance)
				{
					MinDistance = Distance;
					IsRight = false;
					break;
				}

				std::cout << "Distance between ( ";
				for (unsigned int i=0; i<vTSVQ.m_Dimension; ++i)
				{
					std::cout << Vectors[i] << " "; 
				}
				std::cout << ") and ( ";
				for (unsigned int i=0; i<vTSVQ.m_Dimension; ++i)
				{
					std::cout << TempNode->CodeVectors[i] << " "; 
				}
				std::cout << ") is " << Distance << std::endl;
			}
			std::cout << "\nMinDistance = " << MinDistance << std::endl;
			std::cout << "\nIs cluster right: ";
			IsRight ? std::cout << "True" : std::cout << "False";

			std::cout << "\n\n******************************************************************\n" << std::endl;
		}
	}
}