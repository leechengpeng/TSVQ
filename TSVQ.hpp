#pragma once
#include <vector>

namespace LLL
{
	template <typename T, unsigned Dimension>
	struct _SNode
	{
		_SNode() : pLeft(NULL), pRight(NULL) {}
		~_SNode() 
		{
			delete pLeft;
			delete pRight;

			pLeft  = NULL;
			pRight = NULL;
		}

		_SNode* pLeft;
		_SNode* pRight;

		T CodeVectors[Dimension];
		std::vector<T*> VectorsSet;
	};

	template <typename T, unsigned Dimension>
	class TSVQ
	{
	public:
		TSVQ(void);
		~TSVQ(void);
	
		void build(const std::vector<T*>& vVectorSet);
		void build(const std::vector<T*>& vVectorSet, unsigned vCodeVectorsNums);
		void build(const std::vector<T*>& vVectorSet, unsigned vCodeVectorsNums, unsigned vMaxInterations);
		void build(const std::vector<T*>& vVectorSet, unsigned vCodeVectorsNums, unsigned vMaxInterations, double vEpsilon);

		const T* quantizeVector(const T *vVector) const;

	private:
		typedef _SNode<T, Dimension> SNode;

		void __split(std::vector<SNode*> &vSplitNodeSet);
		void __clusterInputVectorsSet(const std::vector<T*>& vVectorSet, std::vector<SNode*>& vSplitNodeSet);
		void __calCentroid(const std::vector<T*> &vVectorsSet, T *voCentroid);
		void __iterate(const std::vector<T*>& vVectorSet, std::vector<SNode*> &vSplitNodeSet, const double vDistortionMeasure);

		double __calDistortionMeasure(const std::vector<SNode*> &vNodeSet, unsigned vNumVectors) const;
		double __distance(const T *vVectorA, const T *vVectorB) const;

		double m_Epsilon;
		SNode* m_RootNode;

		unsigned int m_MaxInterations; 
		unsigned int m_CodeVectorsNums;
	};

	template <typename T, unsigned Dimension>
	TSVQ<T, Dimension>::TSVQ(void) : m_MaxInterations(10), m_CodeVectorsNums(1), m_Epsilon(0.001), m_RootNode(NULL)
	{

	}

	template <typename T, unsigned Dimension>
	TSVQ<T, Dimension>::~TSVQ(void)
	{
		m_Epsilon		  = 0.001;
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
	template <typename T, unsigned Dimension>
	void TSVQ<T, Dimension>::build(const std::vector<T*>& vVectorSet)
	{
		_ASSERT(vVectorSet.size() && Dimension);

		m_RootNode = new SNode();
		m_RootNode->VectorsSet = vVectorSet;
		__calCentroid(vVectorSet, m_RootNode->CodeVectors);

		std::vector<SNode*> NodeSet;
		NodeSet.push_back(m_RootNode);

		auto DistortionMeasure = __calDistortionMeasure(NodeSet, vVectorSet.size());

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

			__iterate(vVectorSet, SplitNodeSet, DistortionMeasure);

			swap(SplitNodeSet, NodeSet);
		}
	}

	//******************************************************************************
	//FUNCTION:
	template <typename T, unsigned Dimension>
	inline void TSVQ<T, Dimension>::build(const std::vector<T*>& vVectorSet, unsigned int vCodeVectorsNums)
	{
		_ASSERT(vCodeVectorsNums != 0);
		m_CodeVectorsNums = vCodeVectorsNums;

		build(vVectorSet);
	}

	//******************************************************************************
	//FUNCTION:
	template <typename T, unsigned Dimension>
	inline void TSVQ<T, Dimension>::build(const std::vector<T*>& vVectorSet, unsigned int vCodeVectorsNums, unsigned int vMaxInterations)
	{
	
		m_MaxInterations = vMaxInterations;

		build(vVectorSet, vCodeVectorsNums);
	}

	//******************************************************************************
	//FUNCTION:
	template <typename T, unsigned Dimension>
	inline void TSVQ<T, Dimension>::build(const std::vector<T*>& vVectorSet, unsigned int vCodeVectorsNums, unsigned int vMaxInterations, double vEpsilon)
	{
		_ASSERT(vEpsilon != 0);
		m_Epsilon = vEpsilon;

		build(vVectorSet, vCodeVectorsNums, vMaxInterations);
	}

	//******************************************************************************
	//FUNCTION:
	template <typename T, unsigned Dimension>
	const T* TSVQ<T, Dimension>::quantizeVector(const T *vVector) const
	{
		_ASSERT(vVector);

		auto Node = m_RootNode;

		// PIXME: if there is no LeafNode?
		while (Node->pLeft != NULL && Node->pRight != NULL)
		{
			auto LeftDistance  = __distance(vVector, Node->pLeft->CodeVectors);
			auto RightDistance = __distance(vVector, Node->pRight->CodeVectors);

			LeftDistance < RightDistance ? Node = Node->pLeft : Node = Node->pRight;
		}

		return Node->CodeVectors;
	}

	//******************************************************************************
	//FUNCTION:
	template <typename T, unsigned Dimension>
	void TSVQ<T, Dimension>::__split(std::vector<SNode*> &vSplitNodeSet)
	{
		for (auto &Node : vSplitNodeSet)
		{
			Node->pLeft  = new SNode();
			Node->pRight = new SNode();

			auto CodeVectors = Node->CodeVectors;
			for (unsigned int Index=0; Index<Dimension; ++Index)
			{
				Node->pLeft->CodeVectors[Index]  = static_cast<T>(CodeVectors[Index] * (1 + m_Epsilon));
				Node->pRight->CodeVectors[Index] = static_cast<T>(CodeVectors[Index] * (1 - m_Epsilon));
			}
		}
	}

	//******************************************************************************
	//FUNCTION:
	template <typename T, unsigned Dimension>
	void TSVQ<T, Dimension>::__clusterInputVectorsSet(const std::vector<T*>& vVectorSet, std::vector<SNode*>& vSplitNodeSet)
	{
		for (const auto &Vectors : vVectorSet)
		{
			unsigned int Index = 0;

			auto MinIndex	 = Index;
			auto MinDistance = __distance(Vectors, vSplitNodeSet.at(Index)->CodeVectors);

			for (++Index; Index<vSplitNodeSet.size(); ++Index)
			{
				auto Distance = __distance(Vectors, vSplitNodeSet.at(Index)->CodeVectors);

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
	template <typename T, unsigned Dimension>
	void TSVQ<T, Dimension>::__iterate(const std::vector<T*>& vVectorSet, std::vector<SNode*> &vSplitNodeSet, const double vDistortionMeasure)
	{
		double CurrentDistortionMeasure  = 0.0;
		double PreviousDistortionMeasure = vDistortionMeasure;

		for (unsigned int IterCounter=0; IterCounter<m_MaxInterations; ++IterCounter)
		{
			__clusterInputVectorsSet(vVectorSet, vSplitNodeSet);

			for (auto &Node : vSplitNodeSet)
			{
				__calCentroid(Node->VectorsSet, Node->CodeVectors);
			}

			CurrentDistortionMeasure = __calDistortionMeasure(vSplitNodeSet, vVectorSet.size());

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
	template <typename T, unsigned Dimension>
	void TSVQ<T, Dimension>::__calCentroid(const std::vector<T*> &vVectorsSet, T *voCentroid)
	{
		_ASSERT(voCentroid);

		if (vVectorsSet.size() == 0) return;

		double *Centroid = new double[Dimension]();

		for (const auto &Vectors : vVectorsSet)
		{
			for (unsigned int Index=0; Index<Dimension; ++Index)
			{
				Centroid[Index] += Vectors[Index];
			}
		}

		for (unsigned int Index=0; Index<Dimension; ++Index)
		{
			voCentroid[Index] = static_cast<T>(Centroid[Index]/vVectorsSet.size());
		}

		delete[] Centroid;
	}

	//******************************************************************************
	//FUNCTION:
	template <typename T, unsigned Dimension>
	double TSVQ<T, Dimension>::__calDistortionMeasure(const std::vector<SNode*> &vNodeSet, unsigned vNumVectors) const
	{
		_ASSERT(vNumVectors && Dimension);

		double EuclideanDistanceSum = 0.0;
		for (auto& Node : vNodeSet)
		{
			for (auto Vector : Node->VectorsSet)
			{
				EuclideanDistanceSum += __distance(Vector, Node->CodeVectors);
			}
		}

		return EuclideanDistanceSum / (vNumVectors * Dimension);
	}

	//******************************************************************************
	//FUNCTION:
	template <typename T, unsigned Dimension>
	double TSVQ<T, Dimension>::__distance(const T *vVectorA, const T *vVectorB) const
	{
		_ASSERT(vVectorA && vVectorB);

		double Sum = 0.0;
		for (unsigned Index=0; Index < Dimension; ++Index)
		{
			Sum += (vVectorA[Index] - vVectorB[Index]) * (vVectorA[Index] - vVectorB[Index]);
		}

		return sqrt(Sum);
	}
}