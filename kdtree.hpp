#ifndef _SPATIAL_KDTREE_H_
#define _SPATIAL_KDTREE_H_

#include <vector>

//#include <mlib/math/glm.hpp>
//#include <win32/Assert.h>



#ifndef ERROR_ASSERT
#ifdef _DEBUG
#define ERROR_ASSERT(X) assert(X)
#else
#define ERROR_ASSERT(X)
#endif
#endif

namespace spatial {

	namespace helper {
		template <typename T>
		int MinIndex(T const* const seq, int const length, int N, int* minIndex, T* min) {
			for (int i = 0; i < N; ++i)
			{
				min[i] = FLT_MAX;
				minIndex[i] = -1;
			}

			for (int i = 0; i < length; ++i) 
			{
				for (int j = 0; j < N; ++j)
				{
					if (seq[i] < min[j])
					{
						for (int h = N-1; h > j; --h)
						{
							min[h] = min[h-1];
							minIndex[h] = minIndex[h-1];
						}
						min[j] = seq[i];
						minIndex[j] = i;
						break;
					}
				}
			}
			return minIndex[0];
		}


		template<typename T>
		void FindNearest(
			T const* const list
			, int const length
			, T const& p
			, int N
			, int* nearestI
			, float* minDist
			, float* distances
			)
		{
			ERROR_ASSERT(length > 0);
			//#pragma omp parallel for schedule(static, 16) //default(none) shared(distances)
			for (int i = 0; i < length; ++i)
				distances[i] = distance(list[i], p);

			MinIndex<float>(distances, length, N, nearestI, minDist);
		}
	}

template<int D, typename T, typename T3D>
class KdTreeSearch;


template<int D = 3, typename T = float, typename T3D = ::glm::vec3>
class KdTree {
	typedef KdTree<D,T,T3D> Tree;

	template <int D, typename T, typename T3D>
	friend class KdTreeSearch;

	class Node {
		friend class Tree;

		virtual bool IsLeaf() const = 0;
		virtual void Add(T3D const& p, Tree & parent) = 0;
	};

	class Stem : public Node {
		friend class Tree;
		friend class KdTreeSearch<D,T,T3D>;

		T splitValue;
		int splitAxis;

		Tree* left;
		Tree* right;


		Stem(Tree* left, Tree* right, T splitValue, int splitAxis)
			: left(left)
			, right(right)
			, splitValue(splitValue)
			, splitAxis(splitAxis)
		{}

		~Stem() {
			delete left;
			delete right;
		}

		inline bool IsLeft(T3D const& p) const { return p[splitAxis] < splitValue; }
		bool IsLeaf() const { return false; }

		void Add(T3D const& p, Tree & tree)
		{
			Add(p);
		}

		inline void Add(T3D const& p)
		{
			Tree* child = IsLeft(p) ? left : right;
			child->Add(p);
		}
	};

	class Leaf : public Node {
		friend class Tree;
		friend class KdTreeSearch<D,T,T3D>;

		// Bounds are only ever used when splitting
		T3D bounds[2];
		std::vector<T3D> bucket;

		Leaf(int bucketsize)
		{
			bounds[0] = bounds[1] = 0;

			bucket.reserve(bucketsize);
		}

		inline int GetSize() const {
			return bucket.size();
		}

		inline void InitBounds(T3D const& p)
		{
			bounds[0] = p;
			bounds[1] = p;
		}
		inline void UpdateBounds(T3D const& p) 
		{
			//other code in Leaf::Add() relies on bounds[1] being >= bounds[0] to save abs
			bounds[0] = ::glm::min(p, bounds[0]);
			bounds[1] = ::glm::max(p, bounds[1]);
		}

		bool IsLeaf() const { return true; }




		inline int GetSplitAxis()
		{
			T3D extend = bounds[1] - bounds[0];
			int axis = 0;
			T range = extend[0];
			for (int i = 1; i < D; ++i) {
				if (extend[i] > range) {
					axis = i;
					range = extend[i];
				}
			}
			if (extend[axis] < 0.001f)
			{
				return -1;
			}
			return axis;
		}

		virtual void Add(T3D const& p, Tree & tree)
		{
			bucket.empty() ?
				InitBounds(p) : UpdateBounds(p);
			bucket.push_back(p);

			//dont need to split
			if (GetSize() < tree.bucketsize)
			{
				return;
			}

			//we need to split
			ERROR_ASSERT(bucket.size() >= tree.bucketsize);

			//double bucket size if no extend in any direction
			//split axis is -1 if we have no extend in any direction
			int splitAxis = GetSplitAxis();
			T const EPSILON = 0.0000001;
			if (splitAxis < 0
				//relies on bounds[1] being >= bounds[0] to save the abs
				|| bounds[1][splitAxis] - bounds[0][splitAxis] <= EPSILON)
			{
				tree.bucketsize *= 2;
				return;
			}

			T splitValue = (bounds[0][splitAxis] + bounds[1][splitAxis]) / 2;
			ERROR_ASSERT(bucket.size() == tree.bucketsize);

			// Don't let the split value be the same as the upper value as
			// can happen due to rounding errors!
			//if (splitValue == bounds[1][splitAxis])
			//{
			//	splitValue = bounds[0][splitAxis];
			//}


			//@todo use this node as new leafnode for left child
			tree.stem = new Stem(new Tree(tree.bucketsize), new Tree(tree.bucketsize), splitValue, splitAxis);

			//move bucket contents into children
			for (int i = 0; i < bucket.size(); ++i) 
			{
				tree.stem->Add(bucket[i], tree);
			}

			//replace ourselves with new stem
			delete tree.leaf;
			tree.leaf = nullptr;
		}

	};


    int bucketsize;


	Leaf* leaf;
	Stem* stem;

public:

/**
 * Constructor for child nodes. Internal use only.
 */
	KdTree(int bucketsize)
		: bucketsize(bucketsize)
		, stem(NULL)
		, leaf(new Leaf(bucketsize))
		{
			ASSERT(bucketsize > 0);
			ASSERT(D > 0);
		}

	KdTree(Leaf* leaf, int bucketsize)
		: bucketsize(bucketsize)
		, stem(NULL)
		, leaf(leaf)
		{
		}


	~KdTree() { 
		if (leaf) { 
			delete leaf;
		}
		else 
		{
			delete stem;
		}
	}

	void Add(T3D const& p)
	{
		if (leaf)
		{
			leaf->Add(p, *this);
		}
		else
		{
			stem->Add(p);
		}

	}

	KdTreeSearch<D,T,T3D> GetSearch()
	{
		return KdTreeSearch<D,T,T3D>(*this);
	}


};

template <int D, typename T, typename T3D>
class KdTreeSearch {
	typedef KdTree<D,T,T3D> Tree;

	//only used for speedup
	std::vector<float> searchMem;
	Tree const& tree;

public:
	KdTreeSearch(Tree const& tree)
		: searchMem(tree.bucketsize)
		, tree(tree)
	{}

	T3D NearestNeighbour(T3D const& p)
	{
		T3D nearest;
		T minDist = FLT_MAX;

		std::vector<Tree const*> stack;
		stack.reserve(256);
		stack.push_back(&tree);
		while (!stack.empty())
		{
			Tree const* const tree = stack.back();
			stack.pop_back();

			if (tree->leaf) 
			{

				float d;
				ERROR_ASSERT(tree->leaf->bucket.size() <= tree->bucketsize);
				int const numPoints = tree->leaf->bucket.size();
				if (searchMem.size() < numPoints)
				{
					searchMem.resize(numPoints);
				}

				int i;
				helper::FindNearest<T3D>(&(tree->leaf->bucket[0]), numPoints, p, 1, &i, &d, &(searchMem[0]));
				T3D tmp = tree->leaf->bucket[i];

				if (d < minDist) {
					minDist = d;
					nearest = tmp;
				}
			}
			else
			{
				auto const* const stem = tree->stem;

				float splitDist = p[stem->splitAxis] - stem->splitValue;
				Tree* nearestChild;
				Tree* furthestChild;
				//IsLeft()
				if (splitDist < 0) {
					nearestChild = stem->left;
					furthestChild = stem->right;
					splitDist = -splitDist;
				}
				else
				{
					nearestChild = stem->right;
					furthestChild = stem->left;
				}

				//discard further child if distance of p to split plane is already larger than the best neighbour found so far
				if (splitDist < minDist) 
				{
					stack.push_back(furthestChild);
				}
				//nearest child is always inside
				stack.push_back(nearestChild);
			}
		}
		return nearest;

#if 0
		auto curNode = this;
		float dist = 0;
		float maxDist = FLT_MAX;
		while (dist < maxDist)
		{
			//down traversal
			while (!curNode.leaf) 
			{
				Stem const* const stem = tree->stem;

				float splitDist = p[stem->splitAxis] - stem->splitValue;
				Tree* nearestChild;
				Tree* furthestChild;
				//IsLeft()
				if (splitDist < 0) {
					curNode = stem->left;
				}
				else
				{
					curNode = stem->right;
				}

				//discard further child if distance of p to split plane is already larger than the best neighbour found so far
				if (splitDist < minDist) 
				{
					stack.push_back(furthestChild);
				}
				//nearest child is always inside
				stack.push_back(nearestChild);

			}
		}
#endif
	}

};
}

#endif // _SPATIAL_KDTREE_H_
