
#include <numeric>
namespace spatial {

	namespace helper {
		template <typename T>
		int MinIndex(T const* const seq, int const length, int N, int* minIndex, T* min) {
			for (int i = 0; i < N; ++i)
			{
				min[i] = ::std::numeric_limits<T>::max();
				minIndex[i] = -1;
			}

			for (int i = 0; i < length; ++i)
			{
				for (int j = 0; j < N; ++j)
				{
					if (seq[i] < min[j])
					{
						for (int h = N - 1; h > j; --h)
						{
							min[h] = min[h - 1];
							minIndex[h] = minIndex[h - 1];
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
			assert(length > 0);
			//#pragma omp parallel for schedule(static, 16) //default(none) shared(distances)
			for (int i = 0; i < length; ++i)
				distances[i] = distance(list[i], p);

			MinIndex<float>(distances, length, N, nearestI, minDist);
		}
	}

	template<int D, typename T, typename T3D>
	class KdTree {
		typedef KdTree<D, T, T3D> Tree;
		typedef Tree Self;

		struct Index
		{
		private:
			int i;
			bool leaf;

		public:
			Index()
			{
			}
		private:
			Index(int i, bool leaf)
				: i(i)
				, leaf(leaf)
			{
			}
		public:

			static inline Index make_leaf(int i)
			{
				return Index(-(i + 1), true);
			}
			static inline Index make_stem(int i)
			{
				return Index(i, false);
			}
			inline int from_leaf() const
			{
				assert(i < 0);
				assert(leaf == i < 0);
				return -i - 1;
			}

			inline bool is_leaf() const
			{
				assert(leaf == i < 0);
				return i < 0;
			}

			inline int from_stem() const
			{
				assert(i >= 0);
				assert(!leaf);
				return i;
			}

		};

		struct Leaf
		{
			::std::vector<T3D> bucket;
			Leaf(int bucketsize)
			{
				bucket.reserve(bucketsize);
			}

			inline int size() const
			{
				return int(bucket.size());
			}

			inline void push_back(T3D const& e)
			{
				bucket.push_back(e);
			}

			inline T3D const& operator[](int i) const
			{
				return const_cast<Leaf*>(this)->operator[](i);
			}

			inline T3D& operator[](int i)
			{
				return bucket[i];
			}

			inline auto begin() const -> decltype(bucket.begin())    { return bucket.begin(); }
			inline auto end()   const -> decltype(bucket.end  ())    { return bucket.end(); }
			inline auto begin()       -> decltype(bucket.begin())    { return bucket.begin(); }
			inline auto end()         -> decltype(bucket.end  ())    { return bucket.end(); }

			Leaf(Leaf&& o)
				: bucket(::std::move(bucket))
			{
			}

			Leaf& operator=(Leaf&& o)
			{
				::std::swap(bucket, o.bucket);
			}
		};

		struct Stem
		{
			int splitAxis;
			T splitValue;

			Index children[2];


		private:
			Stem()
				: children({ 0, 0 })
				, splitValue(0)
				, splitAxis(0)
			{
			}
		public:

			/** split a leaf node */
			Stem(int const splitAxis
				, T const splitValue
				, Leaf& leaf
				, Index const ileaf
				, int const bucketsize
				, ::std::vector<Stem>& stems
				, ::std::vector<Leaf>& leafs
				)
				: splitAxis(splitAxis)
				, splitValue(splitValue)
			{
				assert(ileaf.is_leaf());
				assert(leaf.size() % bucketsize == 1);


				
				//use old leaf as left child
				children[0] = ileaf;
				children[1] = Index::make_leaf(leafs.size());
				leafs.emplace_back(bucketsize);

				auto const old = ::std::move(leaf);


				//move bucket contents into children
				for (auto& p: old)
				{
					auto const c = int(is_upper(p));
					auto& child = children[c];
					leafs[child.from_leaf()].push_back(::std::move(p));
				}

				assert(leafs.size() == stems.size() + 1);
			}

			~Stem() 
			{
			}

			inline Index left() const
			{
				return children[0];
			}
			inline Index right() const
			{
				return children[1];
			}



			inline bool is_upper(T3D const& p) const { return p[splitAxis] > splitValue; }
			bool IsLeaf() const { return false; }

			void Add(T3D const& p, int bucketsize, ::std::vector<Stem>& stems, ::std::vector<Leaf>& leafs)
			{
				auto const c = int(is_upper(p));
				auto& child = children[c];

				add(p, child, bucketsize, stems, leafs);

			}

			static inline int GetSplitAxis(T3D const& bmin, T3D const& bmax, T& range)
			{
				auto extend = bmax - bmin;
				int axis = 0;
				range = extend[0];
				for (int i = 1; i < D; ++i) {
					if (extend[i] > range) {
						axis = i;
						range = extend[i];
					}
				}
				return axis;
			}

			static inline void add(
				T3D const& p
				, Index& child
				, int bucketsize
				, ::std::vector<Stem>& stems
				, ::std::vector<Leaf>& leafs)
			{
				if (child.is_leaf())
				{
					auto& leaf = leafs[child.from_leaf()];
					auto& bucket = leaf;
					//we already increased bucketsize
					while (bucket.size() > bucketsize)
					{
						//so account for that
						bucketsize *= 2;
					}
					bucket.push_back(p);

					auto needs_split = (bucket.size() > bucketsize);
					if (needs_split)
					{
						//we need to split
						assert(bucket.size() >= bucketsize);

						T range;

						T3D bmin = bucket[0];
						T3D bmax = bucket[0];
						for (int i = 1; i < bucket.size(); ++i)
						{
							auto const& x = bucket[i];
							bmin = min(x, bmin);
							bmax = max(x, bmax);
						}

						auto const splitAxis = GetSplitAxis(bmin, bmax, range);
						//only split if points extend in some direction at least EPSILON
						if (range > ::std::numeric_limits<T>::epsilon())
						{
							auto const ileaf = child;
							child = Index::make_stem(stems.size());

							auto splitValue = bmin[splitAxis] + range / T(2);
							assert(splitValue > bmin[splitAxis]);
							assert(splitValue < bmax[splitAxis]);

							stems.emplace_back(splitAxis, splitValue, leaf, ileaf, bucketsize, stems, leafs); //make new stem
						}
						//otherwise just keep this leaf
					}
				}
				else
				{
					stems[child.from_stem()].Add(p, bucketsize, stems, leafs);
				}
			}

		};

		

		int bucketsize;
		::std::vector<Stem> stems;
		::std::vector<Leaf> leafs;
		Index child;

	public:

		KdTree(int bucketsize = 16, int default_nodes = 256)
			: bucketsize(bucketsize)
		{
			assert(bucketsize > 0);
			assert(D > 0);

			stems.reserve(default_nodes);
			leafs.reserve(default_nodes);

			child = Index::make_leaf(leafs.size());
			leafs.emplace_back(bucketsize);
		}


		virtual ~KdTree() 
		{
		}

		void Add(T3D const& p)
		{
			assert(leafs.size() == stems.size() + 1);
			Stem::add(p, child, bucketsize, stems, leafs);
			assert(leafs.size() == stems.size() + 1);
		}

		Self(Self&& o)
			: bucketsize(::std::move( o.bucketsize))
			, leafs(::std::move(o.leafs))
			, stems(::std::move(o.stems))
			, child(::std::move(o.child))
		{
		}

		Self& operator=(Self&& o)
		{
			::std::swap(bucketsize, o.bucketsize);
			::std::swap(leafs, o.leafs);
			::std::swap(stems, o.stems);
			::std::swap(child, o.child);
		}

	private:
		Self(Self const& o);
		Self& operator=(Self const& o);



	public:
		class Search {
			typedef Search Self;

			//only used for speedup
			mutable std::vector<float> searchMem;
			Tree const* tree;

		public:
			Search(Tree const* tree)
				: searchMem(tree->bucketsize)
				, tree(tree)
			{}

			//T3D NearestNeighbour(T3D const& p) const
			//{
			//	T3D r;
			//	return
			//}

			T3D NearestNeighbour(T3D const& p, T minDist = ::std::numeric_limits<T>::max()) const
			{
				T3D nearest;

				if (tree->child.is_leaf())
				{
					int i = search_leaf(p, tree->leafs[tree->child.from_leaf()], minDist, nearest);
					return nearest;
				}

				std::vector<Stem const*> stack;
				stack.reserve(256);

				stack.push_back(&tree->stems[tree->child.from_stem()]);
				while (!stack.empty())
				{
					auto const& stem = *stack.back();
					stack.pop_back();

					{
						float splitDist = p[stem.splitAxis] - stem.splitValue;
						auto nearestChild  = stem.left();
						auto furthestChild = stem.right();
						//IsLeft()
						if (splitDist < 0) {
							::std::swap(nearestChild, furthestChild);
							splitDist = -splitDist;
						}


						if (furthestChild.is_leaf())
						{
							int i = search_leaf(p, tree->leafs[furthestChild.from_leaf()], minDist, nearest);
						}
						else
						{
							//discard further child if distance of p to split plane is already larger than the best neighbour found so far
							//if (splitDist < minDist)
							{
								stack.push_back(&tree->stems[furthestChild.from_stem()]);
							}
						}

						if (nearestChild.is_leaf())
						{
							int i = search_leaf(p, tree->leafs[nearestChild.from_leaf()], minDist, nearest);
						}
						else
						{
							//nearest child is always inside
							stack.push_back(&tree->stems[nearestChild.from_stem()]);
						}
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

			static inline int search_leaf(T3D const& p, Leaf const& leaf, T& minDist, T3D& nearest)
			{
				//assert(leaf.size() <= tree->bucketsize);
				int nearestI = -1;

				int i = 0;
				for (auto& o : leaf)
				{
					auto dist = distance(p, o);

					if (dist < minDist) {
						minDist = dist;
						nearestI = i;
					}

					++i;
				}
				if (nearestI >= 0)
				{
					nearest = leaf[nearestI];
				}
				return nearestI;
			}

			Self(Self&& o)
				: tree(o.tree)
				, searchMem(o.searchMem)
			{
			}

			Self& operator=(Self&& o)
			{
				::std::swap(tree, o.tree);
				::std::swap(searchMem, o.searchMem);
			}

		};

		Search GetSearch()
		{
			return Search(this);
		}


	};


}