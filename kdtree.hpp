
#include <numeric>
namespace spatial {

	namespace detail {
		/** 
			shift a range back by 1, and cut off the last element if necessary 
			returns the new end
		*/
		template <class It>
		It move_backward_cutoff(It const& begin, It const& end, It const& end_allocated)
		{
			auto const& end_ref = (end != end_allocated) ? end + 1 : end_allocated;
			::std::move_backward(begin, end_ref - 1, end_ref);
			return end_ref;
		}

		template<int D, class T, class T3D, class Data, class SearchResult>
		struct Leaf
		{
			//T3D bucket;
			::std::vector<T3D> bucket;
			::std::vector<Data> data_bucket;

			Leaf(int bucketsize)
			{
				bucket.reserve(bucketsize);
				data_bucket.reserve(bucketsize);
			}

			inline int size() const
			{
				assert(bucket.size() == data_bucket.size());
				return int(bucket.size());
			}

			inline void push_back(T3D const& e, Data const& d)
			{
				bucket.push_back(e);
				data_bucket.push_back(d);
			}

			inline T3D const& operator[](int i) const
			{
				return const_cast<Leaf*>(this)->operator[](i);
			}

			inline T3D& operator[](int i)
			{
				return bucket[i];
			}

			inline Data const& data(int i) const
			{
				return const_cast<Leaf*>(this)->data(i);
			}
			inline Data& data(int i)
			{
				return data_bucket[i];
			}

			inline auto begin() const -> decltype(bucket.begin())    { return bucket.begin(); }
			inline auto end()   const -> decltype(bucket.end())    { return bucket.end(); }
			inline auto begin()       -> decltype(bucket.begin())    { return bucket.begin(); }
			inline auto end()         -> decltype(bucket.end())    { return bucket.end(); }

			Leaf(Leaf&& o)
				: bucket(::std::move(o.bucket))
				, data_bucket(::std::move(o.data_bucket))
			{
			}

			Leaf& operator=(Leaf&& o)
			{
				::std::swap(bucket, o.bucket);
				::std::swap(data_bucket, o.data_bucket);
			}

			/**
			DistanceFun receives two T3D indices, and the maximum valid distance
			*/
			template <class DistIt, class MakeData, class DistanceFun>
			/*inline */ void search(
				T3D const& p
				, T& minDist
				, DistIt const& begin
				, DistIt      & end
				, DistIt const& end_allocated
				, MakeData make_data
				, DistanceFun distance
				) const
			{
				auto const num_allocated = ::std::distance(begin, end_allocated);
				auto const& leaf = *this;

				for (int i = 0; i < leaf.size(); ++i)
				{
					auto const& o = leaf[i];

					float dist = distance(p, o, minDist);

					if (dist <= minDist) {
						//search for lowbound in valid matches
#if 0 //1.61m
						auto const lower = ::std::upper_bound(begin, end, dist);
#elif 1 //1.55m
						//linear search
						auto lower = begin;
						for (; lower != end && *lower < dist; ++lower)
						{
						}
#else //1.54m
						//linear search with one step of binary search
						auto const n = ::std::distance(begin, end);
						auto const half = (n / 2);
						auto lower = *(begin + half) <= dist ? begin + half : begin;
						for (; lower != end && *lower <= dist; ++lower)
						{
						}
#endif
						if (lower != end_allocated)
						{
							auto const num_lower = ::std::distance(begin, lower);
							end = move_backward_cutoff(lower, end, end_allocated);
							*lower = dist;
							
							make_data(dist, num_lower, o, data_bucket[i]);
						}
						//if (end == end_allocated)
						//{
						//	minDist = *(end_allocated - 1);
						//}
					}
				}
			}
		};

		
	}

	template<int D, typename T, typename T3D, typename Data, typename SearchResult>
	class KdTree {
		typedef KdTree<D, T, T3D, Data, SearchResult> Tree;
		typedef Tree Self;

	public:
		typedef detail::Leaf<D, T, T3D, Data, SearchResult> Leaf;

		struct Index
		{
		private:
			int i;

		public:
			Index()
			{
			}
		private:
			Index(int i)
				: i(i)
			{
			}
		public:

			static inline Index make_leaf(int i)
			{
				return Index(-(i + 1));
			}
			static inline Index make_stem(int i)
			{
				return Index(i);
			}
			inline int from_leaf() const
			{
				assert(i < 0);
				return -i - 1;
			}

			inline bool is_leaf() const
			{
				return i < 0;
			}

			inline int from_stem() const
			{
				assert(i >= 0);
				return i;
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
				, Index const ioldleaf
				, int const bucketsize
				, ::std::vector<Stem>& stems
				, ::std::vector<Leaf>& leafs
				)
				: splitAxis(splitAxis)
				, splitValue(splitValue)
			{
				assert(ioldleaf.is_leaf());
				//use old leaf as lower child
				children[0] = ioldleaf;
				//create new leaf
				children[1] = Index::make_leaf(leafs.size());
				leafs.emplace_back(bucketsize);
				//need to take new pointers
				auto& lessleaf    = leafs[children[0].from_leaf()];
				auto& greaterleaf = leafs[children[1].from_leaf()];
				assert(lessleaf.size() % bucketsize == 1);
				assert(greaterleaf.size() == 0);

				auto const old = ::std::move(lessleaf);
				assert(old.size() % bucketsize == 1);
				assert(lessleaf.size() == 0);
				assert(greaterleaf.size() == 0);


				//move bucket contents into children
				for (int i = 0; i < old.size(); ++i)
				{
					auto& p = old[i];
					auto const less = !is_upper(p);
					(less ? lessleaf : greaterleaf).push_back(::std::move(p), old.data(i));
				}
				assert(lessleaf.size());
				assert(greaterleaf.size());
			}

			~Stem() 
			{
			}

			inline Index lower() const
			{
				return children[0];
			}
			inline Index upper() const
			{
				return children[1];
			}


			inline Index& get_index(T3D const& p)
			{
				return children[is_upper(p)];
			}

			inline bool is_upper(T3D const& p) const { return p[splitAxis] > splitValue; }

			void Add(
				T3D const& p
				, Data const& d
				, int bucketsize
				, ::std::vector<Stem>& stems
				, ::std::vector<Leaf>& leafs
				)
			{
				add(p, d, get_index(p), bucketsize, stems, leafs);
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
				, Data const& d
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
					bucket.push_back(p, d);

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

							stems.emplace_back(splitAxis, splitValue, ileaf, bucketsize, stems, leafs); //make new stem
							assert(leafs.size() == stems.size() + 1);

						}
						//otherwise just keep this leaf
					}
				}
				else
				{
					stems[child.from_stem()].Add(p, d, bucketsize, stems, leafs);
				}
			}

		};

		

		int bucketsize;
		::std::vector<Stem> stems;
		::std::vector<Leaf> leafs;
		Index child;
	public:
		static int const CACHELINE_SIZE = 64 * sizeof(uint8_t);
		//8: 31.0550s
		//4: 31.9290s
		//12: 35.2570s
		//16: 38.2130s
		KdTree(int bucketsize = CACHELINE_SIZE * 8 / sizeof(T3D) 
			, int default_nodes = 64)
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

		void Add(T3D const& p, Data const& d)
		{
			assert(leafs.size() == stems.size() + 1);
			Stem::add(p, d, child, bucketsize, stems, leafs);
			assert(leafs.size() == stems.size() + 1);
		}

		KdTree(Self&& o)
			: bucketsize(::std::move( o.bucketsize))
			, stems(::std::move(o.stems))
			, leafs(::std::move(o.leafs))
			, child(::std::move(o.child))
		{
		}

		Self& operator=(Self&& o)
		{
			::std::swap(bucketsize, o.bucketsize);
			::std::swap(leafs, o.leafs);
			::std::swap(stems, o.stems);
			::std::swap(child, o.child);
			return *this;
		}

	private:
		KdTree(Self const& o);
		Self& operator=(Self const& o);



	public:
		class Search {
			typedef Search Self;
			static int const search_stack_size = 64;
			Tree const* tree;

		public:
			Search(Tree const* tree)
				: tree(tree)
			{}

			template <class DistIt, class MakeData, class DistanceFun>
			void NearestNeighbour(T3D const& p
				, DistIt const& begin
				, DistIt      & end
				, DistIt const& end_allocated
				, MakeData make_data
				, T const max_distance = ::std::numeric_limits<T>::max()
				, DistanceFun distance = [&] (T3D const& p, T3D const& o, T minDist) {
						float dist = 0;
						for (int j = 0; j < T3D::dim && dist <= minDist; ++j)
						{
							auto const a = p[j] - o[j];
							dist += a * a;
						}
						return dist;
					}
				) const
			{
				auto minDist = max_distance;

				auto const get_leaf = [this](Index const& i) -> Leaf const& {
					return tree->leafs[i.from_leaf()]; 
				};

				if (tree->child.is_leaf())
				{
					get_leaf(tree->child).search(p, minDist, begin, end, end_allocated, make_data, distance);
					return;
				}

				std::vector<Stem const*> stack;
				stack.reserve(search_stack_size);

				stack.push_back(&tree->stems[tree->child.from_stem()]);
				while (!stack.empty())
				{
					auto const stem = *stack.back();
					stack.pop_back();

					{
						//splitdist > 0 -> we're higher than splitplane
						auto nearestChild  = stem.upper();
						auto furthestChild = stem.lower();
						float splitDist = p[stem.splitAxis] - stem.splitValue;
						if (splitDist < 0) {
							::std::swap(nearestChild, furthestChild);
							splitDist = -splitDist;
						}


						if (furthestChild.is_leaf())
						{
							get_leaf(furthestChild).search(p, minDist, begin, end, end_allocated, make_data, distance);
						}
						else
						{
							//discard further child if distance of p to split plane is already larger than the best neighbour found so far
							if (splitDist <= minDist)
							{
								stack.push_back(&tree->stems[furthestChild.from_stem()]);
							}
						}

						if (nearestChild.is_leaf())
						{
							get_leaf(nearestChild).search(p, minDist, begin, end, end_allocated, make_data, distance);
						}
						else
						{
							//nearest child is always inside
							stack.push_back(&tree->stems[nearestChild.from_stem()]);
						}
					}
				}
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


			Search(Self&& o)
				: tree(::std::move(o.tree))
			{
			}

			Self& operator=(Self&& o)
			{
				::std::swap(tree, o.tree);
				return *this;
			}

		};

		Search GetSearch()
		{
			return Search(this);
		}


	};


}
