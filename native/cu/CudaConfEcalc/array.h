
#ifndef CONFECALC_ARRAY_H
#define CONFECALC_ARRAY_H


namespace osprey {

	// a very simple fixed-length array,
	// that can be allocated from either the Java or c++ side
	template<typename T>
	class Array {
		public:

			// only created via malloc-style allocations
			__host__ __device__
			Array(int64_t size) = delete;

			__host__ __device__
			Array(const Array<T> & other) = delete;

			// host/device annotions not needed here for some reason
			~Array() = default;

			// set the size when the Array is in shared or global memory
			__device__
			inline void init(int64_t _size, cg::thread_group threads) {
				if (threads.thread_rank() == 0) {
					size = _size;
				}
				threads.sync();
			}

			// for allocating on a single thread
			__device__
			static inline Array<T> * make(int64_t size) {
				auto out = reinterpret_cast<Array<T> *>(std::malloc(get_bytes(size)));
				out->size = size;
				return out;
			}

			// get the number of items in the array
			__host__ __device__
			inline int64_t get_size() const {
				return size;
			}

			__host__ __device__
			static inline int64_t get_bytes(int64_t size) {
				return cuda::pad_to_alignment(sizeof(Array<T>) + size*sizeof(T), 8);
			}

			// get the total allocated size of the array, in bytes
			__host__ __device__
			inline int64_t get_bytes() const {
				return get_bytes(size);
			}

			__host__ __device__
			inline T & operator[] (int64_t i) {

				// just in case ...
				assert (i >= 0);
				assert (i < size);

				return pointer()[i];
			}

			__host__ __device__
			inline const T & operator [] (int64_t i) const {

				// just in case ...
				assert (i >= 0);
				assert (i < size);

				return pointer()[i];
			}

			__host__
			inline int64_t copy_from_host(const Array<T> & src, int64_t srci, int64_t count, int64_t dsti) {

				// just in case...
				assert (dsti >= 0);
				assert (dsti + count <= size);
				assert (srci >= 0);
				assert (srci + count <= src.size);

				std::copy(src.pointer() + srci, src.pointer() + srci + count, pointer() + dsti);

				return count;
			}

			__host__
			inline int64_t copy_from_host(const Array<T> & src, int64_t dsti) {
				return copy_from_host(src, 0, src.get_size(), dsti);
			}

			__host__
			inline int64_t copy_from_host(const Array<T> & src) {
				return copy_from_host(src, 0);
			}

			__device__
			inline int64_t copy_from_device(const Array<T> & src, int64_t srci, int64_t count, int64_t dsti, cg::thread_group threads) {

				// just in case...
				assert (dsti >= 0);
				assert (dsti + count <= size);
				assert (srci >= 0);
				assert (srci + count <= src.size);

				for (int i=threads.thread_rank(); i<count; i += threads.size()) {
					operator[](dsti + i) = src[srci + i];
				}
				threads.sync();

				return count;
			}

			__device__
			inline int64_t copy_from_device(const Array<T> & src, int64_t dsti, cg::thread_group threads) {
				return copy_from_device(src, 0, src.get_size(), dsti, threads);
			}

			__device__
			inline int64_t copy_from_device(const Array<T> & src, cg::thread_group threads) {
				return copy_from_device(src, 0, threads);
			}

			__device__
			inline int64_t fill_device(int64_t dsti, int64_t count, const T & val, cg::thread_group threads) {

				// just in case...
				assert (dsti >= 0);
				assert (dsti + count <= size);

				for (int i=threads.thread_rank(); i<count; i += threads.size()) {
					operator[](dsti + i) = val;
				}
				threads.sync();

				return count;
			}

			__host__ __device__
			inline void truncate(int64_t smaller_size) {

				assert (smaller_size <= size);

				size = smaller_size;
			}

		private:

			int64_t size;
			int64_t pad; // need to pad to 16 bytes

			__host__ __device__
			inline T * pointer() {
				// the coords follow the class layout
				return reinterpret_cast<T *>(this + 1);
			}

			__host__ __device__
			inline const T * pointer() const {
				// the coords follow the class layout
				return reinterpret_cast<const T *>(this + 1);
			}
	};
	ASSERT_JAVA_COMPATIBLE(Array<int>, 16);
	// NOTE: the array header *must* be a multiple of 16 bytes for float3 alignment to be correct


	// the compiler thinks float3 is 12 bytes, and should be 4-byte aligned
	// it doesn't know that we're actually pretending float3 is 16 bytes, and should be 16-byte aligned
	// so add some specializations to override the default sizes/alignments
	template<>
	__host__
	inline int64_t Array<float3>::get_bytes(int64_t size) {
		return sizeof(Array<float3>) + size*Real3Map<float32_t>::size;
	}

	template<>
	__host__ __device__
	inline float3 & Array<float3>::operator[] (int64_t i) {

		// just in case ...
		assert (i >= 0);
		assert (i < size);

		auto p = reinterpret_cast<int8_t *>(pointer());
		return *reinterpret_cast<float3 *>(p + Real3Map<float32_t>::size*i);
	}

	template<>
	__host__ __device__
	inline const float3 & Array<float3>::operator[] (int64_t i) const {

		// just in case ...
		assert (i >= 0);
		assert (i < size);

		auto p = reinterpret_cast<const int8_t *>(pointer());
		return *reinterpret_cast<const float3 *>(p + Real3Map<float32_t>::size*i);
	}
}


#endif //CONFECALC_ARRAY_H
