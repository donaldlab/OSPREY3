
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
			Array() = delete;

			__host__ __device__
			Array(const Array<T> & other) = delete;

			// host/device annotions not needed here for some reason
			~Array() = default;

			// set the size when the Array is in shared or global memory
			__device__
			inline void init(int64_t _size) {
				if (threadIdx.x == 0) {
					size = _size;
				}
				__syncthreads();
			}

			// for allocating on a single thread
			__host__ __device__
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
				return cuda::pad_to_alignment<16>(sizeof(Array<T>) + size*sizeof(T));
			}

			// get the total allocated size of the array, in bytes
			__host__ __device__
			inline int64_t get_bytes() const {
				return get_bytes(size);
			}

			__host__ __device__
			inline T & operator [] (int64_t i) {

				// just in case ...
				assert (i >= 0);
				assert (i < size);

				return items()[i];
			}

			__host__ __device__
			inline const T & operator [] (int64_t i) const {

				// just in case ...
				assert (i >= 0);
				assert (i < size);

				return items()[i];
			}

			__host__ __device__
			inline T * items() {

				// make sure we're 16-byte aligned
				assert (cuda::is_aligned<16>(this));

				// the coords follow the class layout
				return reinterpret_cast<T *>(this + 1);
			}

			__host__ __device__
			inline const T * items() const {

				// make sure we're 16-byte aligned
				assert (cuda::is_aligned<16>(this));

				// the coords follow the class layout
				return reinterpret_cast<const T *>(this + 1);
			}

			__host__
			inline int64_t copy_from_host(const Array<T> & src, int64_t srci, int64_t count, int64_t dsti) {

				// just in case...
				assert (dsti >= 0);
				assert (dsti + count <= size);
				assert (srci >= 0);
				assert (srci + count <= src.size);

				std::copy(src.items() + srci, src.pointer() + srci + count, items() + dsti);

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
			inline int64_t copy_from_device(const Array<T> & src, int64_t srci, int64_t count, int64_t dsti) {

				// just in case...
				assert (dsti >= 0);
				assert (dsti + count <= size);
				assert (srci >= 0);
				assert (srci + count <= src.size);

				for (uint i=threadIdx.x; i<count; i += blockDim.x) {
					items()[dsti + i] = src.items()[srci + i];
				}
				__syncthreads();

				return count;
			}

			__device__
			inline int64_t copy_from_device(const Array<T> & src, int64_t dsti) {
				return copy_from_device(src, 0, src.get_size(), dsti);
			}

			__device__
			inline int64_t copy_from_device(const Array<T> & src) {
				return copy_from_device(src, 0);
			}

			__device__
			inline int64_t fill_device(int64_t dsti, int64_t count, const T & val) {

				// just in case...
				assert (dsti >= 0);
				assert (dsti + count <= size);

				for (uint i=threadIdx.x; i<count; i += blockDim.x) {
					items()[dsti + i] = val;
				}
				__syncthreads();

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
	};
	ASSERT_JAVA_COMPATIBLE(Array<int>, 16);
	// NOTE: the array header *must* be a multiple of 16 bytes for Real3<T> alignment to be correct
}


#endif //CONFECALC_ARRAY_H
