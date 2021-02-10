
#ifndef CONFECALC_ARRAY_H
#define CONFECALC_ARRAY_H


namespace osprey {

	// a very simple fixed-length array,
	// that can be allocated from either the Java or c++ side
	template<typename T>
	class Array {
		public:

			Array(int64_t size): size(size), things(new T[size]) {}
			Array(const Array<T> & other) = delete;
			~Array() {
				if (things != nullptr) {
					delete[] things;
				}
			}

			inline int64_t get_size() const {
				return size;
			}

			inline T & operator [] (int64_t i) {

				// just in case ...
				assert (i >= 0);
				assert (i < size);

				return pointer()[i];
			}

			inline const T & operator [] (int64_t i) const {

				// just in case ...
				assert (i >= 0);
				assert (i < size);

				return pointer()[i];
			}

			inline int64_t copy_from(const Array<T> & src, int64_t srci, int64_t count, int64_t dsti) {

				// just in case...
				assert(dsti >= 0);
				assert(dsti + count <= size);
				assert(srci >= 0);
				assert(srci + count <= src.size);

				std::copy(src.pointer() + srci, src.pointer() + srci + count, pointer() + dsti);

				return count;
			}

			inline int64_t copy_from(const Array<T> & src, int64_t dsti) {
				return copy_from(src, 0, src.get_size(), dsti);
			}

			inline int64_t copy_from(const Array<T> & src) {
				return copy_from(src, 0);
			}

			inline void truncate(int64_t smaller_size) {

				assert (smaller_size <= size);

				size = smaller_size;
			}

		private:

			int64_t size;
			T * things; // nullptr when created from java

			T * pointer() {
				if (things == nullptr) {
					// when created from java, the coords follow the class layout
					return reinterpret_cast<T *>(this + 1);
				} else {
					// when created from C++, the coords are allocated on the heap
					return things;
				}
			}

			const T * pointer() const {
				if (things == nullptr) {
					// when created from java, the coords follow the class layout
					return reinterpret_cast<const T *>(this + 1);
				} else {
					// when created from C++, the coords are allocated on the heap
					return things;
				}
			}
	};
	ASSERT_JAVA_COMPATIBLE(Array<int>, 16);
}


#endif //CONFECALC_ARRAY_H
