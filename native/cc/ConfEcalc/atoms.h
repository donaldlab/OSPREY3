
#ifndef CONFECALC_ATOMS_H
#define CONFECALC_ATOMS_H


namespace osprey {

	template<typename T>
	class Atoms {

		public:
			explicit Atoms(int32_t _num_atoms)
				: num_atoms(_num_atoms), coords(new Real3<T>[_num_atoms]) {
			}

			~Atoms() {
				if (coords != nullptr) {
					delete[] coords;
				}
			}

			inline int32_t size() const {
				return num_atoms;
			}

			inline Real3<T> & operator [] (int i) {
				return coords_pointer()[i];
			}
			inline const Real3<T> & operator [] (int i) const {
				return coords_pointer()[i];
			}

			inline int64_t copy_from(const Atoms<T> & src, int64_t srci, int64_t count, int64_t dsti) {

				// just in case...
				assert(dsti >= 0);
				assert(dsti + count <= num_atoms);
				assert(srci >= 0);
				assert(srci + count <= src.num_atoms);

				std::copy(src.coords_pointer() + srci, src.coords_pointer() + srci + count, coords_pointer() + dsti);

				return count;
			}

			inline int64_t copy_from(const Atoms<T> & src, int64_t dsti) {
				return copy_from(src, 0, src.num_atoms, dsti);
			}

			inline int64_t copy_from(const Atoms<T> & src) {
				return copy_from(src, 0);
			}

		private:
			const int32_t num_atoms;
			// 4 bytes pad
			Real3<T> * coords; // nullptr when created from java

			Real3<T> * coords_pointer() {
				if (coords == nullptr) {
					// when created from java, the coords follow the class layout
					return reinterpret_cast<Real3<T> *>(this + 1);
				} else {
					// when created from C++, the coords are allocated on the heap
					return coords;
				}
			}

			const Real3<T> * coords_pointer() const {
				if (coords == nullptr) {
					// when created from java, the coords follow the class layout
					return reinterpret_cast<const Real3<T> *>(this + 1);
				} else {
					// when created from C++, the coords are allocated on the heap
					return coords;
				}
			}
	};
	ASSERT_JAVA_COMPATIBLE_REALS(Atoms, 16, 16);

	template<typename T>
	static void write_atom(std::ostream & out, int i, const Real3<T> & atom) {
		out << fmt::make_indents << "[" << std::setw(5) << i << "] " << atom << std::endl;
	}

	template<typename T>
	std::ostream & operator << (std::ostream & out, const Atoms<T> & atoms) {

		out << fmt::set_real(12, 6);

		if (atoms.size() <= 10) {

			// not that many atoms, just write all of them
			for (int i=0; i<atoms.size(); i++) {
				write_atom(out, i, atoms[i]);
			}

		} else {

			// lots of atoms, only write a few at the top and bottom
			for (int i=0; i<3; i++) {
				write_atom(out, i, atoms[i]);
			}
			out << fmt::make_indents << "..." << std::endl;
			for (int i=atoms.size()-3; i<atoms.size(); i++) {
				write_atom(out, i, atoms[i]);
			}
		}

		return out;
	}
}


#endif //CONFECALC_ATOMS_H
