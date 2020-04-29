
#ifndef CONFECALC_ATOMS_H
#define CONFECALC_ATOMS_H


namespace osprey {

	template<typename T>
	static void write_atom(std::ostream & out, int i, const Real3<T> & atom) {
		out << fmt::make_indents << "[" << std::setw(5) << i << "] " << atom << std::endl;
	}

	template<typename T>
	std::ostream & operator << (std::ostream & out, const Array<Real3<T>> & atoms) {

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
