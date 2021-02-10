
#ifndef CONFECALC_FORMATS_H
#define CONFECALC_FORMATS_H


namespace osprey::fmt {

	// some stream magic so callers can set formatting params for real values
	struct RealFormat {
		int width;
		int precision;
	};

	inline RealFormat set_real(int width, int precision) {
		return { width, precision };
	}

	static int real_width_index = std::ios_base::xalloc();

	inline long & real_width(std::ostream & out) {
		return out.iword(real_width_index);
	}

	std::ostream & operator << (std::ostream & out, RealFormat fmt) {
		real_width(out) = fmt.width;
		out << std::fixed << std::setprecision(fmt.precision);
		return out;
	}


	// more stream magic so callers can set indentation
	struct IndentsFormat {
		int n;
	};

	inline IndentsFormat set_indents(int n) {
		return { n };
	}

	static int indents_index = std::ios_base::xalloc();

	inline long & indents(std::ostream & out) {
		return out.iword(indents_index);
	}

	std::ostream & operator << (std::ostream & out, IndentsFormat fmt) {
		indents(out) = fmt.n;
		return out;
	}

	struct MakeIndents {};
	static MakeIndents make_indents;
	std::ostream & operator << (std::ostream & out, MakeIndents ignored) {
		for (int i=0; i<indents(out); i++) {
			out << "\t";
		}
		return out;
	}
}


#endif //CONFECALC_FORMATS_H
