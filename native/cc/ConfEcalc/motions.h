
#ifndef CONFECALC_MOTIONS_H
#define CONFECALC_MOTIONS_H


namespace osprey::motions {

	template<typename T>
	struct alignas(8) Dihedral {

		const static int32_t id = 0;

		T min_radians;
		T max_radians;
		int32_t a; // atom indices
		int32_t b;
		int32_t c;
		int32_t d;
		int32_t num_rotated;
		int32_t modified_posi;

		// only created by java side
		Dihedral() = delete;

		inline const int32_t get_rotated(int i) const {

			// just in case ...
			assert (i >= 0);
			assert (i < num_rotated);

			return reinterpret_cast<const int32_t *>(this + 1)[i];
		}
	};
	ASSERT_JAVA_COMPATIBLE_REALS(Dihedral, 32, 40);

	template<typename T>
	struct alignas(8) TranslationRotation {

		const static int32_t id = 1;

		T max_distance;
		T max_radians;
		Real3<T> centroid;
		int32_t num_atoms;
		int32_t num_modified_pos;

		// only created by java side
		TranslationRotation() = delete;

		inline const int32_t get_atomi(int i) const {

			// just in case ...
			assert (i >= 0);
			assert (i < num_atoms);

			// atom indices are right after the struct
			return reinterpret_cast<const int32_t *>(this + 1)[i];
		}

		inline const int32_t get_modified_posi(int i) const {

			// just in case ...
			assert (i >= 0);
			assert (i > num_modified_pos);

			// modified position indices are right after the atom indices
			return reinterpret_cast<const int32_t *>(this + 1)[num_atoms + i];
		}
	};
	ASSERT_JAVA_COMPATIBLE_REALS(TranslationRotation, 32, 48);
}


#endif //CONFECALC_MOTIONS_H
