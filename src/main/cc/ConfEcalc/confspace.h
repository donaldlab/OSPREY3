
#ifndef CONFECALC_CONFSPACE_H
#define CONFECALC_CONFSPACE_H


namespace osprey {

	const int32_t StaticPos = -1;

	template<typename T>
	struct Conf {
		int64_t atom_coords_offset;
		int64_t atom_molis_offset;
		int32_t frag_index;
		// 4 bytes pad if T = float32_t
		T internal_energy;
		int64_t num_motions;
		int64_t motions_offset;
	};
	ASSERT_JAVA_COMPATIBLE_REALS(Conf, 40, 48);

	struct alignas(8) Pos {
		int32_t num_confs;
		int32_t max_num_atoms;
		int32_t num_frags;
		// 4 bytes pad
	};
	ASSERT_JAVA_COMPATIBLE(Pos, 16);

	template<typename T>
	class ConfSpace {
		public:
			ConfSpace() = delete; // created only on the Java side

			inline const Pos & get_pos(int posi) const {
				auto offsets = reinterpret_cast<const int64_t *>(offset(positions_offset));
				return *reinterpret_cast<const Pos *>(offset(offsets[posi]));
			}

			inline const Conf<T> & get_conf(const Pos & pos, int confi) const {
				// the java side puts the conf offsets right after the Pos struct
				auto conf_offsets = reinterpret_cast<const int64_t *>(&pos + 1);
				return *reinterpret_cast<const Conf<T> *>(offset(conf_offsets[confi]));
			}

			inline const Array<Real3<T>> & get_static_atom_coords() const {
				return *reinterpret_cast<const Array<Real3<T>> *>(offset(static_atom_coords_offset));
			}

			inline const Array<int32_t> & get_static_atom_molis() const {
				return *reinterpret_cast<const Array<int32_t> *>(offset(static_atom_molis_offset));
			}

			inline const Array<Real3<T>> & get_conf_atom_coords(const Conf<T> & conf) const {
				return *reinterpret_cast<const Array<Real3<T>> *>(offset(conf.atom_coords_offset));
			}

			inline const Array<int32_t> & get_conf_atom_molis(const Conf<T> & conf) const {
				return *reinterpret_cast<const Array<int32_t> *>(offset(conf.atom_molis_offset));
			}

			inline const void * get_params() const {
				return reinterpret_cast<const void *>(offset(params_offset));
			}

			inline int64_t index_static_static() const {
				return 0;
			}

			inline int64_t index_static_pos(int posi1) const {
				return 1 + posi1;
			}

			inline int64_t index_pos(int posi1) const {
				return 1 + num_pos + posi1;
			}

			inline int64_t index_pos_pos(int posi1, int posi2) const {
				if (posi2 > posi1) {
					int swap = posi1;
					posi1 = posi2;
					posi2 = swap;
				}
				return 1 + 2*num_pos + posi1*(posi1 - 1)/2 + posi2;
			}

			inline int64_t index(int posi1, int posi2) const {
				if (posi1 == posi2) {
					if (posi1 == StaticPos) {
						return index_static_static();
					} else {
						return index_pos(posi1);
					}
				} else if (posi1 == StaticPos) {
					return index_static_pos(posi2);
				} else if (posi2 == StaticPos) {
					return index_static_pos(posi1);
				} else {
					return index_pos_pos(posi1, posi2);
				}
			}

			inline const void * get_static_static_pair() const {
				return reinterpret_cast<const void *>(offset(pos_pairs()[index_static_static()]));
			}

			inline const void * get_static_pos_pairs(int posi1, int fragi1) const {
				auto frag_offsets = reinterpret_cast<const int64_t *>(offset(pos_pairs()[index_static_pos(posi1)]));
				return offset(frag_offsets[fragi1]);
			}

			inline const void * get_pos_pairs(int posi1, int fragi1) const {
				auto frag_offsets = reinterpret_cast<const int64_t *>(offset(pos_pairs()[index_pos(posi1)]));
				return offset(frag_offsets[fragi1]);
			}

			inline const void * get_pos_pos_pairs(int posi1, int fragi1, int posi2, int fragi2) const {
				if (posi2 > posi1) {
					int swap = posi1;
					posi1 = posi2;
					posi2 = swap;
				}
				auto frag_offsets = reinterpret_cast<const int64_t *>(offset(pos_pairs()[index_pos_pos(posi1, posi2)]));
				return offset(frag_offsets[fragi1*get_pos(posi2).num_frags + fragi2]);
			}

			inline int64_t get_molecule_motion_id(int i) const {

				// just in case ...
				assert (i >= 0);
				assert (i < num_molecule_motions);

				auto offsets = reinterpret_cast<const int64_t *>(offset(molecule_motions_offset));
				return *reinterpret_cast<const int64_t *>(offset(offsets[i]));
			}

			inline const void * get_molecule_motion(int i) const {

				// just in case ...
				assert (i >= 0);
				assert (i < num_molecule_motions);

				auto offsets = reinterpret_cast<const int64_t *>(offset(molecule_motions_offset));
				auto p = reinterpret_cast<const int64_t *>(offset(offsets[i]));
				return reinterpret_cast<const void *>(p + 1);
			}

			inline int64_t get_conf_motion_id(const Conf<T> & conf, int i) const {

				// just in case ...
				assert (i >= 0);
				assert (i < conf.num_motions);

				auto offsets = reinterpret_cast<const int64_t *>(offset(conf.motions_offset));
				return *reinterpret_cast<const int64_t *>(offset(offsets[i]));
			}

			inline const void * get_conf_motion(const Conf<T> & conf, int i) const {

				// just in case ...
				assert (i >= 0);
				assert (i < conf.num_motions);

				auto offsets = reinterpret_cast<const int64_t *>(offset(conf.motions_offset));
				auto p = reinterpret_cast<const int64_t *>(offset(offsets[i]));
				return reinterpret_cast<const void *>(p + 1);
			}

			// memory layout of conf_space
			int32_t num_pos;
			int32_t max_num_conf_atoms;
			int32_t max_num_dofs;
			int32_t num_molecule_motions;
			int64_t positions_offset;
			int64_t static_atom_coords_offset;
			int64_t static_atom_molis_offset;
			int64_t params_offset;
			int64_t pos_pairs_offset;
			int64_t molecule_motions_offset;
			T static_energy;
			// 4 byte pad, if T = float32_t
			//
			// (THIS IS WHERE POSITIONS OFFSET POINTS TO)
			// Int 64 array of position offsets of size confSpace.positions.length
			// For each position:
			//	The position
			//	Int 64 array of conf offsets of size pos.confs.length
			//	For each conf:
			//		The conf 
			//		Int 64 array of atom coords offsets of size conf.coords.size
			//		A struct containing the number of atom coords and a null ptr
			//		For each coordinate:
			//			(X, Y, Z, pad) floats
			//		A struct containing the number of atom molecule indices and a null ptr
			//		For each atom:
			//			an int of atom mol info index
			//		(pad to alignment of 8 bytes)
			//		Int 64 array of conformation motion offsets of size conf.motion.length
			//		For each motion:
			//			an int64 dihedral ID
			//			a dihedral struct
			//			an array of rotated indices of length desc.rotated.length
			//			(pad alignment to 8 bytes)
			//	(THIS IS WHERE STATIC ATOM COORDS OFFSET POINTS TO)
			//	A struct containing the number of static coords and a null ptr
			//	For each static coordinate:
			//		(X, Y, Z, pad) floats
			//	(THIS IS WHERE STATIC ATOM MOLIS OFFSET POINTS TO)
			//	A struct containing the number of static coords atom offsets and a null ptr
			//	For each static coordinate:
			//		An int of the index of the static atom molecule
			//	(pad to alignment of 8 bytes)
			//	(THIS IS WHERE THE PARAMS OFFSET POINTS TO)
			//	A forcefield params struct with the distance_dependent_dielectric (1 byte) and 7 bytes padding
			//	(THIS IS WHERE THE POSITION PAIRS OFFSET POINTS TO)
			//	Int 64 array of position pair offsets of size numPosPairs
			//	The first value in the array is the address after the end of the array
			//	An atom pairs struct
			//	For each amber atom pair:
			//		an amberStruct
			//	For each eff1 atom pair:
			//		an eef1Struct
			//	For each design position:
			//		Int 64 array of fragment pair offsets of size confSpace.numFrag()
			//		For each fragment:
			//			An atom pairs struct
			//			For each amber atom pair:
			//				an amberStruct
			//			For each eff1 atom pair:
			//				an eef1Struct
			//	For each design position:
			//		Int 64 array of fragment pair offsets of size confSpace.numFrag()
			//		For each fragment:
			//			An atom pairs struct
			//			For each amber atom pair:
			//				an amberStruct
			//			For each eff1 atom pair:
			//				an eef1Struct
			//	For each design position pair:
			//		Int 64 array of fragment pair offsets of size confSpace.numFrag(pos1) * confSpace.numFrag(pos2)
			//		For each fragment pair:
			//			An atom pairs struct
			//			For each amber atom pair:
			//				an amberStruct
			//			For each eff1 atom pair:
			//				an eef1Struct
			//	(THIS IS WHERE THE MOLECULE MOTIONS OFFSET POINT TO)
			//	Int 64 array of molecule motions offsets of size numMolMotions
			//	For each molecule:
			//		For each motion:
			//			If translation-rotation:
			//				transRotStruct
			//			If Dihedral motion:
			//				a dihedral struct
			//				an array of rotated indices of length desc.rotated.length
			//				(pad alignment to 8 bytes)

			inline const uint8_t * offset(int64_t offset) const {
				return reinterpret_cast<const uint8_t *>(this) + offset;
			}

			inline const int64_t * pos_pairs() const {
				return reinterpret_cast<const int64_t *>(offset(pos_pairs_offset));
			}
	};
	ASSERT_JAVA_COMPATIBLE_REALS(ConfSpace, 72, 72);

	template<typename T>
	std::ostream & operator << (std::ostream & out, const ConfSpace<T> & conf_space) {

		// show the conf space
		out << "ConfSpace:" << std::endl;
		for (int posi=0; posi<conf_space.get_num_pos(); posi++) {
			const Pos & pos = conf_space.get_pos(posi);
			out << "\tPos[" << posi << "]:" << std::endl;
			out << "\t\tconfs: " << pos.num_confs << std::endl;
			out << "\t\tmax num atoms: " << pos.max_num_atoms << std::endl;

			for (int confi=0; confi<pos.num_confs; confi++) {
				const Conf<T> & conf = conf_space.get_conf(pos, confi);
				const Array<Real3<T>> & atoms = conf_space.get_conf_atoms(conf);
				out << "\t\tConf[" << confi << "]:" << std::endl;
				out << "\t\t\tatoms: " << atoms.size() << std::endl;
				out << fmt::set_indents(4) << atoms << fmt::set_indents(0);
			}
		}

		// show the static atoms
		const Array<Real3<T>> & static_atoms = conf_space.get_static_atoms();
		out << "\tstatic atoms: " << static_atoms.size() << std::endl;
		out << fmt::set_indents(2) << static_atoms << fmt::set_indents(0);

		return out;
	}
}

#endif //CONFECALC_CONFSPACE_H
