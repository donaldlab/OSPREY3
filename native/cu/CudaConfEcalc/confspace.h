
#ifndef CONFECALC_CONFSPACE_H
#define CONFECALC_CONFSPACE_H


namespace osprey {

	const int32_t StaticPos = -1;

	template<typename T>
	struct alignas(16) Conf {
		int64_t atoms_offset;
		int32_t frag_index;
		// 4 bytes pad if T = float64_t
		T internal_energy;
		int64_t num_motions;
		int64_t motions_offset;
		// 8 bytes pad if T = float64_t
	};
	ASSERT_JAVA_COMPATIBLE_REALS(Conf, 32, 48);

	struct alignas(16) Pos {
		int32_t num_confs;
		int32_t max_num_atoms;
		int32_t num_frags;
		// 4 bytes pad
	};
	ASSERT_JAVA_COMPATIBLE(Pos, 16);

	template<typename T>
	class alignas(16) ConfSpace {
		public:
			ConfSpace() = delete; // created only on the Java side

			__host__ __device__
			inline const Pos & get_pos(int posi) const {
				auto offsets = reinterpret_cast<const int64_t *>(offset(positions_offset));
				return *reinterpret_cast<const Pos *>(offset(offsets[posi]));
			}

			__host__ __device__
			inline const Conf<T> & get_conf(const Pos & pos, int confi) const {
				// the java side puts the conf offsets right after the Pos struct
				auto conf_offsets = reinterpret_cast<const int64_t *>(&pos + 1);
				return *reinterpret_cast<const Conf<T> *>(offset(conf_offsets[confi]));
			}

			__host__ __device__
			inline const Array<Real3<T>> & get_static_atoms() const {
				return *reinterpret_cast<const Array<Real3<T>> *>(offset(static_atoms_offset));
			}

			__host__ __device__
			inline const Array<Real3<T>> & get_conf_atoms(const Conf<T> & conf) const {
				return *reinterpret_cast<const Array<Real3<T>> *>(offset(conf.atoms_offset));
			}

			__host__ __device__
			inline const void * get_params() const {
				return reinterpret_cast<const void *>(offset(params_offset));
			}

			__host__ __device__
			inline int64_t index_static_static() const {
				return 0;
			}

			__host__ __device__
			inline int64_t index_static_pos(int posi1) const {
				return 1 + posi1;
			}

			__host__ __device__
			inline int64_t index_pos(int posi1) const {
				return 1 + num_pos + posi1;
			}

			__host__ __device__
			inline int64_t index_pos_pos(int posi1, int posi2) const {
				if (posi2 > posi1) {
					int swap = posi1;
					posi1 = posi2;
					posi2 = swap;
				}
				return 1 + 2*num_pos + posi1*(posi1 - 1)/2 + posi2;
			}

			__host__ __device__
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

			__host__ __device__
			inline const void * get_static_static_pair() const {
				return reinterpret_cast<const void *>(offset(pos_pairs()[index_static_static()]));
			}

			__host__ __device__
			inline const void * get_static_pos_pairs(int posi1, int fragi1) const {
				auto frag_offsets = reinterpret_cast<const int64_t *>(offset(pos_pairs()[index_static_pos(posi1)]));
				return offset(frag_offsets[fragi1]);
			}

			__host__ __device__
			inline const void * get_pos_pairs(int posi1, int fragi1) const {
				auto frag_offsets = reinterpret_cast<const int64_t *>(offset(pos_pairs()[index_pos(posi1)]));
				return offset(frag_offsets[fragi1]);
			}

			__host__ __device__
			inline const void * get_pos_pos_pairs(int posi1, int fragi1, int posi2, int fragi2) const {
				if (posi2 > posi1) {
					int swap = posi1;
					posi1 = posi2;
					posi2 = swap;
				}
				auto frag_offsets = reinterpret_cast<const int64_t *>(offset(pos_pairs()[index_pos_pos(posi1, posi2)]));
				return offset(frag_offsets[fragi1*get_pos(posi2).num_frags + fragi2]);
			}

			__host__ __device__
			inline int64_t get_molecule_motion_id(int i) const {

				// just in case ...
				assert (i >= 0);
				assert (i < num_molecule_motions);

				auto offsets = reinterpret_cast<const int64_t *>(offset(molecule_motions_offset));
				return *reinterpret_cast<const int64_t *>(offset(offsets[i]));
			}

			__host__ __device__
			inline const void * get_molecule_motion(int i) const {

				// just in case ...
				assert (i >= 0);
				assert (i < num_molecule_motions);

				auto offsets = reinterpret_cast<const int64_t *>(offset(molecule_motions_offset));
				auto p = reinterpret_cast<const int64_t *>(offset(offsets[i]));
				return reinterpret_cast<const void *>(p + 1);
			}

			__host__ __device__
			inline int64_t get_conf_motion_id(const Conf<T> & conf, int i) const {

				// just in case ...
				assert (i >= 0);
				assert (i < conf.num_motions);

				auto offsets = reinterpret_cast<const int64_t *>(offset(conf.motions_offset));
				return *reinterpret_cast<const int64_t *>(offset(offsets[i]));
			}

			__host__ __device__
			inline const void * get_conf_motion(const Conf<T> & conf, int i) const {

				// just in case ...
				assert (i >= 0);
				assert (i < conf.num_motions);

				auto offsets = reinterpret_cast<const int64_t *>(offset(conf.motions_offset));
				return reinterpret_cast<const void *>(offset(offsets[i]));
			}

			int32_t num_pos;
			int32_t max_num_conf_atoms;
			int32_t max_num_dofs;
			int32_t num_molecule_motions;
			int64_t size;
			int64_t positions_offset;
			int64_t static_atoms_offset;
			int64_t params_offset;
			int64_t pos_pairs_offset;
			int64_t molecule_motions_offset;
			T static_energy;
			// 4 byte pad, if T = float32_t

			__host__ __device__
			inline const uint8_t * offset(int64_t offset) const {
				return reinterpret_cast<const uint8_t *>(this) + offset;
			}

			__host__ __device__
			inline const int64_t * pos_pairs() const {
				return reinterpret_cast<const int64_t *>(offset(pos_pairs_offset));
			}
	};
	ASSERT_JAVA_COMPATIBLE_REALS(ConfSpace, 80, 80);

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

	struct ConfSpaceSizes {
		int64_t num_pos;
		int64_t max_num_inters;
		int64_t num_atoms;
		int64_t max_num_dofs;
	};
	ASSERT_JAVA_COMPATIBLE(ConfSpaceSizes, 32);
}

#endif //CONFECALC_CONFSPACE_H
