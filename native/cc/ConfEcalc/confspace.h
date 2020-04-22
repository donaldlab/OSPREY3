
#ifndef CONFECALC_CONFSPACE_H
#define CONFECALC_CONFSPACE_H


namespace osprey {

	struct Conf {
		int64_t atoms_offset;
	};
	ASSERT_JAVA_COMPATIBLE(Conf, 8);

	struct Pos {
		int32_t num_confs;
		int32_t max_num_atoms;
		// java puts the confs offsets after the struct
	};
	ASSERT_JAVA_COMPATIBLE(Pos, 8);

	template<typename T>
	class ConfSpace {

		public:
			ConfSpace() = delete; // created only on the Java side

			inline int32_t get_num_pos() const {
				return num_pos;
			}

			inline int32_t get_max_num_conf_atoms() const {
				return max_num_conf_atoms;
			}

			inline const Pos & get_pos(int posi) const {
				auto offsets = reinterpret_cast<const int64_t *>(offset(positions_offset));
				return *reinterpret_cast<const Pos *>(offset(offsets[posi]));
			}

			inline const Conf & get_conf(const Pos & pos, int confi) const {
				auto conf_offsets = reinterpret_cast<const int64_t *>(&pos + 1);
				return *reinterpret_cast<const Conf *>(offset(conf_offsets[confi]));
			}

			inline const Atoms<T> & get_static_atoms() const {
				return *reinterpret_cast<const Atoms<T> *>(offset(static_atoms_offset));
			}

			inline const Atoms<T> & get_conf_atoms(const Conf & conf) const {
				return *reinterpret_cast<const Atoms<T> *>(offset(conf.atoms_offset));
			}

		private:
			int32_t num_pos;
			int32_t max_num_conf_atoms;
			int64_t positions_offset;
			int64_t static_atoms_offset;

			inline const uint8_t * offset(int64_t offset) const {
				return reinterpret_cast<const uint8_t *>(this) + offset;
			}
	};
	ASSERT_JAVA_COMPATIBLE_REALS(ConfSpace, 24, 24);

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
				const Conf & conf = conf_space.get_conf(pos, confi);
				const Atoms<T> & atoms = conf_space.get_conf_atoms(conf);
				out << "\t\tConf[" << confi << "]:" << std::endl;
				out << "\t\t\tatoms: " << atoms.size() << std::endl;
				out << fmt::set_indents(4) << atoms << fmt::set_indents(0);
			}
		}

		// show the static atoms
		const Atoms<T> & static_atoms = conf_space.get_static_atoms();
		out << "\tstatic atoms: " << static_atoms.size() << std::endl;
		out << fmt::set_indents(2) << static_atoms << fmt::set_indents(0);

		return out;
	}
}

#endif //CONFECALC_CONFSPACE_H
