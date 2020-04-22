
#ifndef CONFECALC_ASSIGNMENT_H
#define CONFECALC_ASSIGNMENT_H


namespace osprey {

	template<typename T>
	class Assignment {

		public:
			Assignment(const ConfSpace<T> & _conf_space, const int32_t _conf[])
				: conf_space(_conf_space), conf(_conf), atoms(_conf_space.get_max_num_conf_atoms()) {

				int32_t offset = 0;

				// copy the static atoms
				offset += atoms.copy_from(conf_space.get_static_atoms());

				// copy the atoms for each pos
				for (int posi=0; posi<conf_space.get_num_pos(); posi++) {
					const Pos & pos = conf_space.get_pos(posi);
					const Conf & pconf = conf_space.get_conf(pos, conf[posi]);

					// copy the atoms
					int64_t num_copied = atoms.copy_from(conf_space.get_conf_atoms(pconf), offset);
					offset += num_copied;

					// zero out the rest of the space for this pos
					int64_t atoms_remaining = pos.max_num_atoms - num_copied;
					std::memset(&atoms[offset], 0, sizeof(Real3<T>)*atoms_remaining);
					offset += atoms_remaining;
				}
			}

			const ConfSpace<T> & conf_space;
			const int32_t * conf;
			Atoms<T> atoms;
	};
}


#endif //CONFECALC_ASSIGNMENT_H
