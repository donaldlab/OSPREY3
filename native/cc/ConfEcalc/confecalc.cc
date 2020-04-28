
#include "global.h"
#include "formats.h"
#include "real3.h"
#include "rotation.h"
#include "atoms.h"
#include "confspace.h"
#include "assignment.h"
#include "energy.h"
#include "energy_ambereef1.h"
#include "motions.h"
#include "minimization.h"


#define API extern "C" [[maybe_unused]]


API int version_major() {
	return ConfEcalc_VERSION_MAJOR;
}

API int version_minor() {
	return ConfEcalc_VERSION_MINOR;
}

template<typename T>
static void assign(const osprey::ConfSpace<T> & conf_space, const int32_t conf[], osprey::Real3<T> out[]) {
	osprey::Assignment<T> assignment(conf_space, conf);
	std::copy(&assignment.atoms[0], &assignment.atoms[assignment.atoms.size() - 1], out);
}

API void assign_f32(const osprey::ConfSpace<float32_t> & conf_space, const int32_t conf[], osprey::Real3<float32_t> out[]) {
	assign(conf_space, conf, out);
}
API void assign_f64(const osprey::ConfSpace<float64_t> & conf_space, const int32_t conf[], osprey::Real3<float64_t> out[]) {
	assign(conf_space, conf, out);
}


template<typename T>
static T calc(const osprey::ConfSpace<T> & conf_space, const int32_t conf[], const osprey::PosInter<T> inters[],
              int64_t inters_size, osprey::EnergyFunction<T> efunc, osprey::Real3<T> out[]) {

	osprey::Assignment<T> assignment(conf_space, conf);
	T energy = efunc(assignment, inters, inters_size);
	if (out != nullptr) {
		std::copy(&assignment.atoms[0], &assignment.atoms[assignment.atoms.size() - 1], out);
	}
	return energy;
}

API float32_t calc_amber_eef1_f32(const osprey::ConfSpace<float32_t> & conf_space, const int32_t conf[],
                                  const osprey::PosInter<float32_t> inters[], int64_t inters_size,
                                  osprey::Real3<float32_t> out_coords[]) {
	return calc<float32_t>(conf_space, conf, inters, inters_size, osprey::ambereef1::calc_energy, out_coords);
}
API float64_t calc_amber_eef1_f64(const osprey::ConfSpace<float64_t> & conf_space, const int32_t conf[],
                                  const osprey::PosInter<float64_t> inters[], int64_t inters_size,
                                  osprey::Real3<float64_t> out_coords[]) {
	return calc<float64_t>(conf_space, conf, inters, inters_size, osprey::ambereef1::calc_energy, out_coords);
}

template<typename T>
static T minimize(const osprey::ConfSpace<T> & conf_space, const int32_t conf[],
                  const osprey::PosInter<T> inters[], int64_t inters_size,
                  osprey::EnergyFunction<T> efunc,
                  osprey::Real3<T> out_coords[], uint8_t * out_dofs) {

	osprey::Assignment<T> assignment(conf_space, conf);
	osprey::Minimization<T> minimization(assignment, inters, inters_size);
	minimization.minimize(efunc);
	if (out_coords != nullptr) {
		std::copy(&assignment.atoms[0], &assignment.atoms[assignment.atoms.size() - 1], out_coords);
	}
	if (out_dofs != nullptr) {
		int64_t & out_dofs_size = *reinterpret_cast<int64_t *>(out_dofs);
		T * out_dofs_values = reinterpret_cast<T *>(out_dofs + sizeof(int64_t));

		for (int i=0; i<minimization.get_num_dofs(); i++) {
			out_dofs_values[i] = minimization.get_dof(i);
		}
		out_dofs_size = minimization.get_num_dofs();
	}
	return minimization.energy;
}

API float32_t minimize_amber_eef1_f32(const osprey::ConfSpace<float32_t> & conf_space, const int32_t conf[],
                                      const osprey::PosInter<float32_t> inters[], int64_t inters_size,
                                      osprey::Real3<float32_t> out_coords[], uint8_t * out_dofs) {
	return minimize<float32_t>(conf_space, conf, inters, inters_size, osprey::ambereef1::calc_energy, out_coords, out_dofs);
}
API float64_t minimize_amber_eef1_f64(const osprey::ConfSpace<float64_t> & conf_space, const int32_t conf[],
                                      const osprey::PosInter<float64_t> inters[], int64_t inters_size,
                                      osprey::Real3<float64_t> out_coords[], uint8_t * out_dofs) {
	return minimize<float64_t>(conf_space, conf, inters, inters_size, osprey::ambereef1::calc_energy, out_coords, out_dofs);
}
