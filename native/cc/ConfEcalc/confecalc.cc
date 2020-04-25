
#include "global.h"
#include "formats.h"
#include "real3.h"
#include "atoms.h"
#include "confspace.h"
#include "assignment.h"
#include "energy.h"
#include "energy_ambereef1.h"


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
static T calc_energy(const osprey::ConfSpace<T> & conf_space, const int32_t conf[], const osprey::PosInter<T> inters[], int64_t inters_size) {
	osprey::Assignment<T> assignment(conf_space, conf);
	return osprey::ambereef1::calc_energy(assignment, inters, inters_size);
}

API float32_t calc_energy_amber_eef1_f32(const osprey::ConfSpace<float32_t> & conf_space, const int32_t conf[], const osprey::PosInter<float32_t> inters[], int64_t inters_size) {
	return calc_energy<float32_t>(conf_space, conf, inters, inters_size);
}
API float64_t calc_energy_amber_eef1_f64(const osprey::ConfSpace<float64_t> & conf_space, const int32_t conf[], const osprey::PosInter<float64_t> inters[], int64_t inters_size) {
	return calc_energy<float64_t>(conf_space, conf, inters, inters_size);
}
