
#include "global.h"
#include "array.h"
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
	return CudaConfEcalc_VERSION_MAJOR;
}

API int version_minor() {
	return CudaConfEcalc_VERSION_MINOR;
}

template<typename T>
static void assign(const osprey::ConfSpace<T> & conf_space, const int32_t conf[], osprey::Array<osprey::Real3<T>> & out_coords) {
	osprey::Assignment<T> assignment(conf_space, conf);
	out_coords.copy_from(assignment.atoms);
}

API void assign_f32(const osprey::ConfSpace<float32_t> & conf_space, const int32_t conf[], osprey::Array<osprey::Real3<float32_t>> & out_coords) {
	assign(conf_space, conf, out_coords);
}
API void assign_f64(const osprey::ConfSpace<float64_t> & conf_space, const int32_t conf[], osprey::Array<osprey::Real3<float64_t>> & out_coords) {
	assign(conf_space, conf, out_coords);
}


template<typename T>
static T calc(const osprey::ConfSpace<T> & conf_space, const int32_t conf[], const osprey::Array<osprey::PosInter<T>> & inters,
              osprey::EnergyFunction<T> efunc, osprey::Array<osprey::Real3<T>> * out_coords) {

	osprey::Assignment<T> assignment(conf_space, conf);
	T energy = efunc(assignment, inters);
	if (out_coords != nullptr) {
		out_coords->copy_from(assignment.atoms);
	}
	return energy;
}

API float32_t calc_amber_eef1_f32(const osprey::ConfSpace<float32_t> & conf_space, const int32_t conf[],
                                  const osprey::Array<osprey::PosInter<float32_t>> & inters,
                                  osprey::Array<osprey::Real3<float32_t>> * out_coords) {
	return calc<float32_t>(conf_space, conf, inters, osprey::ambereef1::calc_energy, out_coords);
}
API float64_t calc_amber_eef1_f64(const osprey::ConfSpace<float64_t> & conf_space, const int32_t conf[],
                                  const osprey::Array<osprey::PosInter<float64_t>> & inters,
                                  osprey::Array<osprey::Real3<float64_t>> * out_coords) {
	return calc<float64_t>(conf_space, conf, inters, osprey::ambereef1::calc_energy, out_coords);
}

template<typename T>
static T minimize(const osprey::ConfSpace<T> & conf_space, const int32_t conf[],
                  const osprey::Array<osprey::PosInter<T>> & inters,
                  osprey::EnergyFunction<T> efunc,
                  osprey::Array<osprey::Real3<T>> * out_coords, osprey::Array<T> * out_dofs) {

	// make the coords and the degrees of freedom
	osprey::Assignment<T> assignment(conf_space, conf);
	osprey::Dofs<T> dofs(assignment, inters, efunc);

	// init the dofs to the center of the voxel
	osprey::DofValues<T> vals(dofs.get_size());
	for (int d=0; d<dofs.get_size(); d++) {
		vals.x[d] = dofs[d].center();
	}

	// use CCD to do the minimization
	osprey::minimize_ccd(dofs, vals);

	// set out the results
	if (out_coords != nullptr) {
		out_coords->copy_from(assignment.atoms);
	}
	if (out_dofs != nullptr) {
		out_dofs->copy_from(vals.x, dofs.get_size());
		out_dofs->truncate(dofs.get_size());
	}
	return vals.f;
}

API float32_t minimize_amber_eef1_f32(const osprey::ConfSpace<float32_t> & conf_space, const int32_t conf[],
                                      const osprey::Array<osprey::PosInter<float32_t>> & inters,
                                      osprey::Array<osprey::Real3<float32_t>> * out_coords, osprey::Array<float32_t> * out_dofs) {
	return minimize<float32_t>(conf_space, conf, inters, osprey::ambereef1::calc_energy, out_coords, out_dofs);
}
API float64_t minimize_amber_eef1_f64(const osprey::ConfSpace<float64_t> & conf_space, const int32_t conf[],
                                      const osprey::Array<osprey::PosInter<float64_t>> & inters,
                                      osprey::Array<osprey::Real3<float64_t>> * out_coords, osprey::Array<float64_t> * out_dofs) {
	return minimize<float64_t>(conf_space, conf, inters, osprey::ambereef1::calc_energy, out_coords, out_dofs);
}
