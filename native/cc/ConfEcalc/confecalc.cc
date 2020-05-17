
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
	return ConfEcalc_VERSION_MAJOR;
}

API int version_minor() {
	return ConfEcalc_VERSION_MINOR;
}


namespace osprey {

	template<typename T>
	static void assign(const ConfSpace<T> & conf_space, const int32_t conf[], Array<Real3<T>> & out_coords) {
		Assignment<T> assignment(conf_space, conf);
		out_coords.copy_from(assignment.atoms);
	}
}

API void assign_f32(const osprey::ConfSpace<float32_t> & conf_space, const int32_t conf[], osprey::Array<osprey::Real3<float32_t>> & out_coords) {
	osprey::assign(conf_space, conf, out_coords);
}
API void assign_f64(const osprey::ConfSpace<float64_t> & conf_space, const int32_t conf[], osprey::Array<osprey::Real3<float64_t>> & out_coords) {
	osprey::assign(conf_space, conf, out_coords);
}


namespace osprey {

	template<typename T>
	static T calc(const ConfSpace<T> & conf_space,
	              const int32_t conf[],
	              const Array<PosInter<T>> & inters,
	              EnergyFunction<T> efunc,
	              Array<Real3<T>> * out_coords) {

		Assignment<T> assignment(conf_space, conf);
		T energy = efunc(assignment, inters);
		if (out_coords != nullptr) {
			out_coords->copy_from(assignment.atoms);
		}
		return energy;
	}
}

API float32_t calc_amber_eef1_f32(const osprey::ConfSpace<float32_t> & conf_space,
                                  const int32_t conf[],
                                  const osprey::Array<osprey::PosInter<float32_t>> & inters,
                                  osprey::Array<osprey::Real3<float32_t>> * out_coords) {
	return osprey::calc<float32_t>(conf_space, conf, inters, osprey::ambereef1::calc_energy, out_coords);
}
API float64_t calc_amber_eef1_f64(const osprey::ConfSpace<float64_t> & conf_space,
                                  const int32_t conf[],
                                  const osprey::Array<osprey::PosInter<float64_t>> & inters,
                                  osprey::Array<osprey::Real3<float64_t>> * out_coords) {
	return osprey::calc<float64_t>(conf_space, conf, inters, osprey::ambereef1::calc_energy, out_coords);
}


namespace osprey {

	template<typename T>
	static T minimize(const ConfSpace<T> & conf_space,
	                  const int32_t conf[],
	                  const Array<PosInter<T>> & inters,
	                  EnergyFunction<T> efunc,
	                  Array<Real3<T>> * out_coords,
	                  Array<T> * out_dofs) {

		// make the coords and the degrees of freedom
		Assignment<T> assignment(conf_space, conf);
		Dofs<T> dofs(assignment, inters, efunc);

		// init the dofs to the center of the voxel
		DofValues<T> vals(dofs.get_size());
		for (int d=0; d<dofs.get_size(); d++) {
			vals.x[d] = dofs[d].center();
		}

		// use CCD to do the minimization
		minimize_ccd(dofs, vals);

		// set out the results
		if (out_coords != nullptr) {
			out_coords->copy_from(assignment.atoms);
		}
		if (out_dofs != nullptr) {
			out_dofs->copy_from(vals.x);
			out_dofs->truncate(dofs.get_size());
		}
		return vals.f;
	}
}

API float32_t minimize_amber_eef1_f32(const osprey::ConfSpace<float32_t> & conf_space,
                                      const int32_t conf[],
                                      const osprey::Array<osprey::PosInter<float32_t>> & inters,
                                      osprey::Array<osprey::Real3<float32_t>> * out_coords,
                                      osprey::Array<float32_t> * out_dofs) {
	return osprey::minimize<float32_t>(conf_space, conf, inters, osprey::ambereef1::calc_energy, out_coords, out_dofs);
}
API float64_t minimize_amber_eef1_f64(const osprey::ConfSpace<float64_t> & conf_space,
                                      const int32_t conf[],
                                      const osprey::Array<osprey::PosInter<float64_t>> & inters,
                                      osprey::Array<osprey::Real3<float64_t>> * out_coords,
                                      osprey::Array<float64_t> * out_dofs) {
	return osprey::minimize<float64_t>(conf_space, conf, inters, osprey::ambereef1::calc_energy, out_coords, out_dofs);
}
