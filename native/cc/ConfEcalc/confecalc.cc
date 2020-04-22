
#include "global.h"
#include "formats.h"
#include "real3.h"
#include "atoms.h"
#include "confspace.h"
#include "assignment.h"


#define API extern "C" [[maybe_unused]]


API int version_major() {
	return ConfEcalc_VERSION_MAJOR;
}

API int version_minor() {
	return ConfEcalc_VERSION_MINOR;
}

template<typename T>
void assign(const osprey::ConfSpace<T> & conf_space, const int32_t conf[], osprey::Real3<T> out[]) {
	osprey::Assignment<T> assignment(conf_space, conf);
	std::copy(&assignment.atoms[0], &assignment.atoms[assignment.atoms.size() - 1], out);
}

API void assign_f32(const osprey::ConfSpace<float32_t> & conf_space, const int32_t conf[], osprey::Real3<float32_t> out[]) {
	assign(conf_space, conf, out);
}
API void assign_f64(const osprey::ConfSpace<float64_t> & conf_space, const int32_t conf[], osprey::Real3<float64_t> out[]) {
	assign(conf_space, conf, out);
}
