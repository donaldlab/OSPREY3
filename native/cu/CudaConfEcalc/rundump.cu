
#include <fstream>
#include <endian.h>

#include "global.h"
#include "cuda.h"
#include "array.h"
#include "formats.h"
#include "rotation.h"
#include "atoms.h"
#include "confspace.h"
#include "assignment.h"
#include "energy.h"
#include "energy_ambereef1.h"
#include "motions.h"
#include "minimization.h"
#include "confecalc.h"


int32_t read_be32(std::ifstream & in) {
	int32_t i;
	in.read((char *)&i, sizeof(int32_t));
	i = be32toh(i);
	return i;
}

enum Precision {
	Float32 = 0,
	Float64 = 1
};

class Buf {
	public:

		explicit Buf(std::ifstream & in) {
			size = read_be32(in);
			buf = new char[size];
			in.read(buf, size);
		}

		~Buf() {
			delete[] buf;
		}

		int32_t get_size() const {
			return size;
		}

		const char * ptr() const {
			return buf;
		}

	private:
		int32_t size;
		char * buf;
};


namespace osprey {

	template <typename T, EnergyFunction<T> efunc>
	void run(const Buf & conf_space, const Buf & conf_space_sizes, const Buf & jobs) {

		const ConfSpace<T> & h_conf_space = *reinterpret_cast<const ConfSpace<T> *>(conf_space.ptr());
		auto h_conf_space_sizes = reinterpret_cast<const ConfSpaceSizes *>(conf_space_sizes.ptr());
		auto h_jobs = reinterpret_cast<const MinimizationJobs<T> *>(jobs.ptr());

		int64_t n = h_jobs->get_size();
		Array<T> * h_energies = Array<T>::make(n);

		// init the device
		int device = 0;
		MinimizationBuffers<T> buffers(h_conf_space_sizes, n);
		Stream stream(buffers.get_bytes_host(), buffers.get_bytes_device());
		auto d_conf_space = reinterpret_cast<ConfSpace<T> *>(alloc_conf_space<T>(device, h_conf_space));

		// minimize it!
		std::cout << "Minimizing " << n << " jobs ..." << std::endl;
		minimize_batch<T,efunc>(device, &stream, d_conf_space, h_conf_space_sizes, h_jobs, h_energies);

		// show the energies
		for (int i=0; i<n; i++) {
			std::cout << "energy[" << i << "] = " << (*h_energies)[i] << std::endl;
		}

		// cleanup
		free_conf_space(device, d_conf_space);
		delete h_energies;
	}
}

int main(int argc, char * argv[]) {

	// handle args
	if (argc != 2) {
		std::cout << "Usage:\n\trundump dumpfile\n" << std::endl;
		return -1;
	}

	// get the dump file
	char * path = argv[1];
	std::cout << "Reading dump file: " << path << std::endl;
	std::ifstream file(path, std::ios::binary);
	if (!file.is_open()) {
		std::cout << "can't read dump file!" << std::endl;
		return -2;
	}

	// read the dump file
	int32_t precision = read_be32(file);
	const Buf conf_space = Buf(file);
	const Buf conf_space_sizes = Buf(file);
	const Buf jobs = Buf(file);

	// run the minimizations
	switch (precision) {
		case Float32:
			osprey::run<float32_t,osprey::ambereef1::calc_energy>(conf_space, conf_space_sizes, jobs);
		break;
		case Float64:
			osprey::run<float64_t,osprey::ambereef1::calc_energy>(conf_space, conf_space_sizes, jobs);
		break;
		default:
			std::cout << "file is not a dump file!" << std::endl;
			return -3;
	}

	return 0;
}
