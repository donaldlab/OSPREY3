
#include "global.h"
#include "math.h"
#include "confSpace.h"


typedef struct PosInter {
	uint32_t posi1;
	uint32_t posi2;
	real_t weight;
	real_t offset;
} ALIGN_8 PosInter;

void print_PosInter(const PosInter * inter) {
	printf("[posi1=%d  posi2=%d  weight=%f  offset=%f]",
		inter->posi1, inter->posi2, inter->weight, inter->offset
	);
}

void print_PosInters(const PosInter * inters, uint32_t inters_size) {
	printf("PosInter[%d]:\n", inters_size);
	for (int i=0; i<inters_size; i++) {
		printf("\t");
		print_PosInter(&inters[i]);
		printf("\n");
	}
}

void print_conf(const uint32_t * conf, uint32_t size) {
	printf("[");
	for (int i=0; i<size; i++) {
		if (i > 0) {
			printf(", ");
		}
		printf("%d", conf[i]);
	}
	printf("]");
}

void print_info() {
	printf("ConfEcalc " VAL(REAL) "  v" VAL(ConfEcalc_VERSION_MAJOR) "." VAL(ConfEcalc_VERSION_MINOR) "\n");
	fflush(stdout);
}

void assign(const ConfSpace * conf_space, const int32_t conf[], real3 out[]) {
	AssignedCoords * assigned_coords = conf_space_assign(conf_space, conf);
	coords_copy(out, assigned_coords->atoms->coords, conf_space->max_num_conf_atoms);
	assigned_coords_free(assigned_coords);
}

real_t calc(const ConfSpace * conf_space, const int32_t conf[], const PosInter inters[], uint32_t inters_size) {

	// TODO: refactor these into helper functions and move them into a lib somewhere

	// TEMP
	printf("conf space:\n");
	printf("\tmax_num_conf_atoms: %d\n", conf_space->max_num_conf_atoms);
	printf("conf space positions: %d\n", conf_space->num_pos);
	for (int posi=0; posi<conf_space->num_pos; posi++) {
		const Pos * pos = conf_space_pos(conf_space, posi);
		printf("\tpos %d, confs=%d, max_num_atoms=%d\n", posi, pos->num_confs, pos->max_num_atoms);

		for (int confi=0; confi<pos->num_confs; confi++) {
			const Conf * conf = conf_space_conf(conf_space, pos, confi);
			const Coords * atoms = conf_space_conf_atoms(conf_space, conf);
			printf("\t\tconf %d, atoms = %d\n", confi, atoms->num_atoms);

			/* TEMP
			for (int atomi=0; atomi<atoms->num_atoms; atomi++) {
				printf("\t\t\t");
				print_real3(&atoms->coords[atomi]);
				printf("\n");
			}
			*/
		}
	}

	// TEMP
	const Coords * static_atoms = conf_space_static_atoms(conf_space);
	printf("static atoms: %d\n", static_atoms->num_atoms);
	/* TEMP
	for (int atomi=0; atomi<static_atoms->num_atoms; atomi++) {
		printf("\t");
		print_real3(&static_atoms->coords[atomi]);
		printf("\n");
	}
	*/

	// TEMP
	printf("conf: ");
	print_conf(conf, conf_space->num_pos);
	printf("\n");

	// TEMP
	printf("inters: ");
	print_PosInters(inters, inters_size);

	// TEMP
	AssignedCoords * assigned_coords = conf_space_assign(conf_space, conf);
	printf("AssignedCoords: %d atoms", assigned_coords->atoms->num_atoms);
	assigned_coords_free(assigned_coords);

	return RL(5.0);
};

