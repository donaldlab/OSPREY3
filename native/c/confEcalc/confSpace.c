
#include "global.h"
#include "math.h"
#include "confSpace.h"


// the FN macros really confused the IDE though, so turn off some warnings
#pragma ide diagnostic ignored "InfiniteRecursion"


const int32_t UNASSIGNED = -1;
const int32_t STATIC_POS = -1;

Coords * coords_malloc(uint32_t num_atoms) {
	Coords * coords = (Coords *)malloc(sizeof(Coords) + sizeof(real3)*num_atoms);
	coords->num_atoms = num_atoms;
	return coords;
}

void coords_copy(real3 * dst, const real3 * src, uint64_t num_atoms) {
    memcpy(dst, src, num_atoms*sizeof(real3));
}

void * conf_space_offset(const ConfSpace * conf_space, uint64_t offset) {
	return ((void *)conf_space) + offset;
}

uint64_t * conf_space_pos_offsets(const ConfSpace * conf_space) {
	return (uint64_t *)conf_space_offset(conf_space, conf_space->positions_offset);
}

const Pos * conf_space_pos(const ConfSpace * conf_space, int posi) {
	uint64_t * pos_offsets = conf_space_pos_offsets(conf_space);
	return (Pos *)conf_space_offset(conf_space, pos_offsets[posi]);
}

const Coords * conf_space_static_atoms(const ConfSpace * conf_space) {
	return (Coords *)conf_space_offset(conf_space, conf_space->static_atoms_offset);
}

const Conf * conf_space_conf(const ConfSpace * conf_space, const Pos * pos, int confi) {
	return (Conf *)conf_space_offset(conf_space, pos->conf_offsets[confi]);
}

const Coords * conf_space_conf_atoms(const ConfSpace * conf_space, const Conf * conf) {
	return (Coords *)conf_space_offset(conf_space, conf->coords_offset);
}


void assigned_coords_free(AssignedCoords * coords) {
	free(coords->atoms);
	free(coords);
}

AssignedCoords * conf_space_assign(const ConfSpace * conf_space, const int32_t conf[]) {
	
	// allocate the coords
	AssignedCoords * coords = (AssignedCoords *)malloc(sizeof(AssignedCoords) + sizeof(real3 *)*conf_space->num_pos);
	coords->conf_space = conf_space;
	coords->conf = conf;
	coords->atoms = coords_malloc(conf_space->max_num_conf_atoms);

	// copy the static atoms
	const Coords * static_atoms = conf_space_static_atoms(conf_space);
	coords_copy(coords->atoms->coords, static_atoms->coords, static_atoms->num_atoms);
	uint32_t atoms_offset = static_atoms->num_atoms;

	for (int posi=0; posi<conf_space->num_pos; posi++) {
		const Pos * pos = conf_space_pos(conf_space, posi);
		int32_t confi = conf[posi];

		// set the position atom pointers
		coords->pos_atoms[posi] = &coords->atoms->coords[atoms_offset];
        atoms_offset += pos->max_num_atoms;

		// get the assigned conformation at this position, if any
		if (confi != UNASSIGNED) {

            // copy the atom coords
            const Conf * posConf = conf_space_conf(conf_space, pos, confi);
            const Coords * atoms = conf_space_conf_atoms(conf_space, posConf);
            coords_copy(coords->pos_atoms[posi], atoms->coords, atoms->num_atoms);

            // zero out any remaining space for this pos
            memset(&coords->pos_atoms[posi][atoms->num_atoms], 0, (pos->max_num_atoms - atoms->num_atoms)*sizeof(real3));

        } else {

		    // otherwise, just zero out the coords for this pos
            memset(coords->pos_atoms[posi], 0, pos->max_num_atoms*sizeof(real3));
		}
	}

	return coords;
}

real3 * assigned_coords_atoms(AssignedCoords * assigned_coords, uint32_t posi) {
	if (posi == STATIC_POS) {
		return assigned_coords->atoms->coords;
	} else {
		return assigned_coords->pos_atoms[posi];
	}
}
