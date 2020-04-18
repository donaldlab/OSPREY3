
#ifndef _CONFSPACE_
#define _CONFSPACE_


typedef struct Coords {
	uint32_t num_atoms;
	// 4 bytes pad
	real3 coords[];
} ALIGN_8 Coords;

Coords * coords_malloc(uint32_t num_atoms);
void coords_copy(real3 * dst, const real3 * src, uint64_t num_atoms);


typedef struct Conf {
	uint64_t coords_offset;
} ALIGN_8 Conf;

typedef struct Pos {
	uint32_t num_confs;
	uint32_t max_num_atoms;
	uint64_t conf_offsets[];
} ALIGN_8 Pos;

typedef struct ConfSpace {
	uint32_t num_pos;
	uint32_t max_num_conf_atoms;
	uint64_t positions_offset;
	uint64_t static_atoms_offset;
} ALIGN_8 ConfSpace;


const Pos * conf_space_pos(const ConfSpace * conf_space, int posi);
const Coords * conf_space_static_atoms(const ConfSpace * conf_space);
const Conf * conf_space_conf(const ConfSpace * conf_space, const Pos * pos, int confi);
const Coords * conf_space_conf_atoms(const ConfSpace * conf_space, const Conf * conf);


extern const int32_t UNASSIGNED;
extern const int32_t STATIC_POS;


typedef struct AssignedCoords {
	const ConfSpace * conf_space;
	const int32_t * conf;
	Coords * atoms;
	real3 * pos_atoms[]; // indexed by posi
} AssignedCoords;


void assigned_coords_free(AssignedCoords * coords);
AssignedCoords * conf_space_assign(const ConfSpace * conf_space, const int32_t conf[]);
real3 * assigned_coords_atoms(AssignedCoords * assigned_coords, uint32_t posi);


#endif // _CONFSPACE_

