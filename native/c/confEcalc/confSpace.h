
#ifndef _CONFSPACE_
#define _CONFSPACE_


typedef struct Coords {
	uint32_t num_atoms;
	// 4 bytes pad
	real3 coords[];
} ALIGN_8 Coords;

Coords * FN(coords_malloc)(uint32_t num_atoms);
void FN(coords_copy)(real3 * dst, const real3 * src, uint64_t num_atoms);

#define coords_malloc(...) FN(coords_malloc)(__VA_ARGS__)
#define coords_copy(...) FN(coords_copy)(__VA_ARGS__)

// TODO: NEXTTIME: this REAL/FN stuff is getting out of hand... maybe switch to C++ ?


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


const Pos * FN(conf_space_pos)(const ConfSpace * conf_space, int posi);
const Coords * FN(conf_space_static_atoms)(const ConfSpace * conf_space);
const Conf * FN(conf_space_conf)(const ConfSpace * conf_space, const Pos * pos, int confi);
const Coords * FN(conf_space_conf_atoms)(const ConfSpace * conf_space, const Conf * conf);

#define conf_space_pos(...) FN(conf_space_pos)(__VA_ARGS__)
#define conf_space_static_atoms(...) FN(conf_space_static_atoms)(__VA_ARGS__)
#define conf_space_conf(...) FN(conf_space_conf)(__VA_ARGS__)
#define conf_space_conf_atoms(...) FN(conf_space_conf_atoms)(__VA_ARGS__)


extern const int32_t UNASSIGNED;
extern const int32_t STATIC_POS;


typedef struct AssignedCoords {
	const ConfSpace * conf_space;
	const int32_t * conf;
	Coords * atoms;
	real3 * pos_atoms[]; // indexed by posi
} AssignedCoords;


void FN(assigned_coords_free)(AssignedCoords * coords);
AssignedCoords * FN(conf_space_assign)(const ConfSpace * conf_space, const int32_t conf[]);
real3 * FN(assigned_coords_atoms)(AssignedCoords * assigned_coords, uint32_t posi);

#define assigned_coords_free(...) FN(assigned_coords_free)(__VA_ARGS__)
#define conf_space_assign(...) FN(conf_space_assign)(__VA_ARGS__)
#define assigned_coords_atoms(...) FN(assigned_coords_atoms)(__VA_ARGS__)

#endif // _CONFSPACE_

