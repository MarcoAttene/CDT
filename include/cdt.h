#pragma once

#include "PLC.h"

#pragma intrinsic(fabs)

#include "inputPLC.h"
#include "logger.h"

// createSteinerCDT
// 
// 'plc' is a valid input PLC to the process. Validity is assumed but not verified!
// 'options' is a (possibly empty) string of characters, each controlling
// one option as follows:
// l: log results to cdt_log.csv
// b: add eight vertices to enclose everything in a box
// v: verbose mode
// f: try to make the output representable using floating point
// w: log to screen

TetMesh* createSteinerCDT(inputPLC& plc, const char* options) {
	bool log = false, bbox = false, verbose = false, snap = false, logscreen = false;
	//bool optimize = false;

	for (int i = 0; i < strlen(options); i++) switch (options[i]) {
	case 'l':
		log = true; break;
	case 'b':
		bbox = true; break;
	case 'v':
		verbose = true; break;
	case 'w':
		logscreen = true; break;
		//case 'f':
		//	snap = true; break;
		//case 'o':
		//	optimize = true; break;
	} // Just ignore unknown options

	if (bbox) plc.addBoundingBoxVertices();

	if (logscreen) {
		log = true;
		startLogging(NULL);
	}
	else if (log) startLogging(plc.input_file_name);

	// Build a delaunay tetrahedrization of the vertices
	TetMesh* tin = new TetMesh;
	tin->init_vertices(plc.coordinates.data(), plc.numVertices());
	tin->tetrahedrize();

	if (verbose) printf("DT of the vertices built\n");

	if (log) logTimeChunk();

	// Build a structured PLC linked to the Delaunay tetrahedrization
	PLCx Steiner_plc(*tin, plc.triangle_vertices.data(), plc.numTriangles());

	// Recover segments by inserting Steiner points in both the PLC and the tetrahedrization
	Steiner_plc.segmentRecovery_HSi(!verbose);

	if (log) logTimeChunk();

	// Recover PLC faces by locally remeshing the tetrahedrization
	bool sisMethodWorks = Steiner_plc.faceRecovery(!verbose);

	if (log) logTimeChunk();

	// Mark the tets which are bounded by the PLC.
	// If the PLC is not a valid polyhedron (i.e. it has odd-valency edges)
	// all the tets but the ghosts are marked as "internal".
	uint32_t num_inner_tets = (uint32_t)Steiner_plc.markInnerTets();

	if (log) logTimeChunk();

	if (log) {
		logMemInfo();
		logBoolean(Steiner_plc.is_polyhedron);
		logInteger(plc.numVertices());
		logInteger(Steiner_plc.input_nt);
		logInteger(Steiner_plc.numSteinerVertices());
		logInteger(tin->countNonGhostTets());
		logInteger(num_inner_tets);
		size_t nflip, nflat;
		tin->hasBadSnappedOrientations(nflip, nflat);
		logInteger((uint32_t)nflat);
		logInteger((uint32_t)nflip);
		logBoolean(sisMethodWorks);
		finishLogging();
	}

	return tin;
}
