#ifdef _MSC_VER // Workaround for known bug on MSVC
#define _HAS_STD_BYTE 0  // https://developercommunity.visualstudio.com/t/error-c2872-byte-ambiguous-symbol/93889
#endif

#include <iostream>
#include <fstream>
#include "delaunay.h"
#include "inputPLC.h"
#include "PLC.h"

using namespace std;

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

TetMesh* createSteinerCDT(inputPLC& plc, const char *options) {
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
	case 'f':
		snap = true; break;
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
	TetMesh  *tin = new TetMesh;
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

	if (snap) {
		if (!tin->optimizeNearDegenerateTets(verbose)) {
			std::cerr << "Could not force FP representability.\n";
		}
	}

	//if (optimize) tin->optimizeMesh();

	return tin;
}

// saveOutputFile
// 
// 'tin' is a characterized tet mesh produced by the function above.
// 'filename' is the name of the output file without extension.
// The file produced will be called 'filename.tet' and/or
// 'filename.off' (if 's' option is used).
// 'options' is a (possibly empty) string of characters, each controlling
// one option as follows:
// q: rational output
// n: binary output
// r: remove outer tetrahedra from output (if input is closed)
// s: saves skin to an ASCII OFF file (triangles between IN and OUT)
// m: saves mesh to MEDIT format instead of TET

bool saveOutputFile(TetMesh& tin, const char* filename, const char* options) {
	bool rational = false, binary = false, erode = false, skin = false, medit = false;
	for (int i = 0; i < strlen(options); i++) switch (options[i]) {
	case 'q':
		rational = true; break;
	case 'n':
		binary = true; break;
	case 'r':
		erode = true; break;
	case 's':
		skin = true; break;
	case 'm':
		medit = true; break;
	}

	char tetfilename[2048], offfilename[2048];

	bool ret = true;

	if (medit) {
		if (rational || binary) {
			std::cerr << "Rational and binary modes are not available when saving to MEDIT format\n";
			ret = false;
		}
		else {
			sprintf(tetfilename, "%s.mesh", filename);
			ret &= tin.saveMEDIT(tetfilename, erode);
		}
	}
	else {
		sprintf(tetfilename, "%s.tet", filename);
		if (!rational && !binary) ret &= tin.saveTET(tetfilename, erode);
		if (!rational && binary) ret &= tin.saveBinaryTET(tetfilename, erode);
		if (rational && !binary) ret &= tin.saveRationalTET(tetfilename, erode);
		if (rational && binary) {
			std::cerr << "Save to rational is supported only in ASCII mode\n";
			ret = false;
		}
	}

	if (skin) {
		sprintf(offfilename, "%s.off", filename);
		ret &= tin.saveBoundaryToOFF(offfilename);
	}

	return ret;
}

int main(int argc, char* argv[])
{
	initFPU();

	if (argc < 2) {
		std::cout << "CDT - Create a constrained Delaunay tetrahedrization out of a triangulated OFF file.\n";
		std::cout << "USAGE: CDT [-lbvfqnrs] filename.off\n";
		std::cout << "Example 1: CDT -bv test.off\n";
		std::cout << "Example 2: CDT -b -v test.off\n";
		std::cout << "OPTIONS:\n";
		std::cout << "-l: log results to cdt_log.csv\n";
		std::cout << "-b: add eight vertices to enclose everything in a box\n";
		std::cout << "-v: verbose mode\n";
		std::cout << "-w: log on screen instead of file (implies -l)\n";
		std::cout << "-f: try to make the output representable using floating point\n";
		std::cout << "-q: rational output\n";
		std::cout << "-n: binary output\n";
		std::cout << "-m: use MEDIT format instead of TET\n";
		std::cout << "-r: remove outer tetrahedra from output (if input is closed)\n";
		std::cout << "-s: saves skin to an ASCII OFF file (triangles between IN and OUT)\n";
		std::cout << "OUTPUT:\n";
		std::cout << "Output has same name (and path) as input with an extension appended.\n";
		std::cout << "E.g. CDT my_dir/test.off produces my_dir/test.off.tet\n";
		std::cout << "E.g. CDT -s my_dir/test.off produces my_dir/test.off.tet and my_dir/test.off.off\n";
		return 0;
	}

	char filename[2048] = "..\\Input_file\\bracket.off";

	std::string options = "";

	for (int i = 1; i < argc; i++)
		if (argv[i][0] == '-') {
			for (int j = 1; j < strlen(argv[i]); j++) options += argv[i][j];
		}
		else memcpy(filename, argv[i], strlen(argv[i])+1);

	// Load a valid PLC from file
	inputPLC plc;
	plc.initFromFile(filename, options.find('v') != std::string::npos);

	TetMesh* tin = createSteinerCDT(plc, options.c_str());

	if (saveOutputFile(*tin, filename, options.c_str()))
		printf("Finished\n");

	return 0;
}
