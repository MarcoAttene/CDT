#include <chrono>

FILE* log_fp;
std::chrono::steady_clock::time_point time_point;

inline void startLogging(const char* fn) {
    if (fn != NULL) {
        log_fp = fopen("cdt_log.csv", "r");
        if (log_fp == NULL) {
            log_fp = fopen("cdt_log.csv", "w");
            fprintf(log_fp, "Input_File, Delaunay(ms), Segment_Rec(ms), Face_Rec(ms), "
                "Erosion(ms), Memory(MB), Is_Polyhedron, Input_Verts, Input_Tris, Steiner, Tot_Tets, In_Tets, Flat_Tets, Flipped_Tets, Si's method works\n");
        }
        else {
            fclose(log_fp);
            log_fp = fopen("cdt_log.csv", "a");
        }
        if (log_fp == NULL) ip_error("Can't open cdt_log.csv for logging!\n");

        size_t i;
        for (i = strlen(fn); i > 0; i--) if (fn[i - 1] == '\\' || fn[i - 1] == '/') break;
        fprintf(log_fp, "%s", fn + i);
    }
    else {
        log_fp = stdout;
        fprintf(log_fp, " , Delaunay(ms), Segment_Rec(ms), Face_Rec(ms), "
            "Erosion(ms), Memory(MB), Is_Polyhedron, Input_Verts, Input_Tris, Steiner, Tot_Tets, In_Tets, Flat_Tets, Flipped_Tets, Si's method works\n");
    }

    time_point = std::chrono::steady_clock::now();
}

inline void logTimeChunk() {
    std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
    uint64_t ms = std::chrono::duration_cast<std::chrono::milliseconds>(now - time_point).count();
    time_point = now;

    fprintf(log_fp, ", %zu", ms);
}

inline void logBoolean(bool b) {
    fprintf(log_fp, ", %s", b ? "True" : "False");
}

inline void logInteger(uint32_t n) {
    fprintf(log_fp, ", %u", n);
}

inline void finishLogging() {
    fprintf(log_fp, "\n");
    if (log_fp != stdout) fclose(log_fp);
}

#ifdef _MSC_VER
#include <windows.h>
#include <psapi.h>

// To ensure correct resolution of symbols, add Psapi.lib to TARGETLIBS
// and compile with -DPSAPI_VERSION=1

double getPeakMegabytesUsed()
{
    HANDLE hProcess = OpenProcess(PROCESS_QUERY_INFORMATION | PROCESS_VM_READ, FALSE, GetCurrentProcessId());
    if (NULL == hProcess) return 0;

    PROCESS_MEMORY_COUNTERS pmc;
    double mem = 0;
    if (GetProcessMemoryInfo(hProcess, &pmc, sizeof(pmc)))
    {
        mem = pmc.PeakWorkingSetSize / 1048576.0;
    }

    CloseHandle(hProcess);
    return mem;
}

// Mem info in Mb
inline void logMemInfo()
{
    fprintf(log_fp, ", %.2f", getPeakMegabytesUsed());
}
#else

#include <sys/time.h>
#include <sys/resource.h>

inline void logMemInfo() {
    struct rusage r_usage;
    getrusage(RUSAGE_SELF, &r_usage);
    fprintf(log_fp, ", %.2f", r_usage.ru_maxrss / 1000.0);
}
#endif
