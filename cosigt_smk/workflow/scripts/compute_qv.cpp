#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "edlib.h"

std::string readFasta(const std::string& filename) {
    std::ifstream file(filename);
    std::string line, sequence;
    while (std::getline(file, line)) {
        if (!line.empty() && line[0] != '>') {
            sequence += line;
        }
    }
    return sequence;
}

int main(int argc, char* argv[]) {

    std::string seq1 = readFasta(argv[1]);
    std::string seq2 = readFasta(argv[2]);

    EdlibAlignResult result = edlibAlign(
        seq1.c_str(), seq1.size(),
        seq2.c_str(), seq2.size(),
        edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0)
    );

    int delta = result.editDistance;
    int alignment_length = result.alignmentLength;

    double delta_max = delta < 0.5 ? 0.5 : delta;
    double ratio = delta_max / alignment_length;
    double qv = -10.0 * log10(ratio);

    std::cout << qv << std::endl;

    edlibFreeAlignResult(result);
    return 0;
}

