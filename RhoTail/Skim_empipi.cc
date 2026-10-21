/*
 * File:   Skim_empipi.cc
 * Author: rafopar
 *
 * Skims hipo files, keeping only events for which REC::Particle has:
 *   - at least one e-  (pid ==  11) with abs(status) >= 2000 && abs(status) < 4000, AND
 *   - at least one pi- (pid == -211) with abs(status) >= 2000 && abs(status) < 4000, AND
 *   - at least one pi+ (pid ==  211) with abs(status) >= 2000 && abs(status) < 4000.
 *
 * Usage: Skim_empipi.exe <input1.hipo> [input2.hipo ...]
 *
 * For each input file the run number and the file index range are deduced
 * from its name, e.g. for
 *   rec_clas_005163.evio.00055-00059.hipo
 * the run is 5163, the low file index is 55 and the high file index is 59.
 * The skimmed output is written to
 * OutSkims/RUN/empipiSkim_RUN_FILE_LOW_FILE_HIGH.hipo, e.g.
 * OutSkims/005163/empipiSkim_005163_00055_00059.hipo
 */

#include <cstdlib>
#include <iostream>
#include <string>
#include <regex>

#include <TSystem.h>

// ===== Hipo headers =====
#include <reader.h>
#include <writer.h>
#include <dictionary.h>

using namespace std;

namespace {

    const int PID_EM = 11;
    const int PID_PIP = 211;
    const int PID_PIM = -211;
    const int STATUS_MIN = 2000;
    const int STATUS_MAX = 4000;
    const char* OUT_DIR = "OutSkims";

    // Deduces the run number and the low/high file indices from an input
    // hipo file name such as rec_clas_005163.evio.00055-00059.hipo
    bool parseRunAndFileRange(const std::string& inputFile, std::string& run, std::string& fileLow, std::string& fileHigh) {

        std::string base = inputFile;
        size_t slashPos = base.find_last_of('/');
        if (slashPos != std::string::npos) {
            base = base.substr(slashPos + 1);
        }

        static const std::regex pattern(R"((\d+)\.evio\.(\d+)-(\d+)\.hipo$)");
        std::smatch match;

        if (!std::regex_search(base, match, pattern)) {
            return false;
        }

        run = match[1];
        fileLow = match[2];
        fileHigh = match[3];

        return true;
    }

    bool isGoodStatus(int status) {
        int absStatus = std::abs(status);
        return absStatus >= STATUS_MIN && absStatus < STATUS_MAX;
    }
}

int main(int argc, char** argv) {

    if (argc < 2) {
        cout << "Usage: " << argv[0] << " <input1.hipo> [input2.hipo ...]" << endl;
        cout << "Keeps only events with a good-status e-, pi- and pi+ in REC::Particle" << endl;
        exit(1);
    }

    for (int iFile = 1; iFile < argc; iFile++) {

        const std::string inputFile = argv[iFile];

        std::string run, fileLow, fileHigh;
        if (!parseRunAndFileRange(inputFile, run, fileLow, fileHigh)) {
            cerr << "Could not deduce the run/file range from the file name \"" << inputFile << "\", skipping." << endl;
            continue;
        }

        const std::string outDir = std::string(OUT_DIR) + "/" + run;
        gSystem->Exec(Form("mkdir -p %s", outDir.c_str()));

        const std::string outputFile = outDir + "/empipiSkim_" + run + "_" + fileLow + "_" + fileHigh + ".hipo";

        cout << "Skimming " << inputFile << " -> " << outputFile << endl;

        hipo::reader reader;
        reader.open(inputFile.c_str());

        hipo::dictionary factory;
        reader.readDictionary(factory);

        hipo::writer writer;
        for (const auto& schemaName : factory.getSchemaList()) {
            writer.getDictionary().addSchema(factory.getSchema(schemaName.c_str()));
        }
        writer.open(outputFile.c_str());

        hipo::event event;
        hipo::bank bRecPart(factory.getSchema("REC::Particle"));

        long evCounter = 0;
        long nKept = 0;

        while (reader.next() == true) {
            reader.read(event);
            evCounter = evCounter + 1;

            event.getStructure(bRecPart);

            bool hasEm = false;
            bool hasPiPlus = false;
            bool hasPiMinus = false;
            int nPart = bRecPart.getRows();

            for (int iPart = 0; iPart < nPart && !(hasEm && hasPiPlus && hasPiMinus); iPart++) {
                int pid = bRecPart.getInt("pid", iPart);
                int status = bRecPart.getInt("status", iPart);

                if (!isGoodStatus(status)) {
                    continue;
                }

                if (pid == PID_EM) {
                    hasEm = true;
                } else if (pid == PID_PIP) {
                    hasPiPlus = true;
                } else if (pid == PID_PIM) {
                    hasPiMinus = true;
                }
            }

            if (hasEm && hasPiPlus && hasPiMinus) {
                writer.addEvent(event);
                nKept = nKept + 1;
            }

            if (evCounter % 50000 == 0) {
                cout.flush() << "Processed " << evCounter << " events, kept " << nKept << " \r";
            }
        }

        writer.close();

        cout << endl << "Done. Processed " << evCounter << " events, kept " << nKept
                << " -> " << outputFile << endl;
    }

    return 0;
}
