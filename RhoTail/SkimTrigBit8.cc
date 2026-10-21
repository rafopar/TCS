/*
 * File:   SkimTrigBit8.cc
 * Author: rafopar
 *
 * Skims hipo files, keeping only events for which:
 *   - trigger bit 8 (counting from bit 0) is set in the "trigger" variable
 *     of the RUN::config bank, AND
 *   - there is at least one REC::Particle with pid == 211 and at least one
 *     REC::Particle with pid == -211.
 *
 * Usage: SkimTrigBit8.exe <input1.hipo> [input2.hipo ...]
 *
 * For each input file the run number and the file index range are deduced
 * from its name, e.g. for
 *   rec_clas_005163.evio.00055-00059.hipo
 * the run is 5163, the low file index is 55 and the high file index is 59.
 * The skimmed output is written to
 * OutSkims/RUN/PionPairSkim_RUN_FILE_LOW_FILE_HIGH.hipo, e.g.
 * OutSkims/005163/PionPairSkim_005163_00055_00059.hipo
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

    const int TRIGGER_BIT = 8; // Counting from bit 0
    const int PID_PIP = 211;
    const int PID_PIM = -211;
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
}

int main(int argc, char** argv) {

    if (argc < 2) {
        cout << "Usage: " << argv[0] << " <input1.hipo> [input2.hipo ...]" << endl;
        cout << "Keeps only events with trigger bit " << TRIGGER_BIT << " set (RUN::config.trigger)" << endl;
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

        const std::string outputFile = outDir + "/PionPairSkim_" + run + "_" + fileLow + "_" + fileHigh + ".hipo";

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
        hipo::bank bRunConfig(factory.getSchema("RUN::config"));
        hipo::bank bRecPart(factory.getSchema("REC::Particle"));

        long evCounter = 0;
        long nKept = 0;

        while (reader.next() == true) {
            reader.read(event);
            evCounter = evCounter + 1;

            event.getStructure(bRunConfig);

            if (bRunConfig.getRows() == 0) {
                continue;
            }

            long trigger = bRunConfig.getLong("trigger", 0);
            bool isBit8Set = ((trigger >> TRIGGER_BIT) & 1L) != 0;

            if (!isBit8Set) {
                continue;
            }

            event.getStructure(bRecPart);

            bool hasPiPlus = false;
            bool hasPiMinus = false;
            int nPart = bRecPart.getRows();

            for (int iPart = 0; iPart < nPart && !(hasPiPlus && hasPiMinus); iPart++) {
                int pid = bRecPart.getInt("pid", iPart);

                if (pid == PID_PIP) {
                    hasPiPlus = true;
                } else if (pid == PID_PIM) {
                    hasPiMinus = true;
                }
            }

            if (hasPiPlus && hasPiMinus) {
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
