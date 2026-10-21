set(CMAKE_VERBOSE_MAKEFILE ON)

# clas12AnaTools' RecEvent/RecPart classes embed hipo::bank objects by value and
# were built against clas12AnaTools's own vendored hipo4 checkout (hipo/hipo4),
# whose hipo::bank/structure layout has diverged from the standalone
# /home/rafopar/work/git/hipo checkout used by the plain hipo-reading targets
# below. Any target using RecEvent/RecPart must build against that same
# vendored hipo4/libHipo4.a, or the ABI mismatch corrupts memory at runtime
# (crashes in hipo::bank's destructor).
set(CLAS12ANATOOLS_HIPO4_INCLUDE_DIR /home/rafopar/work/git/clas12AnaTools/hipo/hipo4)
find_library(CLAS12ANATOOLS_HIPO4_LIBRARY NAMES Hipo4 PATHS /home/rafopar/work/git/clas12AnaTools/lib)

################################  SkimTrigBit8  ##################################################
add_executable(SkimTrigBit8 SkimTrigBit8.cc)
set_target_properties(SkimTrigBit8 PROPERTIES OUTPUT_NAME "SkimTrigBit8.exe")
target_include_directories(SkimTrigBit8 PRIVATE ${ROOT_INCLUDE_DIRS} ${Hipo4_INCLUDE_DIRS})
target_link_libraries(SkimTrigBit8 ${Hipo4_LIBRARY} ${LZ4_LIBRARY} ${ROOT_LIBRARIES})
###################################################################################################

################################  Skim_empipi  ####################################################
add_executable(Skim_empipi Skim_empipi.cc)
set_target_properties(Skim_empipi PROPERTIES OUTPUT_NAME "Skim_empipi.exe")
target_include_directories(Skim_empipi PRIVATE ${ROOT_INCLUDE_DIRS} ${Hipo4_INCLUDE_DIRS})
target_link_libraries(Skim_empipi ${Hipo4_LIBRARY} ${LZ4_LIBRARY} ${ROOT_LIBRARIES})
###################################################################################################

################################  AnaPiPiProt  ####################################################
add_executable(AnaPiPiProt AnaPiPiProt.cc)
set_target_properties(AnaPiPiProt PROPERTIES OUTPUT_NAME "AnaPiPiProt.exe")
target_include_directories(AnaPiPiProt PRIVATE ${ROOT_INCLUDE_DIRS} ${CLAS12ANATOOLS_HIPO4_INCLUDE_DIR} ${CLAS12AnaTools_INCLUDE_DIRS})
target_link_libraries(AnaPiPiProt ${CLAS12ANATOOLS_HIPO4_LIBRARY} ${CLAS12AnaTools_LIBRARY} ${LZ4_LIBRARY} ${ROOT_LIBRARIES})
###################################################################################################

################################  CompareHists  ###################################################
add_executable(CompareHists CompareHists.cc)
set_target_properties(CompareHists PROPERTIES OUTPUT_NAME "CompareHists.exe")
target_include_directories(CompareHists PRIVATE ${ROOT_INCLUDE_DIRS})
target_link_libraries(CompareHists ${ROOT_LIBRARIES})
###################################################################################################

################################ Installing Executable targets #############################
install( TARGETS SkimTrigBit8 Skim_empipi AnaPiPiProt CompareHists
        DESTINATION ${CMAKE_INSTALL_PREFIX}/RhoTail/bin )

################################ Installing driver scripts #################################
install( FILES RunSkimTrigBit8.py RunAllRuns.py RunAnaPiPiProt.py RunAnaPiPiProt_MC.py
        DESTINATION ${CMAKE_INSTALL_PREFIX}/RhoTail
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                    GROUP_READ GROUP_EXECUTE
                    WORLD_READ WORLD_EXECUTE)

################################ Installing run lists #######################################
install( FILES F18_in_RunList.dat
        DESTINATION ${CMAKE_INSTALL_PREFIX}/RhoTail)

################################ Installing ROOT macros #####################################
# RhoFitGUI.cc is an interactive ROOT GUI macro (no main()), meant to be run as
# `root RhoFitGUI.cc` or `root RhoFitGUI.cc+`, not compiled by CMake.
install( FILES RhoFitGUI.cc
        DESTINATION ${CMAKE_INSTALL_PREFIX}/RhoTail)

################################ Installing needed directores ##############################
install( DIRECTORY
        DESTINATION ${CMAKE_INSTALL_PREFIX}/RhoTail/OutSkims )
install( DIRECTORY
        DESTINATION ${CMAKE_INSTALL_PREFIX}/RhoTail/MCin )
install( DIRECTORY
        DESTINATION ${CMAKE_INSTALL_PREFIX}/RhoTail/Hists )
install( DIRECTORY
        DESTINATION ${CMAKE_INSTALL_PREFIX}/RhoTail/Figs )
