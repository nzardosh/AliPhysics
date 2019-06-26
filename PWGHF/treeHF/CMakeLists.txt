# **************************************************************************
# * Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
# *                                                                        *
# * Author: The ALICE Off-line Project.                                    *
# * Contributors are mentioned in the code where appropriate.              *
# *                                                                        *
# * Permission to use, copy, modify and distribute this software and its   *
# * documentation strictly for non-commercial purposes is hereby granted   *
# * without fee, provided that the above copyright notice appears in all   *
# * copies and that both the copyright notice and this permission notice   *
# * appear in the supporting documentation. The authors make no claims     *
# * about the suitability of this software for any purpose. It is          *
# * provided "as is" without express or implied warranty.                  *
# **************************************************************************/

#Module
set(MODULE PWGHFtreeHF)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliPhysics_SOURCE_DIR}/PWGHF/treeHF)

# Additional includes - alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
                    ${AliPhysics_SOURCE_DIR}/CORRFW
		    ${AliPhysics_SOURCE_DIR}/OADB
		    ${AliPhysics_SOURCE_DIR}/OADB/COMMON/MULTIPLICITY
                    ${AliPhysics_SOURCE_DIR}/PWGPP/EVCHAR/FlowVectorCorrections/QnCorrections
		    ${AliPhysics_SOURCE_DIR}/PWGPP/EVCHAR/FlowVectorCorrections/QnCorrectionsInterface/
                    ${AliPhysics_SOURCE_DIR}/PWG/EMCAL/EMCALbase
                    ${AliPhysics_SOURCE_DIR}/PWG/FLOW/Base
                    ${AliPhysics_SOURCE_DIR}/PWG/FLOW/Tasks
                    ${AliPhysics_SOURCE_DIR}/PWG/JETFW
                    ${AliPhysics_SOURCE_DIR}/PWG/muon
                    ${AliPhysics_SOURCE_DIR}/PWG/Tools
                    ${AliPhysics_SOURCE_DIR}/PWG/TRD
                    ${AliPhysics_SOURCE_DIR}/PWGHF/vertexingHF
                   # ${AliRoot_SOURCE_DIR}/ANALYSIS/ANALYSIS		
                   # ${AliRoot_SOURCE_DIR}/ANALYSIS/ANALYSISalice	
                    ${AliRoot_SOURCE_DIR}/EMCAL/EMCALbase
                    ${AliRoot_SOURCE_DIR}/EMCAL/EMCALUtils
		    ${AliPhysics_SOURCE_DIR}/PWGJE/EMCALJetTasks

  )

# Sources - alphabetical order
set(SRCS

  AliAnalysisTaskSEHFTreeCreator.cxx
  AliHFTreeHandler.cxx
  AliHFTreeHandlerD0toKpi.cxx
  AliHFTreeHandlerDplustoKpipi.cxx
  AliHFTreeHandlerDstoKKpi.cxx
  AliHFTreeHandlerLctopKpi.cxx
  AliHFTreeHandlerBplustoD0pi.cxx
  AliHFTreeHandlerDstartoKpipi.cxx
  AliHFTreeHandlerLc2V0bachelor.cxx
  AliJetTreeHandler.cxx
  AliParticleTreeHandler.cxx

  )

if(FASTJET_FOUND)


    include_directories(#${AliPhysics_SOURCE_DIR}/HLT/FJWrapper
                        ${AliPhysics_SOURCE_DIR}/JETAN/JETAN
                       )
    include_directories(SYSTEM ${FASTJET_INCLUDE_DIR})
    link_directories(${FASTJET_LIBS_DIR})
    add_definitions(${FASTJET_DEFINITIONS})

endif(FASTJET_FOUND)

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")
list (APPEND HDRS AliHFJetsUtils.h)

if(FASTJET_FOUND)
  list(APPEND HDRS ../../PWGJE/EMCALJetTasks/FJ_includes.h)
endif(FASTJET_FOUND)

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

set(ROOT_DEPENDENCIES Core EG Gpad Hist MathCore Physics Tree)
set(ALIROOT_DEPENDENCIES ANALYSIS ANALYSISalice AOD ESD PWGJETFW STEERBase PWGTools PWGHFvertexingHF)

# Add a shared library
add_library_tested(${MODULE} SHARED  ${SRCS} G__${MODULE}.cxx)

# Generate the ROOT map
# Dependecies
set(LIBDEPS ANALYSISalice PWGflowTasks PWGTRD PWGPPevcharQn PWGPPevcharQnInterface PWGHFvertexingHF PWGEMCALbase PWGJETFW ${ALIROOT_DEPENDENCIES} ${ROOT_DEPENDENCIES} PWGTools)
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

# Generate a PARfile target for this library
 add_target_parfile(${MODULE} "${SRCS}" "${HDRS}" "${MODULE}LinkDef.h" "${LIBDEPS}")

# Linking the library
# add_library_tested(${MODULE} SHARED $<TARGET_OBJECTS:${MODULE}-object>)
 target_link_libraries(${MODULE} ${LIBDEPS} ${ALIROOT_DEPENDENCIES} ${ROOT_DEPENDENCIES})

# target_include_directories(${MODULE}-object PUBLIC $<TARGET_PROPERTY:${MODULE},INCLUDE_DIRECTORIES>)
# set_target_properties(${MODULE}-object PROPERTIES COMPILE_DEFINITIONS $<TARGET_PROPERTY:${MODULE},COMPILE_DEFINITIONS>)


# set_target_properties(${MODULE}-object PROPERTIES COMPILE_FLAGS "${MODULE_COMPILE_FLAGS}")
# Setting the linking flags for the library
# set_target_properties(${MODULE} PROPERTIES LINK_FLAGS
#"${MODULE_LINK_FLAGS}")


# Public include folders that will be propagated to the dependencies
target_include_directories(${MODULE} PUBLIC ${incdirs})


 set(MODULE_COMPILE_FLAGS)
 set(MODULE_LINK_FLAGS)

 if(DATE_FOUND)
    set(MODULE_COMPILE_FLAGS "${DATE_CFLAGS}")
    set(MODULE_LINK_FLAGS "${DATE_LDFLAGS} ${DATE_LIBS}")
 endif(DATE_FOUND)

# Additional compilation and linking flags
 set(MODULE_COMPILE_FLAGS " ${MODULE_COMPILE_FLAGS}")


# System dependent: Modify the way the library is build
if(${CMAKE_SYSTEM} MATCHES Darwin)
    set_target_properties(${MODULE} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif(${CMAKE_SYSTEM} MATCHES Darwin)

# Installation
install(TARGETS ${MODULE}
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)
install(FILES ${HDRS} DESTINATION include)

install(DIRECTORY
            macros
        DESTINATION PWGHF/treeHF)

