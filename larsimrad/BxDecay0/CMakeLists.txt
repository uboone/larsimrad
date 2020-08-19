include_directories ( ENV BXDECAY0_INC )
cet_find_library( BXDECAY0 NAMES BxDecay0 PATHS ENV BXDECAY0_LIB NO_DEFAULT_PATH )

art_make(MODULE_LIBRARIES
         larsim_Simulation
         lardataalg_DetectorInfo
         lardataobj_Simulation
         nusimdata_SimulationBase
         larcoreobj_SummaryData
         larcorealg_Geometry
         larcore_Geometry_Geometry_service
         nurandom_RandomUtils_NuRandomService_service
         ${BXDECAY0}
         ${ART_FRAMEWORK_SERVICES_REGISTRY}
         ${ART_ROOT_IO_TFILE_SUPPORT}
         ${ART_ROOT_IO_TFILESERVICE_SERVICE}
         art_Persistency_Provenance
         ${MF_MESSAGELOGGER}
         ${CLHEP}
         ROOT::Core
         ROOT::EG
         ROOT::Geom
         ${ROOT_LIBRARIES}
        )

#         ${BxDecay0_LIBRARIES}
install_headers()
install_fhicl()
install_source()
