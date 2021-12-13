include_directories ( $ENV{NUGEN_INC} )
cet_find_library( BXDECAY0 NAMES BxDecay0 PATHS ENV BXDECAY0_LIB NO_DEFAULT_PATH )

art_make(MODULE_LIBRARIES
         larsim_Simulation
         lardataalg_DetectorInfo
         lardataobj_Simulation
         nusimdata::SimulationBase
         larcoreobj_SummaryData
         larcorealg_Geometry
         larcore_Geometry_Geometry_service
         nurandom::RandomUtils_NuRandomService_service
         BXDECAY0
         art::Framework_Services_Registry
         art_root_io::tfile_support
         art_root_io::TFileService_service
         art::Persistency_Provenance
         messagefacility::MF_MessageLogger
         ROOT::Core
         ROOT::EG
         ROOT::Geom
        )

install_headers()
install_fhicl()
install_source()
