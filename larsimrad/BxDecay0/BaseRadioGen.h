/// \file  BaseRadioGenHelper.cc
/// \brief Base module for radiogen/decay0gen for sampling etc
/// \author  plasorak@FNAL.GOV
///          June 2020 PLasorak
#pragma once

// C++ includes.
#include <string>
#include <regex>
#include <cmath>
#include <math.h>
#include <memory>
#include <iterator>
#include <sys/stat.h>
#include <TGeoManager.h>
#include <TGeoMaterial.h>
#include <TGeoNode.h>

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "cetlib/search_path.h"
#include "cetlib/exempt_ptr.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
// nurandom includes
#include "nurandom/RandomUtils/NuRandomService.h"

// nusimdata, nugen includes
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nugen/EventGeneratorBase/evgenbase.h"

// lar includes
#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SummaryData/RunData.h"

// root includes

#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TVector3.h"
#include "TRandom.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoisson.h"



namespace simb { class MCTruth; }
namespace evgen { class clhep_random; }

namespace evgen {
  class BaseRadioGen: public art::EDProducer {
  public:
    explicit BaseRadioGen(fhicl::ParameterSet const& pset);
    void produce(art::Event& evt);
    void beginRun(art::Run& run);
    void beginJob();
    void endJob();

    virtual void produce_radio(art::Event& evt) = 0;
    virtual void beginRun_radio(art::Run&) {};
    virtual void beginJob_radio() {}
    virtual void endJob_radio() {}

  protected:
    int GetNDecays();
    bool GetGoodPositionTime(TLorentzVector& position);
    TLorentzVector dirCalc(double p, double m);

    void FillHistos(simb::MCParticle& part);
    std::vector<std::string> m_isotope; ///< isotope to simulate.  Example:  "Ar39"

    CLHEP::HepRandomEngine& GetRandomEngine() { return m_engine; }
    int GetNEvents() { return m_nevent; }
    void GetTs(double& T0, double& T1) {T0=m_T0; T1=m_T1;}
    void GetXs(double& X0, double& X1) {X0=m_X0; X1=m_X1;}
    void GetYs(double& Y0, double& Y1) {Y0=m_Y0; Y1=m_Y1;}
    void GetZs(double& Z0, double& Z1) {Z0=m_Z0; Z1=m_Z1;}

  private:
    void CalculateActiveVolumeFromAllNodes();
    void CalculateActiveVolumeFromXYZ();
    void SimplePDG(int pdg, int& simple, std::string& name);
    void DeclareOutputHistos();
    std::set<const TGeoNode*> m_all_nodes;
    void FillAllNodes(const TGeoNode* curnode);
    void GetNodeXYZMinMax(const TGeoNode* from, const TGeoNode* to,
                          double& x_min, double &x_max, double& y_min, double& y_max, double& z_min, double& z_max);
    bool IsDaughterNode(const TGeoNode* daughter_node, const TGeoNode* mother_node);
    bool findNode(const TGeoNode* curnode, std::string& tgtnname,
                  const TGeoNode* & targetnode);
    bool findMotherNode(const TGeoNode* cur_node, std::string& daughter_name,
                        const TGeoNode* & mother_node);

    std::string m_material;    ///< regex of materials in which to generate the decays.  Example: "LAr"
    std::string m_volume_rand; ///< The volume in which to generate the decays
    std::string m_volume_gen;  ///< The volume in which to generate the decays
    double      m_Bq;          ///< Radioactivity in Becquerels (decay per sec) per cubic cm.
    double      m_rate;        ///< Radioactivity in Becquerels (decay per sec) use either of this of Bq
    double      m_T0;          ///< Beginning of time window to simulate in ns
    double      m_T1;          ///< End of time window to simulate in ns
    double      m_X0;          ///< Bottom corner x position (cm) in world coordinates
    double      m_Y0;          ///< Bottom corner y position (cm) in world coordinates
    double      m_Z0;          ///< Bottom corner z position (cm) in world coordinates
    double      m_X1;          ///< Top corner x position (cm) in world coordinates
    double      m_Y1;          ///< Top corner y position (cm) in world coordinates
    double      m_Z1;          ///< Top corner z position (cm) in world coordinates

    art::ServiceHandle<geo::Geometry const> m_geo_service;
    CLHEP::HepRandomEngine& m_engine;

    TGeoManager* m_geo_manager;
    double m_volume_cc;

    std::unique_ptr<CLHEP::RandFlat   > m_random_flat;
    std::unique_ptr<CLHEP::RandPoisson> m_random_poisson;

    bool m_geo_volume_mode;
    bool m_rate_mode;
    bool m_volume_rand_present;
    
    size_t m_max_tries_event;
    size_t m_max_tries_rate_calculation;
    size_t m_target_n_point_rate_calculation;
    
    std::regex  m_regex_material;
    std::regex  m_regex_volume;
    std::map<const TGeoNode*,double> m_good_nodes     = {};
    std::vector<const TGeoVolume  *> m_good_volumes   = {};
    std::vector<const TGeoMaterial*> m_good_materials = {};

    
    bool  m_flat_distrib_xpos;
    bool  m_flat_distrib_ypos;
    bool  m_flat_distrib_zpos;
    TF1* m_distrib_xpos;
    TF1* m_distrib_ypos;
    TF1* m_distrib_zpos;

    int m_nevent;
    std::map<int,TH2D*> m_pos_xy_TH2D;
    std::map<int,TH2D*> m_pos_xz_TH2D;
    std::map<int,TH1D*> m_dir_x_TH1D;
    std::map<int,TH1D*> m_dir_y_TH1D;
    std::map<int,TH1D*> m_dir_z_TH1D;
    std::map<int,TH1D*> m_mom_TH1D;
    std::map<int,TH1D*> m_ke_TH1D;
    std::map<int,TH1D*> m_time_TH1D;
    TH1D* m_pdg_TH1D;
  };

  
  void BaseRadioGen::GetNodeXYZMinMax(const TGeoNode* from, const TGeoNode* to,
                                      double& x_min, double &x_max,
                                      double& y_min, double& y_max,
                                      double& z_min, double& z_max) {
    
    std::vector<const TGeoNode*> mother_nodes;
    const TGeoNode* current_node=from;
    std::string daughter_name = from->GetName();
    int nmax = 20;
    int iter=0;
    while (current_node != to and iter++<nmax) {
      const TGeoNode* mother_node = nullptr;
      daughter_name =current_node->GetName();
      bool found_mum = findMotherNode(to, daughter_name, mother_node);
      if(not found_mum) {
        throw cet::exception("BaseRadioGen") << "Didn't find the mum of the following node: " << daughter_name;
      }
      mother_nodes.push_back(mother_node);
      current_node = mother_node;
    }


    TGeoVolume* vol   = from->GetVolume();
    TGeoShape*  shape = vol->GetShape();
    TGeoBBox*   bbox  = (TGeoBBox*)shape;

    double dx = bbox->GetDX();
    double dy = bbox->GetDY();
    double dz = bbox->GetDZ();

    double halfs [3] = { dx, dy, dz };
    double posmin[3] = {  1.0e30,  1.0e30,  1.0e30 };
    double posmax[3] = { -1.0e30, -1.0e30, -1.0e30 };

    const double* origin = bbox->GetOrigin();
    for ( int ix = -1; ix <= 1; ix += 2) {
      for ( int iy = -1; iy <= 1; iy += 2) {
        for ( int iz = -1; iz <= 1; iz += 2) {
          double local[3];
          local[0] = origin[0] + (double)ix*halfs[0];
          local[1] = origin[1] + (double)iy*halfs[1];
          local[2] = origin[2] + (double)iz*halfs[2];
          double master[3];
          from->LocalToMaster(local,master);
          for (auto const& mum: mother_nodes) {
            local[0] = master[0];
            local[1] = master[1];
            local[2] = master[2];
            mum->LocalToMaster(local, master);
          }
          for ( int j = 0; j < 3; ++j ) {
            posmin[j] = TMath::Min(posmin[j],master[j]);
            posmax[j] = TMath::Max(posmax[j],master[j]);
          }
        }
      }
    }
      
    x_min = posmin[0];
    y_min = posmin[1];
    z_min = posmin[2];
    x_max = posmax[0];
    y_max = posmax[1];
    z_max = posmax[2];
  }
  
  BaseRadioGen::BaseRadioGen(fhicl::ParameterSet const& pset):
    EDProducer(pset),
    m_engine(art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "BaseRadioGen", pset, "SeedBaseRadioGen")) {

    m_nevent=0;

    produces<std::vector<simb::MCTruth>>();
    produces<sumdata::RunData, art::InRun>();

    m_max_tries_event = pset.get<size_t>("max_tries_event", 1'000'000);
    m_max_tries_rate_calculation = pset.get<size_t>("max_tries_rate_calculation", 40'000'000);
    m_target_n_point_rate_calculation = pset.get<size_t>("target_n_point_rate_calculation", 100'000);
    
    m_material = pset.get<std::string>("material", ".*");
    m_regex_material = (std::regex)m_material;

    m_volume_gen = pset.get<std::string>("volume_gen", ".*");
    m_regex_volume = (std::regex)m_volume_gen;

    m_rate_mode = pset.get_if_present<double>("rate", m_rate);
    if (not m_rate_mode)
      m_Bq = pset.get<double>("BqPercc");

    bool timed_mode = pset.get_if_present<double>("T0", m_T0);
    if (timed_mode) {
      m_T1 = pset.get<double>("T1");
    } else {
      auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
      auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob(clockData);
      int nsample = detProp.NumberTimeSamples();
      m_T0 = -nsample * sampling_rate(clockData);
      m_T1 = -m_T0;
    }

    m_geo_manager = m_geo_service->ROOTGeoManager();

    m_random_flat    = std::make_unique<CLHEP::RandFlat   >(m_engine);
    m_random_poisson = std::make_unique<CLHEP::RandPoisson>(m_engine);

    if (m_material != ".*" || m_volume_gen != ".*") {
      MF_LOG_INFO("BaseRadioGen") << "Calculating the volume of " << m_material << " and the volume " << m_volume_gen << " in the geometry.\n";
      
      double xyz[3];

      TObjArray* volumes   = m_geo_manager->GetListOfVolumes();
      TList    * materials = m_geo_manager->GetListOfMaterials();
      
      for (int i=0; i<volumes->GetEntries(); ++i) {
        TGeoVolume* volume = (TGeoVolume*)volumes->At(i);
        std::string volume_name = volume->GetName();
        bool good = std::regex_match(volume_name, m_regex_volume);
        if (good) {
          
          if (m_volume_gen != ".*") {
            MF_LOG_INFO("BaseRadioGen") << "  Volume " << volume_name << " is an accepted volume of " << volume->GetShape()->Capacity() << "cm^3.\n";
            if (volume->GetShape()->TestShapeBits(TGeoShape::kGeoBox)) {
              TGeoBBox* box = dynamic_cast<TGeoBBox*>(volume->GetShape());
              if (box)
                MF_LOG_INFO("BaseRadioGen") << "  It is a box of size " << 2.*box->GetDX() << " x " << 2.*box->GetDY() << " x " << 2.*box->GetDZ() << " cm^3\n";
            }
          }
          
          m_good_volumes.push_back(volume);
        }
      }
      
      for (int i=0; i<materials->GetEntries(); ++i) {
        TGeoMaterial* material = (TGeoMaterial*)materials->At(i);
        std::string material_name = material->GetName();
        bool good = std::regex_match(material_name, m_regex_material);
        if (good) m_good_materials.push_back(material);
      }
    }
    MF_LOG_INFO("BaseRadioGen") << m_good_volumes  .size() << " volumes correspond to the regex \""   << m_volume_gen << "\".\n";
    MF_LOG_INFO("BaseRadioGen") << m_good_materials.size() << " materials correspond to the regex \"" << m_material   << "\".\n";
    
    double dummy;
    if (pset.get_if_present<double>("X0", dummy) or
        pset.get_if_present<double>("X1", dummy) or
        pset.get_if_present<double>("Y0", dummy) or
        pset.get_if_present<double>("Y1", dummy) or
        pset.get_if_present<double>("Z0", dummy) or
        pset.get_if_present<double>("Z1", dummy)) {
      /// If we specified X, Y, Z we will throw with what the user provided
      /// This will throw an error if somebody specifies only one X0.
      
      m_geo_volume_mode = false;
      m_volume_rand_present = false;
      
      m_X0 = pset.get<double>("X0");
      m_Y0 = pset.get<double>("Y0");
      m_Z0 = pset.get<double>("Z0");
      m_X1 = pset.get<double>("X1");
      m_Y1 = pset.get<double>("Y1");
      m_Z1 = pset.get<double>("Z1");

      CalculateActiveVolumeFromXYZ();
      
      
    } else {
      const TGeoNode* world = gGeoManager->GetNode(0);

      if (pset.get_if_present<std::string>("volume_rand", m_volume_rand)) {
        m_geo_volume_mode = false;
        m_volume_rand_present = true;

        /// If we specified volume_rand we will throw with what the coordinates of the volume
        const TGeoNode* node = nullptr;
        findNode(world, m_volume_rand, node);
        GetNodeXYZMinMax(node, world,
                         m_X0, m_X1,
                         m_Y0, m_Y1,
                         m_Z0, m_Z1);

        CalculateActiveVolumeFromXYZ();
        
      } else {
        /// Else we throw in the whole geometry but only in the nodes who have the specified volume/material
        m_geo_volume_mode = true;
        m_volume_rand_present = false;
        CalculateActiveVolumeFromAllNodes();
      }
    }
      
    m_flat_distrib_xpos = pset.get<double>("flat_distribution_x", 1);
    m_flat_distrib_ypos = pset.get<double>("flat_distribution_y", 1);
    m_flat_distrib_zpos = pset.get<double>("flat_distribution_z", 1);

    m_distrib_xpos = nullptr;
    m_distrib_ypos = nullptr;
    m_distrib_zpos = nullptr;

    auto rand = CLHEP::RandFlat(m_engine);
    gRandom->SetSeed(rand.fireInt(std::numeric_limits<long>::max()));
    
    if (not m_flat_distrib_xpos) m_distrib_xpos = new TF1("distrib_x", pset.get<std::string>("distrib_x").c_str(), m_X0, m_X1);
    if (not m_flat_distrib_ypos) m_distrib_ypos = new TF1("distrib_y", pset.get<std::string>("distrib_y").c_str(), m_Y0, m_Y1);
    if (not m_flat_distrib_zpos) m_distrib_zpos = new TF1("distrib_z", pset.get<std::string>("distrib_z").c_str(), m_Z0, m_Z1);

  }
  
  void BaseRadioGen::CalculateActiveVolumeFromAllNodes() {
    FillAllNodes(gGeoManager->GetTopNode());
    m_volume_cc = 0;
      
    for (auto const& node: m_all_nodes) {
      auto volume   = node->GetVolume();
      auto material = node->GetMedium()->GetMaterial();
        
      auto good_vol = std::find(m_good_volumes  .begin(), m_good_volumes  .end(), volume  );
      auto good_mat = std::find(m_good_materials.begin(), m_good_materials.end(), material);
        
      bool good = good_mat != m_good_materials.end() and good_vol != m_good_volumes.end();
        
      if (good) {
        double capa = volume->GetShape()->Capacity();
        m_volume_cc += capa;
        m_good_nodes[node] = m_volume_cc; 
      }
        
    }
    MF_LOG_INFO("BaseRadioGen") << m_good_nodes.size() << " nodes (i.e. instance of the volumes) satisfy both the regexes.\n";
      
    if (m_good_nodes.size()==0) 
      throw cet::exception("BaseRadioGen") << "Didn't find an instance of material " << m_material << " and the volume " << m_volume_gen << " in the geometry.\n";
  }
  
  void BaseRadioGen::CalculateActiveVolumeFromXYZ() {
    
    m_volume_cc = (m_X1 - m_X0)*(m_Y1 - m_Y0)*(m_Z1 - m_Z0);
    
    if (m_material != ".*" || m_volume_gen != ".*") {
      MF_LOG_INFO("BaseRadioGen") << "Calculating the proportion of " << m_material << " and the volume " << m_volume_gen << " in the specified volume " << m_volume_rand << ".\n";
      size_t nfound=0;
      size_t ntries=0;
      
      double xyz[3];
      TGeoNode* node = nullptr;

      while (nfound<m_target_n_point_rate_calculation and ntries<m_max_tries_rate_calculation) {
        ntries++;
        node = nullptr;
        xyz[0] = m_X0 + (m_X1 - m_X0) * m_random_flat->fire(0, 1.);
        xyz[1] = m_Y0 + (m_Y1 - m_Y0) * m_random_flat->fire(0, 1.);
        xyz[2] = m_Z0 + (m_Z1 - m_Z0) * m_random_flat->fire(0, 1.);
        m_geo_manager->SetCurrentPoint(xyz);
        node = m_geo_manager->FindNode();
        if (!node) continue;
        if (node->IsOverlapping()) continue;

        auto good_mat = std::find(m_good_materials.begin(), m_good_materials.end(), node->GetMedium()->GetMaterial());
        auto good_vol = std::find(m_good_volumes  .begin(), m_good_volumes  .end(), node->GetVolume());
        bool good = good_mat != m_good_materials.end() and good_vol != m_good_volumes.end();

        if (!good) continue;
        nfound++;
      }

      if (nfound==0) {
        throw cet::exception("BaseRadioGen") << "Didn't find the material " << m_material << " and the volume " << m_volume_gen << " in the specified volume " << m_volume_rand << ".\n"
                                             << "Position of the box:\n"
                                             << m_X0 << " " << m_X1 <<"\n"
                                             << m_Y0 << " " << m_Y1 <<"\n"
                                             << m_Z0 << " " << m_Z1 <<"\n";
      }

      double proportion = (double) nfound / ntries;
      double proportion_error = proportion*sqrt(1. / nfound + 1. / ntries);
      m_volume_cc *= proportion;

      MF_LOG_INFO("BaseRadioGen") << "There is " << proportion*100. << "% (+/- " << proportion_error*100. << "%) of " << m_material << " and " << m_volume_gen << " in the specified volume ("
                                  << 100.*proportion_error/proportion << "% relat. uncert.).\n"
                                  << "If the uncertainty is too big, crank up the parameters \"max_tries_rate_calculation\" (default=40,000,000) and/or \"target_n_point_rate_calculation\" (default=100,000) in your fhicl.\n\n\n";
    }
  }

  
  void BaseRadioGen::produce(art::Event& evt) {
    m_nevent++;
    produce_radio(evt);
  }

  void BaseRadioGen::beginRun(art::Run& run) {
    art::ServiceHandle<geo::Geometry const> geo;
    run.put(std::make_unique<sumdata::RunData>(m_geo_service->DetectorName()));
    beginRun_radio(run);
  }

  void BaseRadioGen::beginJob() {
    if (m_isotope.empty()) {
      throw cet::exception("BaseRadioGen") << "m_isotope is empty, you need to fill it yourself in the constructor of your module\n";
    }

    DeclareOutputHistos();
    beginJob_radio();
    m_nevent=0;
  }

  void BaseRadioGen::endJob() {
    if (m_nevent) {
      m_pdg_TH1D->Scale(1./m_nevent);

      for (auto histo : m_pos_xy_TH2D) {
        m_pos_xy_TH2D[histo.first]->Scale(1./m_nevent);
        m_pos_xz_TH2D[histo.first]->Scale(1./m_nevent);
        m_dir_x_TH1D [histo.first]->Scale(1./m_nevent);
        m_dir_y_TH1D [histo.first]->Scale(1./m_nevent);
        m_dir_z_TH1D [histo.first]->Scale(1./m_nevent);
        m_mom_TH1D   [histo.first]->Scale(1./m_nevent);
        m_ke_TH1D    [histo.first]->Scale(1./m_nevent);
        m_time_TH1D  [histo.first]->Scale(1./m_nevent);
      }
    }
    endJob_radio();
  }

  bool BaseRadioGen::GetGoodPositionTime(TLorentzVector& position) {

    /// Deal with the time first
    double time = m_T0 + m_random_flat->fire()*(m_T1 - m_T0);

    
    if (not m_geo_volume_mode) {
      /// We use all the XYZ that have been set earlier, until we find the correct volume and material
      
      size_t n_tries=0;
    
      double xpos = std::numeric_limits<double>::signaling_NaN();
      double ypos = std::numeric_limits<double>::signaling_NaN();
      double zpos = std::numeric_limits<double>::signaling_NaN();
      
      while (n_tries++<m_max_tries_event) {
      
        if (m_flat_distrib_xpos) xpos = m_X0 + m_random_flat->fire()*(m_X1 - m_X0);
        else                     xpos = m_distrib_xpos->GetRandom(m_X0, m_X1);
      
        if (m_flat_distrib_ypos) ypos = m_Y0 + m_random_flat->fire()*(m_Y1 - m_Y0);
        else                     ypos = m_distrib_ypos->GetRandom(m_Y0, m_Y1);
      
        if (m_flat_distrib_zpos) zpos = m_Z0 + m_random_flat->fire()*(m_Z1 - m_Z0);
        else                     zpos = m_distrib_zpos->GetRandom(m_Z0, m_Z1);

        auto node = gGeoManager->FindNode(xpos, ypos, zpos);
        auto volume   = node->GetVolume();
        auto material = node->GetMedium()->GetMaterial();
        
        auto good_vol = std::find(m_good_volumes  .begin(), m_good_volumes  .end(), volume  );
        auto good_mat = std::find(m_good_materials.begin(), m_good_materials.end(), material);
        
        bool good = good_mat != m_good_materials.end() and good_vol != m_good_volumes.end();
        
        if (good) break;
      }
      
      if (std::isnan(xpos) or std::isnan(ypos) or std::isnan(zpos)) {
        MF_LOG_ERROR("BaseRadioGen") << "Error in generation of random position!";
      }
      
      position.SetXYZT(xpos, ypos, zpos, time);
      return true;
      
    } else {
      /// We use the list of nodes

      if (m_good_nodes.empty() or m_volume_cc == 0)
        MF_LOG_ERROR("BaseRadioGen") << "There is no node to throw events in!";
      
      double which_vol = m_random_flat->fire()*m_volume_cc;
    
      const TGeoNode* node = nullptr;
      size_t i=0;
      
      for (auto const& nd: m_good_nodes) {
        if (which_vol < nd.second) {
          node = nd.first;
          break;
        }
        i++;
      }
    
      const TGeoNode* world = gGeoManager->GetNode(0);
      double x_min, x_max;
      double y_min, y_max;
      double z_min, z_max;
      GetNodeXYZMinMax(node, world,
                       x_min, x_max,
                       y_min, y_max,
                       z_min, z_max);
      
      size_t n_tries=0;
      

      while (n_tries++<m_max_tries_event) {
        double xpos = std::numeric_limits<double>::signaling_NaN();
        double ypos = std::numeric_limits<double>::signaling_NaN();
        double zpos = std::numeric_limits<double>::signaling_NaN();
        
        if (m_flat_distrib_xpos) xpos = x_min + m_random_flat->fire()*(x_max - x_min);
        else                     xpos = m_distrib_xpos->GetRandom(x_min, x_max);
      
        if (m_flat_distrib_ypos) ypos = y_min + m_random_flat->fire()*(y_max - y_min);
        else                     ypos = m_distrib_ypos->GetRandom(y_min, y_max);
      
        if (m_flat_distrib_zpos) zpos = z_min + m_random_flat->fire()*(z_max - z_min);
        else                     zpos = m_distrib_zpos->GetRandom(z_min, z_max);

        if (std::isnan(xpos) or std::isnan(ypos) or std::isnan(zpos)) {
          MF_LOG_ERROR("BaseRadioGen") << "Error in generation of random position!";
        }
        position.SetXYZT(xpos, ypos, zpos, time);
      
        const TGeoNode* node_generated = gGeoManager->FindNode(xpos, ypos, zpos);
        if (node_generated == node)
          return true;
      
        if (node->GetNdaughters() and IsDaughterNode(node_generated, node))
          return true;
      }
    }
    
    return false;
    
  }

  //____________________________________________________________________________
  // Generate radioactive decays per isotope per volume according to the FCL parameters
  int BaseRadioGen::GetNDecays() {

    if (m_rate_mode) {
      return m_random_poisson->fire(m_rate);
    } else {
      double rate = abs(m_Bq * (m_T1-m_T0) * m_volume_cc / 1.0E9);
      int n_ev = m_random_poisson->fire(rate);
      return n_ev;
    }

  }

  bool BaseRadioGen::IsDaughterNode(const TGeoNode* daughter_node, const TGeoNode* mother_node) {
    if (mother_node == daughter_node)
      return true;

    TObjArray* daunodes = mother_node->GetNodes();
    if (!daunodes)
      return false;
    
    TIter next(daunodes);
    const TGeoNode* anode = 0;
    
    while ( (anode = (const TGeoNode*)next()) ) {
      bool found = IsDaughterNode(daughter_node,anode);
      if (found)
        return true;
    }

    return false;
    
  }
  
  void BaseRadioGen::FillAllNodes(const TGeoNode* curnode) {
    if (!curnode) return;
    TObjArray* daunodes = curnode->GetNodes();
    if (!daunodes) return;
  
    TIter next(daunodes);

    const TGeoNode* anode = 0;
    while ( (anode = (const TGeoNode*)next()) ) {
      m_all_nodes.insert(anode);
      FillAllNodes(anode);
    }
  }
  
  bool BaseRadioGen::findNode(const TGeoNode* curnode, std::string& tgtnname,
                              const TGeoNode* & targetnode)
  /// Shamelessly stolen from here: https://cdcvs.fnal.gov/redmine/attachments/6719/calc_bbox.C
  {
    std::string nname = curnode->GetName();
    std::string vname = curnode->GetVolume()->GetName();
    if ( nname == tgtnname || vname == tgtnname ) {
      targetnode = curnode;
      return true;
    }

    TObjArray* daunodes = curnode->GetNodes();
    if ( ! daunodes ) return false;
    TIter next(daunodes);
    const TGeoNode* anode = 0;
    while ( (anode = (const TGeoNode*)next()) ) {
      bool found = findNode(anode,tgtnname,targetnode);
      if ( found ) return true;
    }

    return false;
  }


  bool BaseRadioGen::findMotherNode(const TGeoNode* cur_node, std::string& daughter_name,
                                    const TGeoNode* & mother_node)
  /// Adapted from above
  {
    TObjArray* daunodes = cur_node->GetNodes();

    if ( ! daunodes ) return false;

    TIter next(daunodes);

    const TGeoNode* anode = 0;
    bool found = 0;
    while ( (anode = (const TGeoNode*)next()) ) {
      std::string nname = anode->GetName();
      std::string vname = anode->GetVolume()->GetName();

      if ( nname == daughter_name || vname == daughter_name ) {
        mother_node = cur_node;
        return true;
      }
      found = findMotherNode(anode, daughter_name, mother_node);

    }

    return found;
  }

  void BaseRadioGen::DeclareOutputHistos() {
    art::ServiceHandle<art::TFileService> tfs;
    auto pdgs = {1000020040,11,-11,22,2112,2212,9999};
    m_pdg_TH1D = tfs->make<TH1D>("PDG", ";PDG;n particles/event", pdgs.size(), 0, pdgs.size());

    for (auto pdg : pdgs) {

      std::string part="";
      int simple_pdg;
      SimplePDG(pdg, simple_pdg, part);
      art::TFileDirectory dir = tfs->mkdir(part.c_str());

      m_pos_xy_TH2D[simple_pdg] = dir.make<TH2D>("posXY"   , ";X [cm];Y [cm]"             , 100, m_X0, m_X1, 100, m_Y0, m_Y1);
      m_pos_xz_TH2D[simple_pdg] = dir.make<TH2D>("posXZ"   , ";X [cm];Z [cm]"             , 100, m_X0, m_X1, 100, m_Z0, m_Z1);
      m_dir_x_TH1D [simple_pdg] = dir.make<TH1D>("dirX"    , ";X momentum projection"     , 100,   -1,   1);
      m_dir_y_TH1D [simple_pdg] = dir.make<TH1D>("dirY"    , ";Y momentum projection"     , 100,   -1,   1);
      m_dir_z_TH1D [simple_pdg] = dir.make<TH1D>("dirZ"    , ";Z momentum projection"     , 100,   -1,   1);
      m_mom_TH1D   [simple_pdg] = dir.make<TH1D>("Momentum", ";Momentum [MeV];n particles/event", 5000,   0, 500);
      m_ke_TH1D    [simple_pdg] = dir.make<TH1D>("KE"      , ";KE [MeV];n particles/event"      , 5000,   0, 500);
      m_time_TH1D  [simple_pdg] = dir.make<TH1D>("Time"    , ";Time[ns];n particles/event"      , 100, m_T0, m_T1);

      m_pdg_TH1D->GetXaxis()->SetBinLabel(simple_pdg+1, part.c_str());
    }
  }

  void BaseRadioGen::SimplePDG(int pdg, int& simple, std::string& name) {
    switch (pdg) {
    case 1000020040:
      name = "alpha";
      simple = 0;
      return;
    case 22:
      name = "gamma";
      simple = 1;
      return;
    case -11:
      name = "positron";
      simple = 2;
      return;
    case 11:
      name = "electron";
      simple = 3;
      return;
    case 2112:
      name = "neutron";
      simple = 4;
      return;
    case 2212:
      name = "proton";
      simple = 5;
      return;
    default:
      name = "other";
      simple = 6;
      return;
    }
  }


  void BaseRadioGen::FillHistos(simb::MCParticle& part) {
    int pdg = part.PdgCode();
    int simple_pdg = 0;
    std::string dummy="";
    SimplePDG(pdg, simple_pdg, dummy);

    TLorentzVector position = part.Position();
    TLorentzVector mom = part.Momentum();
    double mass = mom.M();
    m_pos_xy_TH2D[simple_pdg]->Fill(position.X(), position.Y());
    m_pos_xz_TH2D[simple_pdg]->Fill(position.X(), position.Z());
    m_dir_x_TH1D [simple_pdg]->Fill(mom.Px()/mom.P());
    m_dir_y_TH1D [simple_pdg]->Fill(mom.Py()/mom.P());
    m_dir_z_TH1D [simple_pdg]->Fill(mom.Pz()/mom.P());
    double ke = (sqrt(mom.P()*mom.P()+mass*mass)-mass)*1000.;
    m_ke_TH1D    [simple_pdg]->Fill(ke);
    m_mom_TH1D   [simple_pdg]->Fill(mom.P()*1000.);
    m_time_TH1D  [simple_pdg]->Fill(position.T());
    m_pdg_TH1D->Fill(simple_pdg);
  }

  TLorentzVector BaseRadioGen::dirCalc(double p, double m) {
    // isotropic production angle for the decay product
    double costheta = (2.0*m_random_flat->fire() - 1.0);

    if (costheta < -1.0) costheta = -1.0;
    if (costheta > 1.0) costheta = 1.0;

    double const sintheta = sqrt(1.0-costheta*costheta);
    double const phi = 2.0*M_PI*m_random_flat->fire();

    return TLorentzVector{p*sintheta*std::cos(phi),
        p*sintheta*std::sin(phi),
        p*costheta,
        std::sqrt(p*p+m*m)};
  }

}//end namespace evgen
