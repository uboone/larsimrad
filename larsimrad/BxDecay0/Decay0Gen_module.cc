/// \file  Decay0Gen_module.cc
/// \brief Generator for radiological decays
/// Module designed to produce a set list of particles for MC to model radiological decays
/// \author  plasorak@FNAL.GOV
///          April 2020 PLasorak

#include "larsimrad/BxDecay0/BaseRadioGen.h"

#include <bxdecay0/i_random.h>
#include <bxdecay0/event.h>            // Decay event data model
#if defined __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wmismatched-tags"
#include <bxdecay0/decay0_generator.h> // Decay0 generator with OOP interface
#pragma clang diagnostic pop
#else 
#include <bxdecay0/decay0_generator.h> // Decay0 generator with OOP interface
#endif
#include <bxdecay0/particle.h>



namespace evgen {
  /// Module to generate particles created by radiological decay, patterend off of SingleGen
  /// Currently it generates only in rectangular prisms oriented along the x,y,z axes

  class Decay0Gen : public evgen::BaseRadioGen {
  public:
    explicit Decay0Gen(fhicl::ParameterSet const& pset);
    ~Decay0Gen() {
      for (auto& i: m_decay0_generator) {
        i->reset();
      }
    }

  private:
    // This is called for each event.
    void produce_radio(art::Event& evt);
    void beginJob_radio();
    void endJob_radio();

    std::shared_ptr<clhep_random> m_random_decay0;

    // Declare a Decay0 generator:
    std::vector<std::unique_ptr<bxdecay0::decay0_generator>> m_decay0_generator;
    TH1D* m_timediff_TH1D;
    bool m_single_isotope_mode;
  };


  /// \brief Wrapper functor for a standard random number generator
  class clhep_random : public bxdecay0::i_random{
  public:
    /// Constructor
    clhep_random(CLHEP::HepRandomEngine& gen):
      m_generator(gen),
      m_rand_flat(gen) { }

    /// Main operator
    virtual double operator()(){
      double v = m_rand_flat.fire(0.,1.);
      return v;
    }
    CLHEP::HepRandomEngine& m_generator;
    CLHEP::RandFlat m_rand_flat;
    virtual ~clhep_random() {};
  };

}

namespace evgen{

  Decay0Gen::Decay0Gen(fhicl::ParameterSet const& pset):
    BaseRadioGen(pset) {
    std::string isotope="";
    m_single_isotope_mode = pset.get_if_present<std::string>("isotope", isotope);


    if (not m_single_isotope_mode) {
      fhicl::ParameterSet decay_chain = pset.get<fhicl::ParameterSet>("decay_chain");
      int index=0;

      while (decay_chain.get_if_present<std::string>("isotope_"+std::to_string(index++), isotope)) {
        m_isotope.push_back(isotope);
      }

    } else {
      m_isotope.push_back(isotope);
    }

    m_random_decay0 = std::make_shared<clhep_random>(GetRandomEngine());

    for (auto const& isotope: m_isotope) {
      auto generator = std::make_unique<bxdecay0::decay0_generator>();
      generator->reset();

      // Configure the Decay0 generator:
      generator->set_decay_category(bxdecay0::decay0_generator::DECAY_CATEGORY_BACKGROUND);
      generator->set_decay_isotope(isotope.c_str());
      try{
        generator->initialize(*m_random_decay0);
      } catch (...) {
        throw cet::exception("Decay0Gen") << "The inialisation of Decay0 failed. Maybe the isotope " << isotope << " doesn't exists?\n";
      }
      m_decay0_generator.push_back(std::move(generator));
    }

  }

  void Decay0Gen::beginJob_radio() {
    art::ServiceHandle<art::TFileService> tfs;
    double T0=0, T1=0;
    GetTs(T0, T1);
    m_timediff_TH1D = tfs->make<TH1D>("TimeDiff", ";Time Diff[ns];n particles" , (int)((T1)/10000), 0, T1);
  }

  void Decay0Gen::endJob_radio() {
    if (GetNEvents())
      m_timediff_TH1D->Scale(1./GetNEvents());
  }

  //____________________________________________________________________________
  void Decay0Gen::produce_radio(art::Event& evt) {
    //unique_ptr allows ownership to be transferred to the art::Event after the put statement
    std::unique_ptr< std::vector<simb::MCTruth> > truthcol(new std::vector<simb::MCTruth>);

    simb::MCTruth truth;
    truth.SetOrigin(simb::kSingleParticle);
    int track_id=-1;
    const std::string primary_str("primary");
    int bad_position=0;
    int tentative_decay=0;
    for (auto const& decay0_gen: m_decay0_generator) {
      int n_decay = GetNDecays();
      
      tentative_decay+=n_decay;
      for (int iDecay=0; iDecay<n_decay; ++iDecay) {
        TLorentzVector position;
        if (GetGoodPositionTime(position)) {

          bxdecay0::event gendecay;     // Declare an empty decay event
          decay0_gen->shoot(*m_random_decay0, gendecay); // Randomize the decay event
          std::vector<bxdecay0::particle> part = gendecay.get_particles();

          double t_alpha=0;
          double t_electron=0;
          for (auto const& p: part) {
            // alpha particles need a little help since they're not in the TDatabasePDG table
            // so don't rely so heavily on default arguments to the MCParticle constructor
            simb::MCParticle part;
            if (not p.is_valid()) {
              p.print(std::cout);
              throw cet::exception("Decay0Gen") << "Invalid part generated by Decay0, printed above (no clue what that means so throw)";
            }

            double mass = bxdecay0::particle_mass_MeV(p.get_code());
            int pdg=0;
            //int simple_pdg=0; // simple_pdg is unused here
            
            if      (p.is_alpha   ()) { pdg = 1000020040; part = simb::MCParticle(track_id, pdg, primary_str,-1,mass/1000,1); }
            else if (p.is_gamma   ()) { pdg =         22; part = simb::MCParticle(track_id, pdg, primary_str); }
            else if (p.is_positron()) { pdg =        -11; part = simb::MCParticle(track_id, pdg, primary_str); }
            else if (p.is_electron()) { pdg =         11; part = simb::MCParticle(track_id, pdg, primary_str); }
            else if (p.is_neutron ()) { pdg =       2112; part = simb::MCParticle(track_id, pdg, primary_str); }
            else if (p.is_proton  ()) { pdg =       2212; part = simb::MCParticle(track_id, pdg, primary_str); }
            else {
              p.print(std::cout);
              throw cet::exception("Decay0Gen") << "Particle above is weird, cannot recognise it.";
            }
            
            if ((p.is_positron() or p.is_electron()) and t_electron == 0) {
              t_electron = p.get_time();
            }
            if (p.is_alpha() and t_alpha == 0) {
              t_alpha = p.get_time();
            }
            track_id--;
            TLorentzVector mom(p.get_px()/1000.,
                               p.get_py()/1000.,
                               p.get_pz()/1000.,
                               sqrt(p.get_p()*p.get_p() + mass*mass)/1000.);
            
            TLorentzVector this_part_position = position;
            double t = position.T();
            this_part_position.SetT(t+p.get_time()*1e9);

            part.AddTrajectoryPoint(this_part_position, mom);
            
            truth.Add(part);
            
            FillHistos(part);
          }
          if (abs(t_alpha-t_electron)>1.e-15){
            m_timediff_TH1D->Fill(1e9*abs(t_alpha-t_electron));
          }
          gendecay.reset();
        } else { // !GetGoodPosition
          ++bad_position;
        }
      } // idecay
    } // m_decay_generators

    MF_LOG_DEBUG("Decay0Gen") << truth;
    if (bad_position>0) {
      MF_LOG_ERROR("Decay0Gen") << "There were " << bad_position << " failed attempts to get a good position to generate a decay out of the target " << tentative_decay << " decays.\n"
                                << "If these 2 numbers are close together, it means that the rate of your background is wrong (underestimated).\n"
                                << "You can fix this by increasing the parameter \"max_tries_event\" to a bigger number (default is 1M) in your fhicl.\n"
                                << "Another way is to change the \"volume_rand\" to a smaller one.\n";
      
    }
    truthcol->push_back(truth);
    evt.put(std::move(truthcol));
  }

}//end namespace evgen

DEFINE_ART_MODULE(evgen::Decay0Gen)
