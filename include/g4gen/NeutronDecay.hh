/// Defines classes for generating neutron decays of various types
#ifndef G4GEN_NEUTRON_DECAY_HH
#define G4GEN_NEUTRON_DECAY_HH
#include <map>
#include <memory>
#include <G4String.hh>
#include <G4ThreeVector.hh>
#include <G4LorentzVector.hh>

namespace g4gen {

class Rng;
class Particle;


/// Abstract base class for generic neutron decay generators
class NeutronDecay {
public:
  /// Ctor
	NeutronDecay() { }
	/// Dtor
	virtual ~NeutronDecay() { }
	/// Set input particle (in excited unbound state)
	/** Parameters for initial decaying nucleus are taken from
	 *  the input of this function.
	 *  \attention It shold be called every time, before calling Generate()
	 */
	virtual void SetInputParticle(const Particle* p) = 0;
	/// Create the RNG used for excitation energy generation
	/** Returned instance is wrapped in "smart" pointer class std::unique_ptr
	 */
	virtual std::unique_ptr<Rng> CreateRngEx() = 0;
	/// Set verbosity level
	/** 0: Print nothing
	 *  1: Print fatal errror/warning messages only (to G4cerr)
	 *  2: Print all error/warning messages (to G4cerr)
	 */
	virtual void SetVerboseLevel(G4int level) = 0;
	/// Returns vebosity level
	virtual G4int GetVerboseLevel() const = 0;
	/// Returns number of neutrons emitted in the decay
	virtual G4int GetNumberOfNeutrons() const = 0;
	/// Set some generic decay parameter (id'd by string)
	virtual void SetDecayParam(const G4String& par, G4double val) = 0;
	/// Return some generic decay paramter (id'd by string)
	virtual G4double GetDecayParam(const G4String& par) = 0;
	/// Do the actual neutron decay generation
	/** Needs to be implemented separately for every decay type.
	 *  \returns `true` if decay was successful, `false` if there
	 *            wasn't enough energy in the initial state 
	 *            for the decay to happen
	 */
	virtual G4bool Generate() = 0;
	/// Return the final states, after calling Generate()
	/** \param indx 0: initial state; 1: final fragment; 2,3,...: neutrons 0,1,...
	 */
	virtual const G4LorentzVector& GetFinal(G4int indx) const = 0;
};


/// 'Intermediate' NeutronDecay class, implementing many of the details
/// common to all decay types.
/** Still abstract, need to implement Generate() */
class NeutronDecayIntermediate : public NeutronDecay {
public:
	NeutronDecayIntermediate(G4int number_of_neutrons_emitted);
	virtual ~NeutronDecayIntermediate();
	virtual void SetInputParticle(const Particle* p);
	virtual std::unique_ptr<Rng> CreateRngEx() = 0;
	virtual G4int GetNumberOfNeutrons() const { return mNumberOfNeutronsEmitted; }
	virtual G4bool Generate() = 0;
	virtual void SetDecayParam(const G4String& par, G4double val);
	virtual G4double GetDecayParam(const G4String& par);
	virtual void SetVerboseLevel(G4int level) { fVerb = level; }
	virtual G4int GetVerboseLevel() const { return fVerb; }
protected:
	void SetFinal(G4int indx, const G4LorentzVector& v);
	virtual const G4LorentzVector& GetFinal(G4int indx) const;
private:
	G4int mNumberOfNeutronsEmitted;
	std::map<G4String, G4double> mParams;
	std::vector<G4LorentzVector> mFinal;
	G4int fVerb;
protected:
	G4double mFinalFragMass; // Rest mass of final decay fragment
	const Particle* mInitial;
};


/// Implements CreateRngEx() for symmetric breit-wigner
/// excitation energy distributions
/** Parameters to set are "ex" and "width"
 *  Setting "width" to zero returns a uniform (spike) decay energy
 */
class NeutronDecayBreitWigner : public NeutronDecayIntermediate {
public:
	NeutronDecayBreitWigner(G4int nneut):
		NeutronDecayIntermediate(nneut),
		mRngEx(0) { }
	virtual ~NeutronDecayBreitWigner()  { }
	virtual G4bool Generate() = 0;
	virtual std::unique_ptr<Rng> CreateRngEx();
protected:
	/// Pointer to created RNG (NO OWNERSHIP)
	/** \attention Derived classes have responsibility to ensure
	 *  its validity!
	 */
	const Rng* mRngEx;
};


/// Concrete class for single neutron decay, Breit Wigner
/** Parameters to set are "energy" and "width"
 *  Setting "width" to zero returns a uniform (spike) decay energy
 */
class OneNeutronDecay : public NeutronDecayBreitWigner {
public:
	OneNeutronDecay();
	virtual ~OneNeutronDecay();
	virtual G4bool Generate();
};

/// Concrete class for two neutron  phase space decay. Includes optional
/// final state interaction (FSI).
/**  If used, FSI is taken from the code by
 *   F. Marquis, sent originally by him to J.K. Smith, then to GAC. 
 *   Reference for the FSI calculation is PLB 476, 219 (2000), 
 *   https://doi.org/10.1016/S0370-2693(00)00141-6
 */
class TwoNeutronDecayPhaseSpace : public NeutronDecayBreitWigner {
public:
	/// Ctor
	/** \param [in] fsi If true, include final state interaction in calculation;
	 *   if false, do not include FSI.
	 */
	TwoNeutronDecayPhaseSpace(G4bool fsi = FALSE);
	virtual ~TwoNeutronDecayPhaseSpace();
	virtual G4bool Generate();
private:
	G4bool fFSI;
};


#if 0
/// Concrete class for two neutron 'dineutron' decay
/** Dineutron decay is calculated using the formalism developed
 *  by A. Volya, PRC 76, 064314, 2006 && EPJ Web Conf., 38, 03003, 2012
 */
class TwoNeutronDecayDiNeutron : public NeutronDecayIntermediate {
public:
	TwoNeutronDecayDiNeutron();
	virtual ~TwoNeutronDecayDiNeutron();
	virtual std::unique_ptr<Rng> CreateRngEx();
	virtual G4bool Generate();
private:
	const RngVolyaDiNeutronEx* mRngEx;
};

/// Concrete class for two neutron sequential decay
/** Calculated using the formalism developed
 *  by A. Volya, PRC 76, 064314, 2006 && EPJ Web Conf., 38, 03003, 2012
 */
class TwoNeutronDecaySequential : public NeutronDecayIntermediate {
public:
	TwoNeutronDecaySequential();
	virtual ~TwoNeutronDecaySequential();
	virtual void SetInputParticle(const Particle* p);
	virtual std::unique_ptr<Rng> CreateRngEx();
	virtual G4bool Generate();
public:
	G4double mIntermediateFragMass;
	const RngVolyaSequentialEx* mRngEx;
};
#endif


/// Factory class
/** Create different types of NeutronDecay instances,
 *  based on inputs
 */
class NeutronDecayFactory {
public:
	/// Set the type of decay, determines what child class is created
	void SetDecayType(G4String type) { fDecayType = type; }
	/// Returns the decay type
	const G4String& GetDecayType() const { return fDecayType; }
	/// Set an optional parameter to be passed to the created NeutronDecay instance
	void SetDecayOption(G4String option, G4double value) { mOptions[option] = value; }
	/// Receive some optional decay parameter
	G4double GetDecayOption(G4String option) const;
	/// Create new instance of NeutronDecay
	/** User is responsible for deleting the returned value */
	NeutronDecay* Create();
private:
	G4String fDecayType;
	std::map<G4String, G4double> mOptions;
};

/// Helper class to calculate neutron evaporation
class NeutronEvaporation {
public:
	/// Ctor
	/** \param m0 Mass of initial state (G.S. + Excitation)
	 *  \param mf Mass of final "fragment" state (G.S. + Excitation)
	 *  \param mn "Neutron" mass (can be different, if 2nd particle isn't a neutron)
	 */
	NeutronEvaporation(G4double m0, G4double mf, G4double mn):
		mM0(m0), mMf(mf), mMn(mn) { }
	/// Dtor
	~NeutronEvaporation() {}
	/// Calculates the neutron evaporation in the CENTER OF MASS frame
	/** \param [out] Frag : Final "fragment" 4-vector in COM
	 *  \param [out] Neut : Final "neutron" 4-vector in COM
	 */
	void operator() (G4LorentzVector* Frag, G4LorentzVector* Neut);
private:
	G4double mM0, mMf, mMn;
};



}


#endif
