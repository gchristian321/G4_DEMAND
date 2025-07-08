/// \brief Defines random number generator classes
#ifndef G4GEN_RNG_HH
#define G4GEN_RNG_HH
#include <memory>
#include <vector>
#include <utility>
#include <globals.hh>



namespace g4gen {

/// Set Global RNG seed
extern void SetRngSeed(G4int seed);

/// Get Global RNG seed
extern G4int GetRngSeed();

/// Abstract random number generator class
/** Derived classes implement details of generating random
 *  numbers from set distributions.
 *  \attention All derived classes should use the CLHEP random number
 *   generators ONLY. This ensures consistency throughout the simulation,
 *   in particular it means only one seed value needs to be changed to
 *   ensire an independent simulation.
 */
class Rng {
public:
	/// Ctor
	Rng() { }
	/// Dtor
	virtual ~Rng() { }
	/// Generate a random number, and return the value
	G4double Generate() { return fR = DoGenerate(); }
	/// Generate a random number above some value, and return it
	G4double GenerateAbove(G4double low);
	/// Generate a random number below some value, and return it
	G4double GenerateBelow(G4double high);
	/// Generate a random number between two values, and return it
	G4double GenerateBetween(G4double low, G4double high);
	/// Get the most recently generated random number
	G4double GetLast() const { return fR; }
private:
	/// Implements the actual random number generation
	/** \returns The generated random number
	 */
	virtual G4double DoGenerate() = 0;
	/// Most recently generated random number
	G4double fR;
};

/// Abstract random number generator class, generating two
/// correlated random numbers from a distribution
/** Derived classes implement details of generating random
 *  numbers from set distributions.
 *  \attention All derived classes should use the CLHEP random number
 *   generators ONLY. This ensures consistency throughout the simulation,
 *   in particular it means only one seed value needs to be changed to
 *   ensire an independent simulation.
 */
class Rng2d {
public:
	/// Ctor
	Rng2d() { }
	/// Dtor
	virtual ~Rng2d() { }
	/// Generate two correlated random numbers and return then
	/// as a `std::pair`
	std::pair<G4double, G4double> Generate() { return fR2 = DoGenerate(); }
	/// Get the most recently generated pair of random numbers
	std::pair<G4double, G4double> GetLast() const { return fR2; }
private:
	/// Implements the actual random number generation
	/** \returns The generated random number pair
	 */
	virtual std::pair<G4double, G4double> DoGenerate() = 0;
	/// Most recently generated pair
	std::pair<G4double, G4double> fR2;
};


///  Always returns a constant value
/** Not really a random number generator, but made to look like one
 */
class RngConstant : public Rng {
public:
	///  Ctor
	/** \param [in] val The constant number to return */
	RngConstant(G4double val): fVal(val) { }
private:
	G4double DoGenerate() { return fVal; }
	G4double fVal;
};	

/// Gaussian [normal] distribution
class RngGaus : public Rng {
public:
	/// Ctor
	/** \param [in] mean Mean value of gaussian
	 *  \param [in] sigma Standard deviation of gaussian
	 */
	RngGaus(G4double mean, G4double sigma);
	~RngGaus();
private:
	G4double DoGenerate();
private:
	G4double mMean, mSigma;
};

/// Breit Wigner distribution, as in resonance decay
class RngBreitWigner : public Rng {
public:
	/// Ctor
	/** \param [in] mean Mean value
	 *  \param [in] width Width
	 */
	RngBreitWigner(G4double mean, G4double width);
	~RngBreitWigner();
private:
	G4double DoGenerate();
private:
	G4double mMean, mSigma;
};

/// Uniform distribution
class RngUniform : public Rng {
public:
	/// Ctor
	/** \param [in] low Low value of uniform range (inclusive)
	 *  \param [in] high High value of uniform range (exclusive)
	 */
	RngUniform(G4double low = 0, G4double high = 1);
	~RngUniform();
private:
	G4double DoGenerate();
private:
	G4double mLow, mHigh;
};

/// Uses custom distrubition from a text file
class RngCustom : public Rng {
public:
	/// Ctor
	/** \param [in] filename Name of the file containing the custom
	 *  distribution.
	 *
	 *  Should be structured like a histogram with the format:
	 *  \code
	 *  bin_low <whitespace> bin_value
	 *  \endcode
	 */
	RngCustom(const G4String& filename);
	~RngCustom();
private:
	G4double DoGenerate();
protected:
	std::vector<G4double> mXlow, mCdf;
};

/// Custom angular distribution for nuclear reactions
class RngCustomAngDist : public RngCustom {
public:
	/// Ctor
	/** \param [in] filename Name of the file containing the custom
	 *  distribution.
	 *  
	 *  Should be structured like a histogram with the format:
	 *  \code
	 *  bin_low <whitespace> bin_value
	 *  \endcode
	 *  The bin values should be dSigma/dOmega. Multiplication by
	 *  sin(theta) is taken care of automatically.
	 */
	RngCustomAngDist(const G4String& filename);
	~RngCustomAngDist();
};

/// Two-dimensional Gaussian [normal] with correlated x,y
/** No option to shift the means, so just add an offset value
 *  manually.
 */
class RngGaus2d : public Rng2d {
public:
	/// Ctor
	/** \param [in] sigma_x standard deviation in x
	 *  \param [in] sigma_y standard deviation in y
	 *  \param [in] rho  Correlation coefficient
	 */
	RngGaus2d(G4double sigma_x, G4double sigma_y, G4double rho);
private:
	virtual std::pair<G4double, G4double> DoGenerate();
private:
	G4double mSigmaX, mSigmaY, mRho;
};

/// Class to prevent infinite loops when trying to
/// generate RNGs within a range
class CheckMaxTries {
public:
	/// Ctor
	/** \param [in] max Max number of tries before complaining.
	 */
	CheckMaxTries(G4int max=10000): kMaxTries(max) { }
	/// Check if a loop has reached the maximum number of tries.
	/** If it has, give user the option to exit the program, or to continue
	 */
	void operator() (G4int& n, const char* fct, G4double* low = 0, G4double* high = 0);
private:
	const G4int kMaxTries;
};

}


#endif
