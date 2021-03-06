#ifndef SLIPHARDEN_H
#define SLIPHARDEN_H

#include "sliprules.h"
#include "crystallography.h"

#include "../objects.h"
#include "../history.h"
#include "../interpolate.h"

#include "../math/tensors.h"
#include "../math/rotations.h"

#include "../windows.h"

#include <map>
#include <vector>
#include <string>

namespace neml {

class SlipRule; // Why would we need a forward declaration?

/// ABC for a slip hardening model
class NEML_EXPORT SlipHardening: public NEMLObject
{
 public:
  /// Request whatever history you will need
  virtual void populate_history(History & history) const = 0;
  /// Setup history
  virtual void init_history(History & history) const = 0;

  /// Map the set of history variables to the slip system hardening
  virtual double hist_to_tau(size_t g, size_t i, const History & history,
                             double T, const History & fixed) const = 0;
  /// Derivative of the map wrt to history
  virtual History
      d_hist_to_tau(size_t g, size_t i, const History & history,
                    double T, const History & fixed) const = 0;

  /// The rate of the history
  virtual History hist(const Symmetric & stress,
                     const Orientation & Q, const History & history,
                     Lattice & L, double T, const SlipRule & R,
                     const History & fixed) const = 0;
  /// Derivative of the history wrt stress
  virtual History d_hist_d_s(const Symmetric & stress,
                             const Orientation & Q, const History & history,
                             Lattice & L, double T,
                             const SlipRule & R,
                             const History & fixed) const = 0;
  /// Derivative of the history wrt the history
  virtual History
      d_hist_d_h(const Symmetric & stress,
                 const Orientation & Q,
                 const History & history,
                 Lattice & L,
                 double T, const SlipRule & R,
                 const History & fixed) const = 0;

  /// Whether this particular model uses the Nye tensor
  virtual bool use_nye() const;
};

/// Slip strength rules where all systems share the same strength
class NEML_EXPORT SlipSingleHardening: public SlipHardening
{
 public:
  /// Map the set of history variables to the slip system hardening
  virtual double hist_to_tau(size_t g, size_t i, const History & history,
                             double T, const History & fixed) const;
  /// Derivative of the map wrt to history
  virtual History
      d_hist_to_tau(size_t g, size_t i, const History & history,
                    double T, const History & fixed) const;

  /// The scalar map
  virtual double hist_map(const History & history, double T, 
                          const History & fixed) const = 0;
  /// The derivative of the scalar map
  virtual History d_hist_map(const History & history, double T,
                             const History & fixed) const = 0;
};

/// Slip strength rule where all systems evolve on a single scalar strength
class NEML_EXPORT SlipSingleStrengthHardening: public SlipSingleHardening
{
 public:
  SlipSingleStrengthHardening(std::string var_name = "strength");

  /// Request whatever history you will need
  virtual void populate_history(History & history) const;
  /// Setup history
  virtual void init_history(History & history) const;

  /// The rate of the history
  virtual History hist(const Symmetric & stress,
                     const Orientation & Q, const History & history,
                     Lattice & L, double T, const SlipRule & R,
                     const History & fixed) const;
  /// Derivative of the history wrt stress
  virtual History d_hist_d_s(const Symmetric & stress,
                             const Orientation & Q, const History & history,
                             Lattice & L, double T,
                             const SlipRule & R,
                             const History & fixed) const;
  /// Derivative of the history wrt the history
  virtual History
      d_hist_d_h(const Symmetric & stress,
                 const Orientation & Q,
                 const History & history,
                 Lattice & L,
                 double T, const SlipRule & R,
                 const History & fixed) const;

  /// The scalar map
  virtual double hist_map(const History & history, double T, 
                          const History & fixed) const;
  /// The derivative of the scalar map
  virtual History d_hist_map(const History & history, double T,
                             const History & fixed) const;

  /// Set the variable's name
  void set_variable(std::string name);

  /// Static (not evolving) strength
  virtual double static_strength(double T) const = 0;

  /// Nye contribution (defaults to zero)
  double nye_contribution(const History & fixed, double T) const;

  /// Actual implementation of any Nye contribution (defaults to zero)
  virtual double nye_part(const RankTwo & nye, double T) const;

  /// Setup the scalar
  virtual double init_strength() const = 0;

  /// Scalar evolution law
  virtual double hist_rate(const Symmetric & stress, const Orientation & Q,
                           const History & history,
                           Lattice & L, double T, const SlipRule & R, 
                           const History & fixed) const = 0;
  /// Derivative of scalar law wrt stress
  virtual Symmetric d_hist_rate_d_stress(const Symmetric & stress, const Orientation & Q,
                                         const History & history,
                                         Lattice & L, double T,
                                         const SlipRule & R, 
                                         const History & fixed) const = 0;
  /// Derivative of scalar law wrt the scalar
  virtual History d_hist_rate_d_hist(const Symmetric & stress, const Orientation & Q,
                                     const History & history,
                                     Lattice & L, double T,
                                     const SlipRule & R, 
                                     const History & fixed) const = 0;

 protected:
  std::string var_name_;
};

/// Sum of individual SlipSingleStrenghHardening models (static strengths also
/// summed)
class NEML_EXPORT SumSlipSingleStrengthHardening: public SlipSingleHardening
{
 public:
  /// Initialize with a list of models
  SumSlipSingleStrengthHardening(std::vector<std::shared_ptr<SlipSingleStrengthHardening>>
                                 models);

  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

  /// Request whatever history you will need
  virtual void populate_history(History & history) const;
  /// Setup history
  virtual void init_history(History & history) const;

  /// The rate of the history
  virtual History hist(const Symmetric & stress,
                     const Orientation & Q, const History & history,
                     Lattice & L, double T, const SlipRule & R,
                     const History & fixed) const;
  /// Derivative of the history wrt stress
  virtual History d_hist_d_s(const Symmetric & stress,
                             const Orientation & Q, const History & history,
                             Lattice & L, double T,
                             const SlipRule & R, 
                             const History & fixed) const;
  /// Derivative of the history wrt the history
  virtual History
      d_hist_d_h(const Symmetric & stress,
                 const Orientation & Q,
                 const History & history,
                 Lattice & L,
                 double T, const SlipRule & R,
                 const History & fixed) const;

  /// The scalar map
  virtual double hist_map(const History & history, double T,
                          const History & fixed) const;
  /// The derivative of the scalar map
  virtual History d_hist_map(const History & history, double T,
                             const History & fixed) const;

  /// Whether this model uses the Nye tensor
  virtual bool use_nye() const;

 private:
  size_t nmodels() const;
  const std::vector<std::shared_ptr<SlipSingleStrengthHardening>> models_;
};

static Register<SumSlipSingleStrengthHardening> regSumSlipSingleStrengthHardening;

/// Slip strength rule where the single strength evolves with sum|dg|
class NEML_EXPORT PlasticSlipHardening: public SlipSingleStrengthHardening
{
 public:
  PlasticSlipHardening(std::string var_name = "strength");

  /// Scalar evolution law
  virtual double hist_rate(const Symmetric & stress, const Orientation & Q,
                           const History & history, Lattice & L, double T,
                           const SlipRule & R, const History & fixed) const;
  /// Derivative of scalar law wrt stress
  virtual Symmetric d_hist_rate_d_stress(const Symmetric & stress, const Orientation & Q,
                                         const History & history, Lattice & L, double T,
                                         const SlipRule & R,
                                         const History & fixed) const;
  /// Derivative of scalar law wrt the scalar
  virtual History d_hist_rate_d_hist(const Symmetric & stress, const Orientation & Q,
                                     const History & history, Lattice & L, double T,
                                     const SlipRule & R, 
                                     const History & fixed) const;

  /// Prefactor
  virtual double hist_factor(double strength, Lattice & L, double T, 
                             const History & fixed) const = 0;
  /// Derivative of the prefactor
  virtual double d_hist_factor(double strength, Lattice & L, double T,
                               const History & fixed) const = 0;

};

/// Everyone's favorite Voce model
class NEML_EXPORT VoceSlipHardening: public PlasticSlipHardening
{
 public:
  /// Initialize with the saturated strength, the rate constant, and a constant
  /// strength
  VoceSlipHardening(std::shared_ptr<Interpolate> tau_sat,
                    std::shared_ptr<Interpolate> b,
                    std::shared_ptr<Interpolate> tau_0,
                    std::shared_ptr<Interpolate> k,
                    std::string var_name = "strength");

  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

  /// Setup the scalar
  virtual double init_strength() const;

  /// Static strength
  virtual double static_strength(double T) const;

  /// Prefactor
  virtual double hist_factor(double strength, Lattice & L, double T,
                             const History & fixed) const;
  /// Derivative of the prefactor
  virtual double d_hist_factor(double strength, Lattice & L, double T,
                               const History & fixed) const;

  /// Dynamically determine if we're going to use the Nye tensor
  virtual bool use_nye() const;

  /// Actual nye contribution
  virtual double nye_part(const RankTwo & nye, double T) const;

 private:
  std::shared_ptr<Interpolate> tau_sat_;
  std::shared_ptr<Interpolate> b_;
  std::shared_ptr<Interpolate> tau_0_;
  std::shared_ptr<Interpolate> k_;
};

static Register<VoceSlipHardening> regVoceSlipHardening;

/// The simplest of all linear hardening models
class NEML_EXPORT LinearSlipHardening: public PlasticSlipHardening
{
 public:
  /// Initialize with the saturated strength, the rate constant, and a constant
  /// strength
  LinearSlipHardening(std::shared_ptr<Interpolate> tau0,
                      std::shared_ptr<Interpolate> k1,
                      std::shared_ptr<Interpolate> k2,
                      std::string var_name = "strength");

  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

  /// Setup the scalar
  virtual double init_strength() const;

  /// Static strength
  virtual double static_strength(double T) const;

  /// Prefactor
  virtual double hist_factor(double strength, Lattice & L, double T,
                             const History & fixed) const;
  /// Derivative of the prefactor
  virtual double d_hist_factor(double strength, Lattice & L, double T,
                               const History & fixed) const;

  /// Dynamically determine if we're going to use the Nye tensor
  virtual bool use_nye() const;

  /// Actual nye contribution
  virtual double nye_part(const RankTwo & nye, double T) const;

 private:
  std::shared_ptr<Interpolate> tau0_;
  std::shared_ptr<Interpolate> k1_;
  std::shared_ptr<Interpolate> k2_;
};

static Register<LinearSlipHardening> regLinearSlipHardening;

}

#endif // SLIPHARDEN_H
