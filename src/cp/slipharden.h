#ifndef SLIPHARDEN_H
#define SLIPHARDEN_H

#include "crystallography.h"
#include "sliprules.h"

#include "../history.h"
#include "../interpolate.h"
#include "../objects.h"

#include "../math/rotations.h"
#include "../math/tensors.h"

#include <map>
#include <string>
#include <vector>

#include <stdexcept>

namespace neml {

class SlipRule; // Why would we need a forward declaration?

/// ABC for a slip hardening model
class SlipHardening : public NEMLObject {
public:
  /// Request whatever history you will need
  virtual void populate_history(History &history) const = 0;
  /// Setup history
  virtual void init_history(History &history) const = 0;

  /// Map the set of history variables to the slip system hardening
  virtual double hist_to_tau(size_t g, size_t i, const History &history,
                             double T) const = 0;
  /// Derivative of the map wrt to history
  virtual History d_hist_to_tau(size_t g, size_t i, const History &history,
                                double T) const = 0;

  /// The rate of the history
  virtual History hist(const Symmetric &stress, const Orientation &Q,
                       const History &history, Lattice &L, double T,
                       const SlipRule &R) const = 0;
  /// Derivative of the history wrt stress
  virtual History d_hist_d_s(const Symmetric &stress, const Orientation &Q,
                             const History &history, Lattice &L, double T,
                             const SlipRule &R) const = 0;
  /// Derivative of the history wrt the history
  virtual History d_hist_d_h(const Symmetric &stress, const Orientation &Q,
                             const History &history, Lattice &L, double T,
                             const SlipRule &R) const = 0;
};

/// Slip strength rules where all systems share the same strength
class SlipSingleHardening : public SlipHardening {
public:
  /// Map the set of history variables to the slip system hardening
  virtual double hist_to_tau(size_t g, size_t i, const History &history,
                             double T) const;
  /// Derivative of the map wrt to history
  virtual History d_hist_to_tau(size_t g, size_t i, const History &history,
                                double T) const;

  /// The scalar map
  virtual double hist_map(const History &history, double T) const = 0;
  /// The derivative of the scalar map
  virtual History d_hist_map(const History &history, double T) const = 0;
};

/// Slip strength rule where all systems evolve on a single scalar strength
class SlipSingleStrengthHardening : public SlipSingleHardening {
public:
  /// Request whatever history you will need
  virtual void populate_history(History &history) const;
  /// Setup history
  virtual void init_history(History &history) const;

  /// The rate of the history
  virtual History hist(const Symmetric &stress, const Orientation &Q,
                       const History &history, Lattice &L, double T,
                       const SlipRule &R) const;
  /// Derivative of the history wrt stress
  virtual History d_hist_d_s(const Symmetric &stress, const Orientation &Q,
                             const History &history, Lattice &L, double T,
                             const SlipRule &R) const;
  /// Derivative of the history wrt the history
  virtual History d_hist_d_h(const Symmetric &stress, const Orientation &Q,
                             const History &history, Lattice &L, double T,
                             const SlipRule &R) const;

  /// The scalar map
  virtual double hist_map(const History &history, double T) const;
  /// The derivative of the scalar map
  virtual History d_hist_map(const History &history, double T) const;

  /// Static (not evolving) strength
  virtual double static_strength(double T) const = 0;

  /// Setup the scalar
  virtual double init_strength() const = 0;

  /// Scalar evolution law
  virtual double hist_rate(const Symmetric &stress, const Orientation &Q,
                           const History &history, Lattice &L, double T,
                           const SlipRule &R) const = 0;
  /// Derivative of scalar law wrt stress
  virtual Symmetric d_hist_rate_d_stress(const Symmetric &stress,
                                         const Orientation &Q,
                                         const History &history, Lattice &L,
                                         double T, const SlipRule &R) const = 0;
  /// Derivative of scalar law wrt the scalar
  virtual double d_hist_rate_d_strength(const Symmetric &stress,
                                        const Orientation &Q,
                                        const History &history, Lattice &L,
                                        double T, const SlipRule &R) const = 0;
};

/// Slip strength rule where all systems evolve on a single scalar strength
class SlipDoubleStrengthHardening : public SlipSingleStrengthHardening {
public:
  /// Request whatever history you will need
  void populate_history(History &history) const override;
  /// Setup history
  void init_history(History &history) const override;
};

/// Slip strength rule where the single strength evolves with sum|dg|
class PlasticSlipHardening : public SlipSingleStrengthHardening {
public:
  /// Scalar evolution law
  virtual double hist_rate(const Symmetric &stress, const Orientation &Q,
                           const History &history, Lattice &L, double T,
                           const SlipRule &R) const;
  /// Derivative of scalar law wrt stress
  virtual Symmetric d_hist_rate_d_stress(const Symmetric &stress,
                                         const Orientation &Q,
                                         const History &history, Lattice &L,
                                         double T, const SlipRule &R) const;
  /// Derivative of scalar law wrt the scalar
  virtual double d_hist_rate_d_strength(const Symmetric &stress,
                                        const Orientation &Q,
                                        const History &history, Lattice &L,
                                        double T, const SlipRule &R) const;

  /// Prefactor
  virtual double hist_factor(double strength, Lattice &L, double T) const = 0;
  /// Derivative of the prefactor
  virtual double d_hist_factor(double strength, Lattice &L, double T) const = 0;
};

/// Everyone's favorite Voce model
class VoceSlipHardening : public PlasticSlipHardening {
public:
  /// Initialize with the saturated strength, the rate constant, and a constant
  /// strength
  VoceSlipHardening(std::shared_ptr<Interpolate> tau_sat,
                    std::shared_ptr<Interpolate> b,
                    std::shared_ptr<Interpolate> tau_0);

  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet &params);
  /// Default parameters
  static ParameterSet parameters();

  /// Setup the scalar
  virtual double init_strength() const;

  /// Static strength
  virtual double static_strength(double T) const;

  /// Prefactor
  virtual double hist_factor(double strength, Lattice &L, double T) const;
  /// Derivative of the prefactor
  virtual double d_hist_factor(double strength, Lattice &L, double T) const;

private:
  std::shared_ptr<Interpolate> tau_sat_;
  std::shared_ptr<Interpolate> b_;
  std::shared_ptr<Interpolate> tau_0_;
};

static Register<VoceSlipHardening> regVoceSlipHardening;

/// Everyone's favorite Voce model
class DoubleVoceSlipHardening : public PlasticSlipHardening {
public:
  DoubleVoceSlipHardening(std::shared_ptr<Interpolate> tau_sat_1,
                          std::shared_ptr<Interpolate> b_1,
                          std::shared_ptr<Interpolate> tau_sat_2,
                          std::shared_ptr<Interpolate> b_2,
                          std::shared_ptr<Interpolate> tau_0);

  /// Request whatever history you will need
  void populate_history(History &history) const override;
  /// Setup history
  void init_history(History &history) const override;

  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet &params);
  /// Default parameters
  static ParameterSet parameters();

  /// Setup the scalar
  virtual double init_strength() const;

  /// Static strength
  virtual double static_strength(double T) const;

  /// The scalar map
  double hist_map(const History &history, double T) const override;
  /// The derivative of the scalar map
  History d_hist_map(const History &history, double T) const override;

  /// Prefactor
  double hist_factor(const unsigned int s_idx, const History &hisotry,
                     Lattice &L, double T) const;
  /// Derivative of the prefactor
  double d_hist_factor(const unsigned int s_idx, const History &hisotry,
                       Lattice &L, double T) const;

  /// Scalar evolution law
  double hist_rate(const unsigned int s_idx, const Symmetric &stress,
                   const Orientation &Q, const History &history, Lattice &L,
                   double T, const SlipRule &R) const;

  /// Derivative of scalar law wrt stress
  Symmetric d_hist_rate_d_stress(const unsigned int s_idx,
                                 const Symmetric &stress, const Orientation &Q,
                                 const History &history, Lattice &L, double T,
                                 const SlipRule &R) const;

  /// Derivative of scalar law wrt the scalar
  double d_hist_rate_d_strength(const unsigned int s_idx,
                                const Symmetric &stress, const Orientation &Q,
                                const History &history, Lattice &L, double T,
                                const SlipRule &R) const;

  /// The rate of the history
  virtual History hist(const Symmetric &stress, const Orientation &Q,
                       const History &history, Lattice &L, double T,
                       const SlipRule &R) const override;

  /// Derivative of the history wrt the history
  History d_hist_d_h(const Symmetric &stress, const Orientation &Q,
                     const History &history, Lattice &L, double T,
                     const SlipRule &R) const override;

  /// Derivative of the history wrt stress
  History d_hist_d_s(const Symmetric &stress, const Orientation &Q,
                     const History &history, Lattice &L, double T,
                     const SlipRule &R) const override;

  /// Scalar evolution law
  double hist_rate(const Symmetric &stress, const Orientation &Q,
                   const History &history, Lattice &L, double T,
                   const SlipRule &R) const override {
    throw std::runtime_error(
        "original hist_rate overloaded, this should not be called");
    return 0;
  };

  /// Prefactor
  double hist_factor(double strength, Lattice &L, double T) const override {
    throw std::runtime_error(
        "original hist_factor overloaded, this should not be called");
    return 0;
  };

  /// Prefactor
  double d_hist_factor(double strength, Lattice &L, double T) const override {
    throw std::runtime_error(
        "original d_hist_factor overloaded, this should not be called");
    return 0;
  };

  /// Derivative of scalar law wrt the scalar
  double d_hist_rate_d_strength(const Symmetric &stress, const Orientation &Q,
                                const History &history, Lattice &L, double T,
                                const SlipRule &R) const override {
    throw std::runtime_error("original d_hist_rate_d_strength overloaded, this "
                             "should not be called");
    return 0;
  };
  /// Derivative of scalar law wrt stress
  Symmetric d_hist_rate_d_stress(const Symmetric &stress, const Orientation &Q,
                                 const History &history, Lattice &L, double T,
                                 const SlipRule &R) const override {
    throw std::runtime_error("original d_hist_rate_d_strength overloaded, this "
                             "should not be called");
    Symmetric a;
    return a;
  };

private:
  std::shared_ptr<Interpolate> tau_sat_1_;
  std::shared_ptr<Interpolate> b_1_;
  std::shared_ptr<Interpolate> tau_sat_2_;
  std::shared_ptr<Interpolate> b_2_;
  std::shared_ptr<Interpolate> tau_0_;
};

static Register<DoubleVoceSlipHardening> regDoubleVoceSlipHardening;
} // namespace neml

#endif // SLIPHARDEN_H
