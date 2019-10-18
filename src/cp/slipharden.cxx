#include "slipharden.h"

namespace neml {

double SlipSingleHardening::hist_to_tau(size_t g, size_t i,
                                        const History &history,
                                        double T) const {
  return hist_map(history, T);
}

History SlipSingleHardening::d_hist_to_tau(size_t g, size_t i,
                                           const History &history,
                                           double T) const {
  return d_hist_map(history, T);
}

void SlipSingleStrengthHardening::populate_history(History &history) const {
  history.add<double>("strength");
}

void SlipSingleStrengthHardening::init_history(History &history) const {
  history.get<double>("strength") = init_strength();
}

History SlipSingleStrengthHardening::hist(const Symmetric &stress,
                                          const Orientation &Q,
                                          const History &history, Lattice &L,
                                          double T, const SlipRule &R) const {
  History res = history.copy_blank();
  res.get<double>("strength") = hist_rate(stress, Q, history, L, T, R);

  return res;
}

History SlipSingleStrengthHardening::d_hist_d_s(const Symmetric &stress,
                                                const Orientation &Q,
                                                const History &history,
                                                Lattice &L, double T,
                                                const SlipRule &R) const {
  History res = history.derivative<Symmetric>();
  res.get<Symmetric>("strength") =
      d_hist_rate_d_stress(stress, Q, history, L, T, R);
  return res;
}

History SlipSingleStrengthHardening::d_hist_d_h(const Symmetric &stress,
                                                const Orientation &Q,
                                                const History &history,
                                                Lattice &L, double T,
                                                const SlipRule &R) const {
  History res = history.derivative<History>();
  res.get<double>("strength_strength") =
      d_hist_rate_d_strength(stress, Q, history, L, T, R);
  return res;
}

double SlipSingleStrengthHardening::hist_map(const History &history,
                                             double T) const {
  return history.get<double>("strength") + static_strength(T);
}

History SlipSingleStrengthHardening::d_hist_map(const History &history,
                                                double T) const {
  History res = history.derivative<double>();
  res.get<double>("strength") = 1.0;
  return res;
}

double PlasticSlipHardening::hist_rate(const Symmetric &stress,
                                       const Orientation &Q,
                                       const History &history, Lattice &L,
                                       double T, const SlipRule &R) const {
  double strength = history.get<double>("strength");

  return hist_factor(strength, L, T) * R.sum_slip(stress, Q, history, L, T);
}

Symmetric PlasticSlipHardening::d_hist_rate_d_stress(const Symmetric &stress,
                                                     const Orientation &Q,
                                                     const History &history,
                                                     Lattice &L, double T,
                                                     const SlipRule &R) const {
  double strength = history.get<double>("strength");
  return hist_factor(strength, L, T) *
         R.d_sum_slip_d_stress(stress, Q, history, L, T);
}

double PlasticSlipHardening::d_hist_rate_d_strength(const Symmetric &stress,
                                                    const Orientation &Q,
                                                    const History &history,
                                                    Lattice &L, double T,
                                                    const SlipRule &R) const {
  double strength = history.get<double>("strength");
  History dhist = R.d_sum_slip_d_hist(stress, Q, history, L, T);
  double dstrength = dhist.get<double>("strength");
  return d_hist_factor(strength, L, T) * R.sum_slip(stress, Q, history, L, T) +
         hist_factor(strength, L, T) * dstrength;
}

VoceSlipHardening::VoceSlipHardening(std::shared_ptr<Interpolate> tau_sat,
                                     std::shared_ptr<Interpolate> b,
                                     std::shared_ptr<Interpolate> tau_0)
    : tau_sat_(tau_sat), b_(b), tau_0_(tau_0) {}

std::string VoceSlipHardening::type() { return "VoceSlipHardening"; }

std::unique_ptr<NEMLObject>
VoceSlipHardening::initialize(ParameterSet &params) {
  return neml::make_unique<VoceSlipHardening>(
      params.get_object_parameter<Interpolate>("tau_sat"),
      params.get_object_parameter<Interpolate>("b"),
      params.get_object_parameter<Interpolate>("tau_0"));
}

ParameterSet VoceSlipHardening::parameters() {
  ParameterSet pset(VoceSlipHardening::type());

  pset.add_parameter<NEMLObject>("tau_sat");
  pset.add_parameter<NEMLObject>("b");
  pset.add_parameter<NEMLObject>("tau_0");

  return pset;
}

double VoceSlipHardening::init_strength() const { return 0.0; }

double VoceSlipHardening::static_strength(double T) const {
  return tau_0_->value(T);
}

double VoceSlipHardening::hist_factor(double strength, Lattice &L,
                                      double T) const {
  double tau_sat = tau_sat_->value(T);
  double b = b_->value(T);

  std::cout << "tau_sat " << tau_sat << " b " << b << " strength " << strength
            << std::endl;
  return b * (tau_sat - strength);
}

double VoceSlipHardening::d_hist_factor(double strength, Lattice &L,
                                        double T) const {
  double b = b_->value(T);

  return -b;
}

DoubleVoceSlipHardening::DoubleVoceSlipHardening(
    std::shared_ptr<Interpolate> tau_sat_1, std::shared_ptr<Interpolate> b_1,
    std::shared_ptr<Interpolate> tau_sat_2, std::shared_ptr<Interpolate> b_2,
    std::shared_ptr<Interpolate> tau_0)
    : tau_sat_1_(tau_sat_1), b_1_(b_1), tau_sat_2_(tau_sat_2), b_2_(b_2),
      tau_0_(tau_0) {}

std::string DoubleVoceSlipHardening::type() {
  return "DoubleVoceSlipHardening";
}

std::unique_ptr<NEMLObject>
DoubleVoceSlipHardening::initialize(ParameterSet &params) {
  return neml::make_unique<DoubleVoceSlipHardening>(
      params.get_object_parameter<Interpolate>("tau_sat_1"),
      params.get_object_parameter<Interpolate>("b_1"),
      params.get_object_parameter<Interpolate>("tau_sat_2"),
      params.get_object_parameter<Interpolate>("b_2"),
      params.get_object_parameter<Interpolate>("tau_0"));
}

ParameterSet DoubleVoceSlipHardening::parameters() {
  ParameterSet pset(DoubleVoceSlipHardening::type());

  pset.add_parameter<NEMLObject>("tau_sat_1");
  pset.add_parameter<NEMLObject>("b_1");
  pset.add_parameter<NEMLObject>("tau_sat_2");
  pset.add_parameter<NEMLObject>("b_2");
  pset.add_parameter<NEMLObject>("tau_0");

  return pset;
}

void DoubleVoceSlipHardening::populate_history(History &history) const {
  history.add<double>("strength1");
  history.add<double>("strength2");
}

void DoubleVoceSlipHardening::init_history(History &history) const {
  history.get<double>("strength1") = init_strength();
  history.get<double>("strength2") = init_strength();
}

double DoubleVoceSlipHardening::init_strength() const { return 0.0; }

double DoubleVoceSlipHardening::static_strength(double T) const {
  return tau_0_->value(T);
}

double DoubleVoceSlipHardening::hist_map(const History &history,
                                         double T) const {
  double s = static_strength(T);
  s += history.get<double>("strength1");
  s += history.get<double>("strength2");
  return s;
}

History DoubleVoceSlipHardening::d_hist_map(const History &history,
                                            double T) const {
  History res = history.derivative<double>();
  res.get<double>("strength1") = 1.0;
  res.get<double>("strength2") = 1.0;
  return res;
}

double DoubleVoceSlipHardening::hist_factor(const unsigned int s_idx,
                                            const History &hisotry, Lattice &L,
                                            double T) const {
  double tau_sat, b, strength;
  if (s_idx == 1) {
    tau_sat = tau_sat_1_->value(T);
    b = b_1_->value(T);
    strength = hisotry.get<double>("strength1");
  } else if (s_idx == 2) {
    tau_sat = tau_sat_2_->value(T);
    b = b_2_->value(T);
    strength = hisotry.get<double>("strength2");
  } else {
    throw std::range_error("hist_factor: value of s_idx out of range");
  }
  std::cout << "tau_sat " << tau_sat << " b " << b << " strength " << strength
            << std::endl;
  return b * (tau_sat - strength);
}

double DoubleVoceSlipHardening::d_hist_factor(const unsigned int s_idx,
                                              const History &hisotry,
                                              Lattice &L, double T) const {
  double b;
  if (s_idx == 1) {
    b = b_1_->value(T);
  } else if (s_idx == 2) {
    b = b_2_->value(T);
  } else {
    throw std::range_error("d_hist_factor: value of s_idx out of range");
  }

  return -b;
}

double DoubleVoceSlipHardening::hist_rate(const unsigned int s_idx,
                                          const Symmetric &stress,
                                          const Orientation &Q,
                                          const History &history, Lattice &L,
                                          double T, const SlipRule &R) const {

  return hist_factor(s_idx, history, L, T) *
         R.sum_slip(stress, Q, history, L, T);
}

Symmetric DoubleVoceSlipHardening::d_hist_rate_d_stress(
    const unsigned int s_idx, const Symmetric &stress, const Orientation &Q,
    const History &history, Lattice &L, double T, const SlipRule &R) const {

  return hist_factor(s_idx, history, L, T) *
         R.d_sum_slip_d_stress(stress, Q, history, L, T);
}

double DoubleVoceSlipHardening::d_hist_rate_d_strength(
    const unsigned int s_idx, const Symmetric &stress, const Orientation &Q,
    const History &history, Lattice &L, double T, const SlipRule &R) const {

  History dhist = R.d_sum_slip_d_hist(stress, Q, history, L, T);
  double dstrength;
  if (s_idx == 1) {
    dstrength = dhist.get<double>("strength1");
  } else if (s_idx == 2) {
    dstrength = dhist.get<double>("strength2");
  } else {
    throw std::range_error("d_hist_factor: value of s_idx out of range");
  }

  return d_hist_factor(s_idx, history, L, T) *
             R.sum_slip(stress, Q, history, L, T) +
         hist_factor(s_idx, history, L, T) * dstrength;
}

History DoubleVoceSlipHardening::hist(const Symmetric &stress,
                                      const Orientation &Q,
                                      const History &history, Lattice &L,
                                      double T, const SlipRule &R) const {
  History res = history.copy_blank();
  res.get<double>("strength1") = hist_rate(1, stress, Q, history, L, T, R);
  res.get<double>("strength2") = hist_rate(2, stress, Q, history, L, T, R);
  return res;
}

History DoubleVoceSlipHardening::d_hist_d_h(const Symmetric &stress,
                                            const Orientation &Q,
                                            const History &history, Lattice &L,
                                            double T, const SlipRule &R) const {
  History res = history.derivative<History>();
  res.get<double>("strength1_strength1") =
      d_hist_rate_d_strength(1, stress, Q, history, L, T, R);
  res.get<double>("strength2_strength2") =
      d_hist_rate_d_strength(2, stress, Q, history, L, T, R);
  res.get<double>("strength1_strength2") = 0;
  res.get<double>("strength2_strength1") = 0;
  return res;
}

History DoubleVoceSlipHardening::d_hist_d_s(const Symmetric &stress,
                                            const Orientation &Q,
                                            const History &history, Lattice &L,
                                            double T, const SlipRule &R) const {
  History res = history.derivative<Symmetric>();
  res.get<Symmetric>("strength1") =
      d_hist_rate_d_stress(1, stress, Q, history, L, T, R);
  res.get<Symmetric>("strength2") =
      d_hist_rate_d_stress(2, stress, Q, history, L, T, R);
  return res;
}

} // namespace neml
