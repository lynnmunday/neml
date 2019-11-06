#!/usr/bin/env python3

import sys
sys.path.append('..')

from neml import models, elasticity, drivers, damage, creep

import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
  # Elasticity
  E = 180000.0
  nu = 0.3

  # Power law creep parameters
  n = 6.0
  A = 1000.0

  # Damage parameters
  xi = 0.5
  phi = 2.5
  B = 1.0e14

  # Loading
  srate = 1000000.0

  # Setup model
  emodel = elasticity.IsotropicLinearElasticModel(E, "youngs", nu, "poissons")
  bmodel = models.SmallStrainElasticity(emodel)
  scmodel = creep.NormalizedPowerLawCreep(A, n)
  cfmodel = creep.J2CreepModel(scmodel)
  cmodel = models.SmallStrainCreepPlasticity(emodel, bmodel, cfmodel)
  model = damage.ModularCreepDamageModel_sd(emodel, B, xi, phi,
      damage.VonMisesEffectiveStress(), cmodel)

  # Computed life
  srange = np.linspace(A/30,A/20, 10)

  tfs = (srange/B)**(-xi) / (1+phi)
  
  slife = []
  erates = []
  for s,tf in zip(srange, tfs):
    res = drivers.creep(model, s, srate, tf * 1.25, verbose = False,
        logspace = False, check_dmg = True, nsteps = 1000)
    slife.append(res['rtime'][-1])
    erates.append([res['rtime'],res['rrate']])
  
  plt.figure()
  plt.plot(srange, tfs, label = "Analytic formula")
  plt.plot(srange, np.array(slife), label = "ModularCreepDamageModel")
  plt.xlabel("Stress")
  plt.ylabel("Failure time")
  plt.legend(loc = 'best')
  plt.tight_layout()
  plt.show()

  plt.figure()
  for stress, (times, rates) in zip(srange,erates):
    l, = plt.loglog(times, rates)
    ws = 1.0 - (1.0 - (1+phi)*(stress/B)**xi * times)**(1.0/(1.0+phi))
    rates2 = (stress/(A*(1.0-ws)))**(n)
    plt.loglog(times, rates2, color = l.get_color(), ls = '--')

  plt.xlabel("Time")
  plt.ylabel("Rate")
  plt.tight_layout()

  plt.show()
