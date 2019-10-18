#!/usr/bin/env python3

import sys
sys.path.append('../..')

import numpy as np

from neml.cp import crystallography, slipharden, sliprules, inelasticity, kinematics, singlecrystal
from neml.math import rotations, tensors, nemlmath
from neml import elasticity

import matplotlib.pyplot as plt

def run(t0, ts1, b1, ts2, b2, nsteps, std_voce = False ):
   L = np.array([[-0.5,0,0],[0,1.0,0],[0,0,-0.5]])
   erate = 1.0e-4
   steps = nsteps
   emax = 0.005

   E = 100000.0
   nu = 0.3





   g0 = 1.0
   n = 12.0

   # Setup
   L *= erate
   dt = emax / steps / erate

   d = nemlmath.sym(0.5*(L+L.T))
   w = nemlmath.skew(0.5*(L-L.T))

   if not std_voce:
     strengthmodel = slipharden.DoubleVoceSlipHardening(ts1, b1, ts2, b2, t0)
   else:
     strengthmodel = slipharden.VoceSlipHardening(ts1, b1,  t0)

   slipmodel = sliprules.PowerLawSlipRule(strengthmodel, g0, n)
   imodel = inelasticity.AsaroInelasticity(slipmodel)
   emodel = elasticity.IsotropicLinearElasticModel(E, "youngs", nu, "poissons")
   kmodel = kinematics.StandardKinematicModel(emodel, imodel)

   lattice = crystallography.CubicLattice(1.0)
   lattice.add_slip_system([1,1,0],[1,1,1])

   model = singlecrystal.SingleCrystalModel(kmodel, lattice, verbose = False)

   T = 300.0

   h_n = model.init_store()

   s_n = np.zeros((6,))

   d_n = np.zeros((6,))
   w_n = np.zeros((3,))

   u_n = 0.0
   p_n = 0.0
   t_n = 0.0
   T_n = T

   e = [0.0]
   s = [0.0]

   for i in range(steps):
     d_np1 = d_n + d * dt
     w_np1 = w_n + w * dt
     t_np1 = t_n + dt
     T_np1 = T_n

     s_np1, h_np1, A_np1, B_np1, u_np1, p_np1 = model.update_ld_inc(d_np1, d_n, w_np1, w_n, T_np1, T_n, t_np1, t_n, s_n, h_n, u_n, p_n)

     e.append(d_np1[1])
     s.append(s_np1[1])

     d_n = np.copy(d_np1)
     w_n = np.copy(w_np1)
     s_n = np.copy(s_np1)
     h_n = np.copy(h_np1)
     t_n = t_np1
     T_n = T_np1
     u_n = u_np1
     p_n = p_np1

   print(s[-1])
   plt.plot(e, s)



if __name__ == "__main__":
  t0 = 50.0


  ts1 = 50.0
  b1 = 100.0
  ts2 = 0.0
  b2 = 0.0
  run(t0, ts1, b1, ts2, b2, 100 )


  ts1 = 0.0
  b1 = 0.0
  ts2 = 50.0
  b2 = 100.0
  run(t0, ts1, b1, ts2, b2, 10000 )

  ts1 = 50.0
  b1 = 100.0
  ts2 = 0.0
  b2 = 0.0
  run(t0, ts1, b1, ts2, b2, 100, std_voce = True )

  plt.show()
