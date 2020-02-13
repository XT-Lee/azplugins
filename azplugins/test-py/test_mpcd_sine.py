# Copyright (c) 2018-2019, Michael P. Howard
# This file is part of the azplugins project, released under the Modified BSD License.

# Maintainer: astatt

import unittest
import numpy as np
import hoomd
from hoomd import md
from hoomd import mpcd

try:
    from hoomd import azplugins
    import hoomd.azplugins.mpcd
except ImportError:
    import azplugins
    import azplugins.mpcd
import unittest

# unit tests for mpcd sine channel  geometry
class mpcd_sine_test(unittest.TestCase):
    def setUp(self):

        hoomd.context.initialize()

        # set the decomposition in z for mpi builds
        if hoomd.comm.get_num_ranks() > 1:
            hoomd.comm.decomposition(nz=2)

        # default testing configuration
        hoomd.init.read_snapshot(hoomd.data.make_snapshot(N=0, box=hoomd.data.boxdim(L=15.)))

        # initialize the system from the starting snapshot
        # test vertical, diagonal, and horizontal collisions to wall
        snap = mpcd.data.make_snapshot(N=3)
        snap.particles.position[:] = [[1.,-3.0,-3.8],[3.5,0.,3.],[-4.2,5.1,-2.2]]
        snap.particles.velocity[:] = [[0.,0.,-1.],[1.,0.,0.],[-1.,-1.,-1.]]
        self.s = mpcd.init.read_snapshot(snap)

        mpcd.integrator(dt=0.1)

    # test creation can happen (with all parameters set)
    def test_create(self):
        azplugins.mpcd.sine(H=4., h=2. ,p=1, V=0.1, boundary="no_slip")

    # test for setting parameters
    def test_set_params(self):
        sine = azplugins.mpcd.sine(H=4.,h=2. ,p=1)
        self.assertAlmostEqual(sine.H, 4.)
        self.assertAlmostEqual(sine.V, 0.)
        self.assertEqual(sine.boundary, "no_slip")
        self.assertAlmostEqual(sine._cpp.geometry.getHwide(), 4.)
        self.assertAlmostEqual(sine._cpp.geometry.getVelocity(), 0.)
        self.assertEqual(sine._cpp.geometry.getBoundaryCondition(), mpcd._mpcd.boundary.no_slip)

        # change H and also ensure other parameters stay the same
        sine.set_params(H=2.)
        self.assertAlmostEqual(sine.H, 2.)
        self.assertAlmostEqual(sine.V, 0.)
        self.assertEqual(sine.boundary, "no_slip")
        self.assertAlmostEqual(sine._cpp.geometry.getHwide(), 2.)
        self.assertAlmostEqual(sine._cpp.geometry.getVelocity(), 0.)
        self.assertEqual(sine._cpp.geometry.getBoundaryCondition(), mpcd._mpcd.boundary.no_slip)

        # change V
        sine.set_params(V=0.1)
        self.assertAlmostEqual(sine.V, 0.1)
        self.assertAlmostEqual(sine._cpp.geometry.getVelocity(), 0.1)

        # change BCs
        sine.set_params(boundary="slip")
        self.assertEqual(sine.boundary, "slip")
        self.assertEqual(sine._cpp.geometry.getBoundaryCondition(), mpcd._mpcd.boundary.slip)

    # test for invalid boundary conditions being set
    def test_bad_boundary(self):
        sine = azplugins.mpcd.sine(H=4.,h=2. ,p=1)
        sine.set_params(boundary="no_slip")
        sine.set_params(boundary="slip")

        with self.assertRaises(ValueError):
            sine.set_params(boundary="invalid")

    # test basic stepping behavior with no slip boundary conditions
    def test_step_noslip(self):
        azplugins.mpcd.sine(H=4.,h=2. ,p=1, boundary='no_slip')

        # take one step, no particle hits the wall
        hoomd.run(1)
        snap = self.s.take_snapshot()
        if hoomd.comm.get_rank() == 0:
            np.testing.assert_array_almost_equal(snap.particles.position[0], [1,-3.0,-3.9])
            np.testing.assert_array_almost_equal(snap.particles.velocity[0], [0.,0.,-1.])
            np.testing.assert_array_almost_equal(snap.particles.position[1], [3.6,0.0,3.0])
            np.testing.assert_array_almost_equal(snap.particles.velocity[1], [1.,0,0])
            np.testing.assert_array_almost_equal(snap.particles.position[2], [-4.3,5.0,-2.3])
            np.testing.assert_array_almost_equal(snap.particles.velocity[2], [-1.,-1.,-1.])

        # take another step where  particle 1 will now hit the wall vertically
        # point of wall contact is z=-(cos(2*pi/15.)+3) = -3.913545, remaining integration time is 0.086455
        # so resulting position is -3.913545+0.086455=-3.82709
        hoomd.run(1)
        snap = self.s.take_snapshot()
        if hoomd.comm.get_rank() == 0:
            np.testing.assert_array_almost_equal(snap.particles.position[0], [1,-3.0,-3.82709])
            np.testing.assert_array_almost_equal(snap.particles.velocity[0], [0.,0.,+1.])
            np.testing.assert_array_almost_equal(snap.particles.position[1], [3.7,0.0,3.0])
            np.testing.assert_array_almost_equal(snap.particles.velocity[1], [1.,0.,0.])
            np.testing.assert_array_almost_equal(snap.particles.position[2], [-4.4,4.9,-2.4])
            np.testing.assert_array_almost_equal(snap.particles.velocity[2], [-1.,-1.,-1.])

        # take another step, where  particle 2 will now hit the wall horizontally
        # dt = 0.05, particle travels exactly 0.05 inside, and then gets projected back right onto the wall
        hoomd.run(1)
        snap = self.s.take_snapshot()
        if hoomd.comm.get_rank() == 0:
            np.testing.assert_array_almost_equal(snap.particles.position[0], [1,-3.0,-3.72709])
            np.testing.assert_array_almost_equal(snap.particles.velocity[0], [0,0,+1.])
            np.testing.assert_array_almost_equal(snap.particles.position[1], [3.7,0.0,3.0])
            np.testing.assert_array_almost_equal(snap.particles.velocity[1], [-1.,0.,0.])
            np.testing.assert_array_almost_equal(snap.particles.position[2], [-4.5,4.8,-2.5])
            np.testing.assert_array_almost_equal(snap.particles.velocity[2], [-1.,-1.,-1.])

        # take another step, no particle collides, check for spurious collisions
        hoomd.run(1)
        snap = self.s.take_snapshot()
        if hoomd.comm.get_rank() == 0:
            np.testing.assert_array_almost_equal(snap.particles.position[0], [1,-3.0,-3.62709])
            np.testing.assert_array_almost_equal(snap.particles.velocity[0], [0,0,+1.])
            np.testing.assert_array_almost_equal(snap.particles.position[1], [3.6,0.0,3.0])
            np.testing.assert_array_almost_equal(snap.particles.velocity[1], [-1.,0.,0.])
            np.testing.assert_array_almost_equal(snap.particles.position[2], [-4.6,4.7,-2.6])
            np.testing.assert_array_almost_equal(snap.particles.velocity[2], [-1.,-1.,-1.])

        # take another step, last particle collides
        # wall intersection: -4.636956 4.663044 -2.63696 (calculated with python) dt = 0.063042
        # position -4.636956+0.06304 4.663044+0.063042 -2.63696+0.063042
        hoomd.run(1)
        snap = self.s.take_snapshot()
        if hoomd.comm.get_rank() == 0:
            np.testing.assert_array_almost_equal(snap.particles.position[0], [1,-3.0,-3.52709])
            np.testing.assert_array_almost_equal(snap.particles.velocity[0], [0,0,+1.])
            np.testing.assert_array_almost_equal(snap.particles.position[1], [3.5,0.0,3.0])
            np.testing.assert_array_almost_equal(snap.particles.velocity[1], [-1.,0.,0.])
            np.testing.assert_array_almost_equal(snap.particles.position[2], [-4.573913,  4.726087, -2.573919])
            np.testing.assert_array_almost_equal(snap.particles.velocity[2], [1.,1.,1.])

    #same as test above except for slip -> velcities differ
    def test_step_slip(self):
        azplugins.mpcd.sine(H=4.,h=2. ,p=1, boundary="slip")

        # take one step, no particle hits the wall
        hoomd.run(1)
        snap = self.s.take_snapshot()
        if hoomd.comm.get_rank() == 0:
            np.testing.assert_array_almost_equal(snap.particles.position[0], [1,-3.0,-3.9])
            np.testing.assert_array_almost_equal(snap.particles.velocity[0], [0.,0.,-1.])
            np.testing.assert_array_almost_equal(snap.particles.position[1], [3.6,0.0,3.0])
            np.testing.assert_array_almost_equal(snap.particles.velocity[1], [1.,0,0])
            np.testing.assert_array_almost_equal(snap.particles.position[2], [-4.3,5.0,-2.3])
            np.testing.assert_array_almost_equal(snap.particles.velocity[2], [-1.,-1.,-1.])

        # take another step where  particle 1 will now hit the wall vertically
        # point of contact with wall same test before, but velocity needs to be reflected.
        # point of wall contact is z=-(cos(2*pi/15.)+3) = -3.913545, remaining integration time is 0.086455
        # so resulting position is -3.913545+0.086455*v_z=-3.82709
        # B for surface normal is -0.17037344664, so v_y = 0,
        # v_x = 0 + 2B/(B^2+1) = 0 -0.33113500075
        # v_z =  -1 + 2/(B^2+1) = -1 + 1.94358338862
        # now, new pos = contact point wall + dt*v
        hoomd.run(1)
        snap = self.s.take_snapshot()
        if hoomd.comm.get_rank() == 0:
            np.testing.assert_array_almost_equal(snap.particles.position[0], [0.971372,-3.0,-3.831968])
            np.testing.assert_array_almost_equal(snap.particles.velocity[0], [-0.331135,0.,0.943583])
            np.testing.assert_array_almost_equal(snap.particles.position[1], [3.7,0.0,3.0])
            np.testing.assert_array_almost_equal(snap.particles.velocity[1], [1.,0.,0.])
            np.testing.assert_array_almost_equal(snap.particles.position[2], [-4.4,4.9,-2.4])
            np.testing.assert_array_almost_equal(snap.particles.velocity[2], [-1.,-1.,-1.])

        # one more step, second particle collides
        # B = 0.418879 ( x0 approx 3.7)
        # velocities: v_y = 0
        # v_x = 1 - 2*B^2/(B^2+1) =  0.70146211038
        # v_z = 0 -2B/(B^2+1) =  -0.71270674733
        hoomd.run(1)
        snap = self.s.take_snapshot()
        if hoomd.comm.get_rank() == 0:
            np.testing.assert_array_almost_equal(snap.particles.position[0], [0.9382585,-3.0,-3.7376097])
            np.testing.assert_array_almost_equal(snap.particles.velocity[0], [-0.331135,0.,0.943583])
            np.testing.assert_array_almost_equal(snap.particles.position[1], [3.785073, 0.,2.964365])
            np.testing.assert_array_almost_equal(snap.particles.velocity[1], [0.70146211038,0.,-0.71270674733])
            np.testing.assert_array_almost_equal(snap.particles.position[2], [-4.5,4.8,-2.5])
            np.testing.assert_array_almost_equal(snap.particles.velocity[2], [-1.,-1.,-1.])

        # two more steps, last particle collides
        # B = 0.390301  (x0 approx -4.6)
        # velocities: v_y = -1
        # v_x = -1 - 2*B(-B-1)/(B^2+1) = -0.05819760480217273
        # v_z = -1 -2(-B-1)/(B^2+1) = 1.4130155833518931
        hoomd.run(2)
        snap = self.s.take_snapshot()
        if hoomd.comm.get_rank() == 0:
            np.testing.assert_array_almost_equal(snap.particles.position[2], [-4.640625,4.600002,-2.547881])
            np.testing.assert_array_almost_equal(snap.particles.velocity[2], [-0.05819760480217273,-1.,1.4130155833518931])

    # test that setting the sine size too large raises an error
    def test_validate_box(self):
        # initial configuration is invalid
        sine = azplugins.mpcd.sine(H=10.,h=2. ,p=1)
        with self.assertRaises(RuntimeError):
            hoomd.run(1)

        # now it should be valid
        sine.set_params(H=4.,h=2. ,p=1)
        hoomd.run(2)

        # make sure we can invalidate it again
        sine.set_params(H=10.,h=2. ,p=1)
        with self.assertRaises(RuntimeError):
            hoomd.run(1)

    # test that particles out of bounds can be caught
    def test_out_of_bounds(self):
        sine = azplugins.mpcd.sine(H=2., h=1. ,p=1)
        with self.assertRaises(RuntimeError):
            hoomd.run(1)

        sine.set_params(H=5.,h=2. ,p=1)
        hoomd.run(1)


    # test that virtual particle filler can be attached, removed, and updated
    def test_filler(self):
        # initialization of a filler
        sine = azplugins.mpcd.sine(H=5.,h=2. ,p=1)
        sine.set_filler(density=5., kT=1.0, seed=42, type='A')
        self.assertTrue(sine._filler is not None)

        # run should be able to setup the filler, although this all happens silently
        hoomd.run(1)

        # changing the geometry should still be OK with a run
        sine.set_params(V=1.0)
        hoomd.run(1)

        # changing filler should be allowed
        sine.set_filler(density=10., kT=1.5, seed=7)
        self.assertTrue(sine._filler is not None)
        hoomd.run(1)

        # assert an error is raised if we set a bad particle type
        with self.assertRaises(RuntimeError):
            sine.set_filler(density=5., kT=1.0, seed=42, type='B')

        # assert an error is raised if we set a bad density
        with self.assertRaises(RuntimeError):
            sine.set_filler(density=-1.0, kT=1.0, seed=42)

        # removing the filler should still allow a run
        sine.remove_filler()
        self.assertTrue(sine._filler is None)
        hoomd.run(1)


    def tearDown(self):
        del self.s

if __name__ == '__main__':
    unittest.main(argv = ['test.py', '-v'])
