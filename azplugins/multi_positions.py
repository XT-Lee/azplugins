import hoomd
import numpy
import hoomd.azplugins as azplugins


def set_multi_positions(traps_pos,k,r_cut,system):
    R""" Add an array of harmonic restraining potentials based on specified positions
    Args:
        traps_pos(float array_like): an array of n*3 which defines n potentials you will apply to particles
        k (float or array_like): Force constant, isotropic or in each of *x*-, *y*-, and *z*-directions
        r_cut (non-negative float): Force constant, a scalar
        system(object): the object created and initialized by hoomd.init

    The Hamiltonian is augmented with a harmonic potential calculated based on the distance between the current and
    initial positions of particles. This effectively allows particles to have nearly fixed position while still retaining
    their degrees of freedom. The potential has the following form:
    .. math::
        :nowrap:
        \begin{equation*}
        V(\mathbf{r}) = \frac{1}{2} \mathbf{k} \mathbf{\Delta r} \mathbf{\Delta r}^T  where \Delta r < r_cut, 
        V(\mathbf{r}) = \frac{1}{2} \mathbf{k} \mathbf{r_cut} \mathbf{r_cut}^T  where \Delta r >= r_cut, 
        \end{equation*}
    The strength of the potential depends on a spring constant :math:`\mathbf{k}` so that the particle position can be
    restrained in any of the three coordinate directions. If :math:`\mathbf{k}` is very large, the particle position is
    essentially constrained. However, because the particles retain their degrees of freedom, shorter integration timesteps
    must be taken for large :math:`\mathbf{k}` to maintain stability.
    .. note::
        The displacement is calculated using the minimum image convention.
    Examples::
        import hoomd.azplugins.multi_positions as mp
        sys=hoomd.init.create_lattice(unitcell=hoomd.lattice.sq(a=5), n=9);
        traps_pos=[[10,10,0],[-10,-10,0]];
        mp.set_multi_positions(traps_pos,k=100.0,r_cut=6.0,system=sys)
    .. warning::
        Virial calculation is not implemented because the particles are tethered to fixed positions. A warning will be raised
        if any calls to :py:class:`hoomd.analyze.log` are made because the logger always requests the virial flags. However,
        this warning can be safely ignored if the pressure (tensor) is not being logged or the pressure is not of interest.
    """
    #create a list of trap positions
    #generate an array for following for loop
    traps_shape = numpy.asarray(traps_pos);
    num_traps=traps_shape.shape[0];
    arr_traps = range(num_traps);
    #generate an array for following for loop
    snap = system.take_snapshot(all=True)
    num_particles=len(snap.particles.position[:,0])
    arr_particles = range(num_particles);

    #generate an array of traps for object 'position'
    group=hoomd.group.all()
    pos1 = azplugins.restrain.position(group=group,k=0.0,r_cut=0.0);
    traps_list=[pos1]*num_traps;

    #define the traps
    #and apply them to all the particles
    for (i) in  arr_traps :
        traps_list[i]=azplugins.restrain.position(group=group,k=k,r_cut=r_cut);
        for (j) in  arr_particles :
            traps_list[i].set_position(int(j),traps_pos[i][:]);