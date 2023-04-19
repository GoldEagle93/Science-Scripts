import parmed as pmd  # type: ignore
from meld.system import patchers
from meld import util  
from parmed import unit as u
import numpy as np  # type: ignore

class DummyPatcher(patchers.PatcherBase):
    def __init__(self, n_dummies, coord_dummies):
        """
        Patch system to include extra dummy atoms

        Parameters
        ----------
        n_dummies: int
            number of unique dummies
        """
        self.n_dummies = n_dummies
        self.coord_dummies = coord_dummies
        self.resids = []

    def patch(self, top_string, crd_string):
        INTOP = "in.top"
        INRST = "in.rst"
        OUTTOP = "out.top"
        OUTRST = "out.rst"

        with util.in_temp_dir():
            with open(INTOP, "wt") as outfile:
                outfile.write(top_string)
            with open(INRST, "wt") as outfile:
                outfile.write(crd_string)

            base = pmd.load_file(INTOP)
            crd = pmd.load_file(INRST)
            base.coordinates = crd.coordinates

            # create a new structure to add our dummy atoms to
            parm = pmd.Structure()

            # add in atom type for our dummy particles
            atype = pmd.AtomType("SDUM", 0, mass=12000000000.0, charge=0.0)
            atype.set_lj_params(0.0, 0.0)

            for i in range(self.n_dummies):
                a1 = pmd.Atom(
                    name="S{}".format(i),
                    atomic_number=atype.atomic_number,
                    type=str(atype),
                    charge=atype.charge,
                    mass=atype.mass,
                    solvent_radius=1.0,
                    screen=0.5,
                )
                a1.atom_type = atype

                parm.add_atom(a1, resname="SDM", resnum=i)

            # we add noise here because we'll get NaN if the particles ever
            # end up exactly on top of each other
            parm.positions = np.array(self.coord_dummies)

            # combine the old system with the new dummy atoms
            comb = base + parm
            last_index = comb.residues[-1].idx
            self.resids = list(range(last_index - self.n_dummies + 1, last_index + 1))

            comb.write_parm(OUTTOP)
            comb.write_rst7(OUTRST)
            with open(OUTTOP, "rt") as infile:
                top_string = infile.read()
            with open(OUTRST, "rt") as infile:
                crd_string = infile.read()
            with open('topology.top','w') as out:
                out.write(top_string)
            with open('topology.crd','w') as out:
                out.write(crd_string)
        return top_string, crd_string
