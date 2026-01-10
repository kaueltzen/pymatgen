from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pymatgen.io.lobster import LobsterNeighbors
    from pymatgen.util.typing import PathLike


class LobsterNeighborsVisualizer:
    """
    Class to create VESTA files with bonds defined from LOBSTER bonding analysis.
    TODO bond_tol file_name in init or plotting function? see Plotter classes of pmg
    """

    def __init__(
        self, lobster_neighbors: LobsterNeighbors, bond_tol: float = 0.1, file_name: str = "output.vesta"
    ) -> None:
        """
        Args:
            lobsterneighbors : pymatgen.io.lobster.LobsterNeighbors
                LobsterNeighbors object used to determine bonding.
            bond_tol : float
                Tolerance added/subtracted from bond lengths for VESTA SBOND entries.
            file_name : str
                Output VESTA file name.
        """
        self.bond_tol = bond_tol
        self.file_name = file_name
        self.lobster_neighbors = lobster_neighbors
        self.bonds = self._get_relevant_bonds()

    @classmethod
    def from_lobster_calc(
        cls,
        lobster_path: PathLike,
        structure_file: str = "POSCAR",
        icohp_file: str = "ICOHPLIST.lobster",
    ):
        """
        Construct a LobsterNeighborsPlotter from a LOBSTER calculation directory.

        Args:
            lobster_path : str or Path
                Path to the LOBSTER calculation directory.
            bond_tol : float
                Tolerance added/subtracted from bond lengths for VESTA SBOND entries.
            file_name : str
                Output VESTA file name.
            structure_file : str, optional
                Structure file used by LOBSTER (default: "POSCAR").
            icohp_file : str, optional
                ICOHP list file produced by LOBSTER (default: "ICOHPLIST.lobster").

        Returns:
            LobsterNeighborsPlotter
        """
        # TODO

    def _get_relevant_bonds(self) -> list:
        bonds = []
        for site_id, site in enumerate(self.lobster_neighbors.structure.sites):
            nbs = self.lobster_neighbors.get_nn_info(structure=self.lobster_neighbors.structure, n=site_id)
            for nb in nbs:
                site_str = f"{site.species_string}{site_id + 1}"
                site_to_str = f"{nb['site'].species_string}{nb['site_index'] + 1}"
                length = round(nb["edge_properties"]["bond_length"], 5)
                if [site_to_str, site_str, length] not in bonds:
                    # TODO advanced: get bond lengths of relevant and not-relevant ICOHPS and
                    # define meaningful distance intervals for vesta
                    # TODO advanced: bond thickness in vesta as a function of ICOHP values?
                    bonds.append([site_str, site_to_str, length])
        return bonds

    def write_vesta(self) -> None:
        """
        Write a VESTA file with bonds defined as per LOBSTER calculation output.

        """
        # partly adapted fr. Janine George:
        # https://github.com/materialsproject/pymatgen/blob/682bfd855dd89264b0ff8d6bee4815f8834cc5ac/src/pymatgen/phonon/thermal_displacements.py#L310
        with open(self.file_name, mode="w", encoding="utf-8") as file:
            file.write("#VESTA_FORMAT_VERSION 3.5.4\n \n \n")
            file.write("CRYSTAL\n\n")
            file.write("TITLE\n")
            file.write("Custom bonds\n\n")
            file.write("GROUP\n")
            file.write("1 1 P 1\n\n")
            file.write("CELLP\n")
            file.write(
                f"{self.lobster_neighbors.structure.lattice.a} "
                f"{self.lobster_neighbors.structure.lattice.b} "
                f"{self.lobster_neighbors.structure.lattice.c} "
                f"{self.lobster_neighbors.structure.lattice.alpha} "
                f"{self.lobster_neighbors.structure.lattice.beta} "
                f"{self.lobster_neighbors.structure.lattice.gamma}\n"
            )
            file.write("  0.000000   0.000000   0.000000   0.000000   0.000000   0.000000\n")  # error on parameters
            file.write("STRUC\n")  # codespell:ignore struc

            for site_idx, site in enumerate(self.lobster_neighbors.structure, start=1):
                file.write(
                    f"{site_idx} {site.species_string} {site.species_string}{site_idx} 1.0000 {site.frac_coords[0]} "
                    f"{site.frac_coords[1]} {site.frac_coords[2]} 1a 1\n"
                )
                file.write(" 0.000000 0.000000 0.000000 0.00\n")

            file.write("  0 0 0 0 0 0 0\n")

            # Assumption: no thermal displacement info available
            for site_idx, site in enumerate(self.lobster_neighbors.structure, start=1):
                file.write(f"{site_idx}  {site.species_string}{site_idx} 0.00000 \n")
                file.write(" 0.000000 0.000000 0.000000 \n")

            # Define bonds
            file.write("SBOND\n")
            file.writelines(
                f"{bond_id}  {bond[0]}  {bond[1]}  {bond[2] - self.bond_tol}  {bond[2] + self.bond_tol} "
                f" 0  1  1  1  1  0.250  2.000 127 127 127\n"
                for bond_id, bond in enumerate(self._bonds, start=1)
            )
            file.write("0 0 0 0\n")

            # Minimal styling
            file.write(
                "STYLE\n"
                "MODEL   2  1  0\n"
                "SURFS   0  1  1\n"
                "FORMS   0  1\n"
                "ATOMS   0  0  1\n"
                "BONDS   1\n"
                "POLYS   1\n"
                "POLYP\n"
                " 204 1  1.000 180 180 180\n"
            )
