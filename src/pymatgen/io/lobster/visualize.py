from __future__ import annotations

from pymatgen.analysis.local_env import NearNeighborsVisualizer
from pymatgen.io.lobster.lobsterenv import LobsterNeighbors


class LobsterNeighborsVisualizer(NearNeighborsVisualizer):
    """
    Class to create VESTA files with bonds defined from LOBSTER bonding analysis.
    """

    def __init__(
        self,
        near_neighbors: LobsterNeighbors,
        bond_tol: float = 0.1,
    ) -> None:
        """
        Args:
            near_neighbors : pymatgen.io.lobster.LobsterNeighbors
                LobsterNeighbors object used to determine bonding.
            bond_tol : float
                Tolerance added/subtracted from actual bond lengths for VESTA SBOND entries.
        """
        super().__init__(near_neighbors=near_neighbors, structure=near_neighbors.structure, bond_tol=bond_tol)

    @classmethod
    def from_lobster_calc(
        cls,
        bond_tol: float = 0.1,
        *lobster_neighbor_args,
        **lobster_neighbor_kwargs,
    ):
        """
        Construct a LobsterNeighborsVisualizer object from a structure and LOBSTER calculation output.

        Args:
             bond_tol : float
                Tolerance added/subtracted from actual bond lengths for VESTA SBOND entries.
            *args, **kwargs
                args and kwargs supported by LobsterNeighbors init
        Returns:
            LobsterNeighborsVisualizer
        """
        lobster_neighbors = LobsterNeighbors(*lobster_neighbor_args, **lobster_neighbor_kwargs)
        return cls(near_neighbors=lobster_neighbors, bond_tol=bond_tol)
