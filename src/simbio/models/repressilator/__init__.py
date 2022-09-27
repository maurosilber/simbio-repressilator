from simbio.model import (
    EmptyCompartment,
    EmptyGroup,
    Parameter,
    SingleReaction,
    Species,
)
from simbio.reactions.single import Destruction


class ProportionalCreation(SingleReaction):
    """
    A -> A + B
    """

    A: Species
    B: Species
    rate: Parameter

    @property
    def reactants(self):
        return (self.A,)

    @property
    def products(self):
        return (self.A, self.B)

    @staticmethod
    def reaction_rate(t, A, rate):
        return rate * A


class Hill(SingleReaction):
    A: Species
    B: Species
    rate: Parameter
    p: Parameter
    n: Parameter

    @property
    def reactants(self):
        return (self.A,)

    @property
    def products(self):
        return (self.A, self.B)

    @staticmethod
    def reaction_rate(t, A, rate, p, n):
        return rate / (1 + p**n)


class Gene(EmptyGroup):
    mRNA: Species = 0
    protein: Species = 0
    mRNA_degradation: Parameter
    protein_to_mRNA_degradation_ratio: Parameter  # beta

    remove_mRNA = Destruction(A=mRNA, rate=mRNA_degradation)
    remove_protein = Destruction(A=protein, rate=protein_to_mRNA_degradation_ratio)
    create_protein = ProportionalCreation(
        A=mRNA,
        B=protein,
        rate=protein_to_mRNA_degradation_ratio,
    )


class Repressilator(EmptyCompartment):
    k: Parameter
    beta: Parameter
    rate: Parameter
    p: Parameter
    n: Parameter

    lacl = Gene(mRNA_degradation=k, protein_to_mRNA_degradation_ratio=beta)
    tetR = Gene(mRNA_degradation=k, protein_to_mRNA_degradation_ratio=beta)
    cl = Gene(mRNA_degradation=k, protein_to_mRNA_degradation_ratio=beta)

    create_lacl_mRNA = Hill(A=lacl.protein, B=tetR.mRNA, rate=rate, p=p, n=n)
    create_tetR_mRNA = Hill(A=tetR.protein, B=cl.mRNA, rate=rate, p=p, n=n)
    create_cl_mRNA = Hill(A=cl.protein, B=lacl.mRNA, rate=rate, p=p, n=n)
