"""
Used for filtering merges according to calculated descriptors.
"""

from joblib import Parallel, delayed
from oddt import toolkit, fingerprints

from filter.config_filter import config_filter
from filter.generic_scoring import Score_generic


class IfpScore(Score_generic):

    def __init__(self, smis: list, synthons, fragmentA, fragmentB, proteinA, proteinB, merge, mols, names, mol_files,
                 holo_files, apo_files):
        super().__init__(smis, synthons, fragmentA, fragmentB, proteinA, proteinB, merge, mols, names, mol_files,
                         holo_files, apo_files)
        self.scores = None

    def _get_protein(self, protein):
        """
        Function loads the protein from the pdb file.

        :param protein: the protein file
        :type protein: pdb file

        :return: protein
        :rtype: ODDT protein
        """
        protein = next(toolkit.readfile('pdb', protein))
        protein.protein = True
        return protein

    def _get_mol(self, mol):
        """
        Function loads the molecule from the mol file.

        :param mol: mol to read (fragment or placed merge)
        :type mol: mol file

        :return: molecule
        :rtype: ODDT molecule
        """
        return next(toolkit.readfile('mol', mol))

    def make_fp(self, protein, mol):
        """
        Function creates interaction fingerprint between the molecule and protein.

        :param protein: protein loaded from file
        :type protein: ODDT protein
        :param mol: mol loaded from file
        :type mol: ODDT mol

        :return: fingerprint
        :rtype: numpy array
        """
        fp = fingerprints.InteractionFingerprint(mol, protein)
        return fp

    def calc_bond_percentage(self, fragment_fp, merge_fp):
        """
        Alternative to Tversky (which requires the fp to be binary).
        Calculate the percentage of the bonds made by the fragment that are preserved by
        the merge. Only checks bits in which the fragment makes an interaction but still
        accounts for the fact that this is a count vector.

        :param fragment_fp: interaction fingerprint of fragment
        :type fp: numpy array
        :param fp: interaction fingerprint of merge
        :type fp: numpy array

        :return: proportion of interactions maintained by the merge
        :rtype: float
        """
        frag_count = 0
        merge_count = 0

        for i, j in zip(fragment_fp, merge_fp):
            if i > 0:  # only check bits in which the fragment makes an interaction
                frag_count += i
                if i >= j:
                    merge_count += j
                else:
                    merge_count += i  # avoid percentage >100%

        perc = merge_count / frag_count

        return perc

    def ifp_score_avg(self, merge, fragmentA, fragmentB, prot):
        """
        Function to calculate the average percentage bonds preserved between the merge and the fragment
        fingerprints. Calculates for each fragment and returns the average.

        :param merge: the file of the merge to filter
        :type merge: mol file
        :param fragmentA: the fragment A file
        :type fragmentA: mol file
        :param fragmentB: the fragment B file
        :type fragmentB: mol file
        :param prot: protein file
        :type prot: pdb file

        :return: average percentage preserved
        :rtype: float
        """
        # load the molecules
        protein = self._get_protein(prot)
        merge_mol = self._get_mol(merge)
        fA_mol = self._get_mol(fragmentA)
        fB_mol = self._get_mol(fragmentB)
        # create all the fingerprints
        merge_fp, fA_fp, fB_fp = self.make_fp(protein, merge_mol), self.make_fp(protein, fA_mol), self.make_fp(protein, fB_mol)

        perc_A = self.calc_bond_percentage(fA_fp, merge_fp)
        perc_B = self.calc_bond_percentage(fB_fp, merge_fp)
        mean = (perc_A + perc_B) / 2

        return mean

    def ifp_score_total(self, merge, fragmentA, fragmentB, prot):
        """
        Function to calculate the percentage bonds preserved between the merge and both fragment
        fingerprints.

        :param merge: the file of the merge to filter
        :type merge: mol file
        :param fragmentA: the fragment A file
        :type fragmentA: mol file
        :param fragmentB: the fragment B file
        :type fragmentB: mol file
        :param prot: protein file
        :type prot: pdb file

        :return: percentage total bonds preserved
        :rtype: float
        """
        # load the molecules
        protein = self._get_protein(prot)
        merge_mol = self._get_mol(merge)
        fA_mol = self._get_mol(fragmentA)
        fB_mol = self._get_mol(fragmentB)
        # create all the fingerprints
        merge_fp, fA_fp, fB_fp = self.make_fp(protein, merge_mol), self.make_fp(protein, fA_mol), self.make_fp(protein, fB_mol)

        frag_count = 0
        merge_count = 0

        comb_frag_fp = [a + b for a, b in zip(fA_fp, fB_fp)]
        perc = self.calc_bond_percentage(comb_frag_fp, merge_fp)

        return perc

    def score_mol(self, merge, fragmentA, fragmentB, prot, type_calculation='avg'):
        """
        Calculates the percentage of interactions made by the original fragments that are
        maintained by the merge. Can either calculate as the average of the score for each
        fragment or as the percentage of all interactions made.
        """
        if type_calculation == 'avg':
            score = self.ifp_score_avg(merge, fragmentA, fragmentB, prot)
        elif type_calculation == 'total':
            score = self.ifp_score_total(merge, fragmentA, fragmentB, prot)

        return score

    def score_all(self, cpus: int = config_filter.N_CPUS_FILTER_PAIR):

        self.scores = Parallel(n_jobs=cpus, backend='multiprocessing') \
            (delayed(self.score_mol)(mol_file, self.fragmentA, self.fragmentB, apo_file) for mol_file, apo_file in
             zip(self.mol_files, self.apo_files))

        return self.scores
