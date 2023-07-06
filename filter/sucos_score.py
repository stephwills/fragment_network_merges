# Copyright <2019> <University of Oxford>
# This code is licensed under MIT license (see LICENSE.txt for details)
# Code from https://github.com/MarcMoesser/SuCOS

import argparse
import gzip
import os
import sys

import rdkit
from filter.config_filter import config_filter
from filter.generic_scoring import Score_generic
from joblib import Parallel, delayed
from rdkit import Chem, RDConfig
from rdkit.Chem import AllChem, rdmolfiles, rdShapeHelpers
from rdkit.Chem.FeatMaps import FeatMaps

#################################################
#### Setting up the features to use in FeatureMap
fdef = AllChem.BuildFeatureFactory(
    os.path.join(RDConfig.RDDataDir, "BaseFeatures.fdef")
)
#    keep = ('Donor','Acceptor','NegIonizable','PosIonizable','Aromatic')

fmParams = {}
for k in fdef.GetFeatureFamilies():
    fparams = FeatMaps.FeatMapParams()
    fmParams[k] = fparams

keep = (
    "Donor",
    "Acceptor",
    "NegIonizable",
    "PosIonizable",
    "ZnBinder",
    "Aromatic",
    "Hydrophobe",
    "LumpedHydrophobe",
)

#################################################


def get_FeatureMapScore(small_m, large_m, score_mode=FeatMaps.FeatMapScoreMode.Best):
    featLists = []
    for m in [small_m, large_m]:
        rawFeats = fdef.GetFeaturesForMol(m)
        # filter that list down to only include the ones we're intereted in
        featLists.append([f for f in rawFeats if f.GetFamily() in keep])
    fms = [
        FeatMaps.FeatMap(feats=x, weights=[1] * len(x), params=fmParams)
        for x in featLists
    ]
    fms[0].scoreMode = score_mode
    fm_score = fms[0].ScoreFeats(featLists[1]) / min(
        fms[0].GetNumFeatures(), len(featLists[1])
    )
    return fm_score


def calc_SuCOS(
    ref_file,
    prb_file,
    score_mode=FeatMaps.FeatMapScoreMode.Best,
    write=False,
    return_all=False,
):
    if type(ref_file) == str:
        if os.path.splitext(ref_file)[-1] == ".sdf":
            reflig = Chem.MolFromMolFile(ref_file, sanitize=True)
        elif os.path.splitext(ref_file)[-1] == ".mol2":
            reflig = Chem.MolFromMol2File(ref_file, sanitize=True)
        elif os.path.splitext(ref_file)[-1] == ".mol":
            reflig = rdmolfiles.MolFromMolFile(ref_file, sanitize=True)
    elif type(ref_file) == rdkit.Chem.rdchem.Mol:
        reflig = ref_file

    if type(prb_file) == str:
        if os.path.splitext(prb_file)[-1] == ".sdf":
            prb_mols = Chem.SDMolSupplier(prb_file, sanitize=True)
        elif os.path.splitext(prb_file)[-1] == ".mol":
            prb_mols = [rdmolfiles.MolFromMolFile(prb_file, sanitize=True)]
        elif os.path.splitext(prb_file)[-1] == ".gz":
            tmp = os.path.splitext(prb_file)[0]
            if os.path.splitext(tmp)[-1] == ".sdf":
                inf = gzip.open(prb_file)
                prb_mols = Chem.ForwardSDMolSupplier(inf, sanitize=True)
    elif type(prb_file) == rdkit.Chem.rdchem.Mol:
        prb_mols = [prb_file]

    try:
        reflig
    except NameError:
        raise ValueError("Incorrect file format for ref lig")
    try:
        prb_mols
    except NameError:
        raise ValueError("Incorrect file format for prb lig")

    if write:
        w = Chem.SDWriter("%s_SuCOS_score.sdf" % os.path.splitext(prb_file)[0])
    prb_mols = [x for x in prb_mols if x]

    for prb_mol in prb_mols:
        ##############################################
        ####### Feature map
        ##############################################

        fm_score = get_FeatureMapScore(reflig, prb_mol, score_mode)
        # fm_score = np.clip(fm_score, 0, 1) #Commented out, no longer necessary if using score_mode=Best (~Marc)
        ##############################################

        # tversky_ind = rdShapeHelpers.ShapeTverskyIndex(reflig, prb_mol, 1.0, 0.0)
        # SuCOS_score = 0.5*fm_score + 0.5*tversky_ind

        protrude_dist = rdShapeHelpers.ShapeProtrudeDist(
            reflig, prb_mol, allowReordering=False
        )
        # protrude_dist = np.clip(protrude_dist, 0, 1)
        SuCOS_score = 0.5 * fm_score + 0.5 * (1 - protrude_dist)

        # For high throughput consider commenting out print statements to deconvolute terminal (~Marc)
        # print("********************************")
        # print("SuCOS score:\t%f" % SuCOS_score)
        # print("Chem features:\t%f" % fm_score)
        # print("ShapeProtrudeDist:\t%f" % protrude_dist)
        # print("********************************")

        prb_mol.SetProp("SuCOS_score", str(SuCOS_score))
        prb_mol.SetProp("Volume_score", str(1 - protrude_dist))
        prb_mol.SetProp("Feature_score", str(fm_score))
        if write:
            w.write(prb_mol)
    if return_all:
        return SuCOS_score, fm_score, (1 - protrude_dist)
    else:
        return SuCOS_score


class SuCOSScore(Score_generic):
    def __init__(
        self,
        smis=None,
        synthons=None,
        fragmentA=None,
        fragmentB=None,
        proteinA=None,
        proteinB=None,
        merge=None,
        mols=None,
        names=None,
        work_pair_dir=None,
        out_pair_dir=None,
        mol_files=None,
        holo_files=None,
        apo_files=None,
    ):
        super().__init__(
            smis,
            synthons,
            fragmentA,
            fragmentB,
            proteinA,
            proteinB,
            merge,
            mols,
            names,
            work_pair_dir,
            out_pair_dir,
            mol_files,
            holo_files,
            apo_files,
        )
        self.scores = None

    def score_mol(
        self, mol_file: str, fragmentA_file: str, fragmentB_file: str
    ) -> float:
        """
        Calculates the average SuCOS score

        :param mol_file: the file of the merge to filter
        :type mol_file: str
        :param fragmentA_file: the fragment A file
        :type fragmentA_file: str
        :param fragmentB_file: the fragment B file
        :type fragmentB_file: str

        :return: avg SuCOS score
        :rtype: float
        """
        scoreA = calc_SuCOS(fragmentA_file, mol_file, FeatMaps.FeatMapScoreMode.Best)
        scoreB = calc_SuCOS(fragmentB_file, mol_file, FeatMaps.FeatMapScoreMode.Best)
        avg = (scoreA * scoreB) ** 0.5
        return avg

    def score_all(self, cpus: int = config_filter.N_CPUS_FILTER_PAIR) -> list:
        self.scores = Parallel(n_jobs=cpus, backend="multiprocessing")(
            delayed(self.score_mol)(mol_file, self.fragmentA, self.fragmentB)
            for mol_file in self.mol_files
        )

        return self.scores

    def filter_mol(
        self,
        mol_file: str,
        fragmentA_file: str,
        fragmentB_file: str,
        score_threshold: float,
    ) -> float:
        """
        Filters using the average SuCOS score
        """
        scoreA = calc_SuCOS(fragmentA_file, mol_file, FeatMaps.FeatMapScoreMode.Best)
        scoreB = calc_SuCOS(fragmentB_file, mol_file, FeatMaps.FeatMapScoreMode.Best)
        avg = (scoreA * scoreB) ** 0.5

        if avg >= score_threshold:
            result = True
        else:
            result = False
        return result

    def filter_all(self, cpus: int = config_filter.N_CPUS_FILTER_PAIR, *args):
        """
        Runs the SuCOS score filter on all molecules in parallel.
        """
        results = Parallel(n_jobs=cpus, backend="multiprocessing")(
            delayed(self.filter_mol)(mol_file, self.fragmentA, self.fragmentB, *args)
            for mol_file in self.mol_files
        )

        return results, self.mols


def parse_args(args):
    """
    Parse command line arguments
    """
    parser = argparse.ArgumentParser(
        epilog="""
    python filter/sucos_score.py --input_file data/toFilter.sdf --output_file results.sdf --score_threshold 0.6
    """
    )
    # command line args definitions
    parser.add_argument(
        "-i",
        "--input_file",
        required=True,
        help="input sdf file containing molecules to filter",
    )
    parser.add_argument(
        "-o",
        "--output_file",
        required=True,
        help="output sdf file to write filtered SMILES",
    )
    parser.add_argument(
        "-a", "--fragmentA_file", required=True, help="fragment A mol file"
    )
    parser.add_argument(
        "-b", "--fragmentB_file", required=True, help="fragment B mol file"
    )
    parser.add_argument(
        "-c",
        "--n_cpus",
        required=False,
        help="number of CPUs",
        type=int,
        default=os.cpu_count(),
    )
    parser.add_argument(
        "-t",
        "--score_threshold",
        required=False,
        help="avg SuCOS score to pass threshold",
        type=float,
        default=0.5,
    )
    return parser.parse_args(args)


def main():
    from filter.generic_squonk import Squonk_generic

    args = parse_args(sys.argv[1:])
    job = Squonk_generic(
        "SuCOSScore",
        args,
        args.input_file,
        args.output_file,
        args.fragmentA_file,
        args.fragmentB_file,
        score=True,
    )
    job.execute_job(args.n_cpus, args.score_threshold)


if __name__ == "__main__":
    main()
