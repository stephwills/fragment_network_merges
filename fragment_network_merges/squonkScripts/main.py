#!/usr/bin/env python
import glob
import json
import os
import shutil
import tempfile
from typing import List, Optional

from rdkit import Chem


def main(fragments:Optional[List[str]], proteins:Optional[List[str]],
         outfile:str, fragIdFieldName:Optional[str], fragIdFieldValue:Optional[str], smilesFieldName:Optional[str],
         proteinFieldName:Optional[str], proteinFieldValue:Optional[str], wdir:Optional[str]):
    """

    :param fragments:
    :param proteins:
    :param outfile:
    :param fragIdFieldName:
    :param fragIdFieldValue:
    :param smilesFieldName:
    :param proteinFieldName:
    :param proteinFieldValue:
    :param wdir:
    """
    from fragment_network_merges.merge.config_merge import config_merge
    from fragment_network_merges.filter.config_filter import config_filter

    os.environ["NEO4J_USER"] = os.environ.get("NEO4J_USERNAME", os.environ.get("NEO4J_USER", config_merge.NEO4J_USER))
    os.environ["NEO4J_PASS"] = os.environ.get("NEO4J_PASSWORD", os.environ.get("NEO4J_PASS", config_merge.NEO4J_PASS))
    os.environ["NEO4J_URI"] = os.environ.get("NEO4J_URI", os.environ.get("NEO4J_SERVER", config_merge.NEO4J_URI))

    # if fragment_protein is not None:
    #     assert proteins is None and fragments is None
    #     fragments, proteins = zip(*fragment_protein)
    #     fragments = list(fragments)
    #     proteins = list(proteins)

    nprots = len(proteins)
    nfrags = len(fragments)
    if nprots >1:
        assert nprots == nfrags
    else:
        proteins = proteins*nfrags

    fragmentIds = []

    #1st, recreate old fragalysis stylestructure.
    with tempfile.TemporaryDirectory() as tmpdir:
        if wdir is not None:
            wdir = os.path.abspath(wdir)
            os.makedirs(wdir, exist_ok=True)
            tmpdir = wdir
        cwdir = os.getcwd()
        os.chdir(tmpdir)
        targetName = os.path.basename(fragments[0]).split("-")[0]
        print("targetName", targetName)
        os.makedirs(os.path.join(targetName, "aligned"), exist_ok=True)
        for i, frag in enumerate(fragments):
            frag_noext = os.path.splitext(frag)[0]
            dirname = "_".join(os.path.basename(frag_noext).split("_")[:2])
            fragId = "_".join(dirname.split("-")[1:])
            print(fragId)
            fragmentIds.append(fragId)
            fragwdir = os.path.join(targetName, "aligned", dirname)
            os.makedirs(fragwdir, exist_ok=True)
            fragFname = os.path.join(fragwdir, dirname+".mol")
            frag = os.path.join(cwdir, frag)
            if frag.endswith(".sdf"):
                sup = Chem.SDMolSupplier(frag)
                mol = sup[0]
                Chem.MolToMolFile(mol, fragFname)
            elif frag.endswith(".mol"):
                shutil.copyfile(frag, fragFname)
                mol = Chem.MolFromMolFile(frag)
            else:
                raise ValueError()
            smilesFname = os.path.join(fragwdir, dirname+"_smiles.txt")
            with open(smilesFname, "w") as f:
                f.write(Chem.MolToSmiles(mol, isomericSmiles=False))
            prot = proteins[i]
            apo_desolv_name = os.path.join(fragwdir, dirname+"_apo-desolv.pdb")
            prot = os.path.join(cwdir, prot)
            shutil.copyfile(prot, apo_desolv_name)
        print(f"Input data prepared at {tmpdir}")

        print("Launching queries")
        output_dir_query = os.path.join(tmpdir, "queries_output")
        os.makedirs(output_dir_query, exist_ok=True)
        wdir = os.path.join(tmpdir, "working_dir")
        os.makedirs(wdir, exist_ok=True)

        config_merge.FRAGALYSIS_DATA_DIR = tmpdir
        config_filter.FRAGALYSIS_DATA_DIR = tmpdir
        config_filter.N_CPUS_FILTER_PAIR= int(os.environ.get("N_CPUS_FILTER_PAIR", 1)) #TODO: Tim, how to ask for a given number of cpus?
        config_filter.SCORING_PIPELINE = [ "SuCOSScore"] #Not working with "IfpScore"

        from fragment_network_merges.merge.query import run_query
        print(targetName)
        print(fragmentIds)
        print(config_merge.FRAGALYSIS_DATA_DIR)
        run_query(fragmentIds, targetName, output_dir_query, remove_similar_fragments=False,  working_dir=wdir)
        # import subprocess; subprocess.call(f"cp -r ~/tmp/steph_trial/* {tmpdir}", shell=True)
        print("Queries were executed!", output_dir_query, wdir)
        wdir = os.path.join(tmpdir, "working_dir_filter")
        output_dir_filter = os.path.join(tmpdir, "filter_output")
        os.makedirs(wdir, exist_ok=True)
        os.makedirs(output_dir_filter, exist_ok=True)

        allPairsFname = os.path.join(tmpdir, "working_dir", f"{targetName}_pairs.json")
        if os.path.exists(allPairsFname):
            with open(allPairsFname) as f:
                all_pairs = json.load(f)
        else:
            all_pairs = []
            print("No valid pairs were found!!")

        to_filter = []
        for pair in all_pairs:
            fname = os.path.join(output_dir_query, "_".join(pair)+".json")
            if os.path.isfile(fname):
                to_filter.append(tuple(pair)+(fname,))

        from fragment_network_merges.filter.filter_pipeline import run_filter_pipeline
        for i, (fA,fB, merge_file) in enumerate(to_filter):
            print(f"Filtering {fA,fB}")
            print("FRAGALYSIS_DATA_DIR", config_filter.FRAGALYSIS_DATA_DIR)
            run_filter_pipeline(fA, fB, targetName, merge_file, working_dir=wdir,
                                output_dir=output_dir_filter, sim_search=False)
            print("Next pair!")

        filtered_jsons = glob.glob(f"{output_dir_filter}/**/*_filtered.json", recursive=True)
        print(filtered_jsons)
        mols = []
        for filename in filtered_jsons:
            with open(filename) as f:
                data = json.load(f)
            print(data)
            if not data:
                continue
            for mol_name, info in data.items():
                print(mol_name)
                # Get the basename of the file (without extension)
                dirname, basename = os.path.split(filename)
                mol_dir = os.path.join(dirname, mol_name.replace("_", "-"))
                mol_fname = glob.glob(f"{mol_dir}/*.mol", recursive=False)[0]

                # # Read the .mol file
                mol = Chem.MolFromMolFile(mol_fname)
                #
                # Set the molecule name
                mol.SetProp('_Name', mol_name)

                def replace_even_dashes(s):
                    # Split the string at every "-"
                    out = ""
                    i = 0
                    for c in s:
                        if c == "-":
                            if i % 2 == 1:
                                c = ","
                            else:
                                c = "_"
                            i += 1
                        out += c
                    return out

                frags = replace_even_dashes(info["pair"])
                if fragIdFieldName:
                    mol.SetProp(fragIdFieldName, frags)
                if smilesFieldName:
                    mol.SetProp(smilesFieldName, Chem.MolToSmiles(mol))
                if "PlipIfpScore" in info:
                    mol.SetDoubleProp("PlipIfpScore", info["PlipIfpScore"])
                mol.SetDoubleProp("SuCOSScore", info["SuCOSScore"])
                if proteinFieldValue:
                    _proteinFieldValue = proteinFieldValue if nprots == 1 else frags.split(",")[0]
                    if proteinFieldName:
                        mol.SetProp(proteinFieldName, _proteinFieldValue)  # TODO: Check if this is what we want
                # # Add the molecule to the list
                mols.append(mol)
        print(f"Found {len(mols)} mols. Saved at {outfile}")
        # Write all molecules to the .sdf file
        os.chdir(cwdir)
        assert outfile
        writer = Chem.SDWriter(outfile)
        for mol in mols:
            writer.write(mol)

        # Close the writer
        writer.close()
        print("Done!")


def _test():
    fragments = glob.glob("/home/sanchezg/oxford/data/fragalysis/nsp13/aligned/nsp13-*/nsp13-*.mol")
    proteins = glob.glob("/home/sanchezg/oxford/data/fragalysis/nsp13/aligned/nsp13-*/nsp13-*apo-desolv.pdb")
    #nsp13 "x0176_0B" vs "x0246_0B" gives nice results
    fragments = [x for x in fragments if "x0176_0B" in x or "x0246_0B" in x]
    proteins = [x for x in proteins if "x0176_0B" in x or "x0246_0B" in x]

    """
    Valid examples
 Generating synthons: x0176_0B and x0246_0B
    """
    outfile = "/tmp/fragnent.sdf"
    fragIdField = "ref_mols"
    smilesFieldName = "SMILES"
    proteinFieldName = "ref_pdb"
    proteinFieldValue = "x0176_0B"

    os.environ["KUBECONFIG"] = "$HOME/.kube/_fragment_network_config"
    os.environ["NEO4J_USER"] = "ruben"
    os.environ["NEO4J_PASS"] = "oJ3%xp0d"
    """
kubectl port-forward pods/graph-0 -n graph-b 7474:7474 &
kubectl port-forward pods/graph-0 -n graph-b 7687:7687 &
    """

    main(fragments, proteins, outfile, fragIdField, smilesFieldName, proteinFieldName, proteinFieldValue)

if __name__ == "__main__":
    # _test()
    import sys
    sys.path.append(os.path.dirname(os.path.dirname(__file__)))
    from argParseFromDoc import parse_function_and_call
    parse_function_and_call(main)