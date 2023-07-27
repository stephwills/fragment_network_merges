import os
import json
import dask.bag as db


def combine_search_jsons(path_to_jsons):
  def combine_two_dict(d1, d2, ):
    new_d = {}
    assert len(d1) == len(d2), "Error, both dicts should be equaly-sized"
    for query_smi, matches_list1 in d1.items():
      matches_list2 = d2[query_smi]
      assert len(matches_list1) == len(matches_list2), "Error, both list should be equaly-sized"
      final_matches = []
      for match1, match2 in zip(matches_list1, matches_list2):
        similarity1 = match1[0]
        similarity2 = match2[0]
        if similarity1 > similarity2:
          final_matches.append(match1)
        else:
          final_matches.append(match2)
      new_d[query_smi] = final_matches
    return new_d

  def load_json(fname):
    with open(fname) as f:
      return json.load(f)

  filenames = [] #[os.path.join(path_to_jsons, fname) for fname in os.listdir(path_to_jsons) if fname.endswith(".json") ]
  if not isinstance(path_to_jsons, (list, tuple)):
    path_to_jsons = [path_to_jsons]
  for name in path_to_jsons:
    if os.path.isdir(name):
      filenames += [os.path.join(name, fname) for fname in os.listdir(name) if fname.endswith(".json") ]
    elif name.endswith(".json"):
      filenames.append(name)

  if len(filenames) >1:
    bag = db.from_sequence(filenames[1:]).map(load_json).fold(combine_two_dict, initial=load_json(filenames[0]))
    final_search = bag.compute()
  else:
    final_search = load_json(filenames[0])
  return final_search

