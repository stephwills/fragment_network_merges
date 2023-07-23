FROM informaticsmatters/squonk2-fragmenstein-base:stable

ENV HOME=/code
WORKDIR ${HOME}

RUN conda install -y -c conda-forge pymol-open-source plip openbabel neo4j-python-driver=4

RUN pip install pebble==4.6.3 argParseFromDoc==0.1


COPY merge ./merge/
COPY filter ./filter/
COPY squonkScripts ./squonkScripts/
COPY utils ./utils/
COPY setup.py ./
COPY README.md ./
# COPY data ./

RUN pwd
RUN ls
RUN pip install .

ENV N_CPUS_FILTER_PAIR=1
ENV USE_FRAGMENSTEIN_WICTOR=1
# ENV NEO4J_URI="bolt://graph.graph-stfc:7687"
