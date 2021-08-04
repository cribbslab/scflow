'''
scflow.py - Single Cell analysis workflows for cribbslab
========================================================

Our single-cell workflows are organised into three different sections.
* main: contains the upstream pipelines for mapping, velocyto, QC and filtering
* seurat: contains Seurat specific workflows
* scanpy: contains scanpy specific workflows


For this message and a list of available keywords type::

    scflow --help


To see the available pipelines within each section type::

    scflow <section>

To run a specific pipeline/workflow type the following::

    scflow <section> <workflow> [workflow options] [workflow arguments]

To get help for a specify workflow, type::

    scflow <section> <workflow> --help
'''

import os
import sys
import re
import glob
import imp
import scpipelines


def printListInColumns(l, ncolumns):
    '''output list *l* in *ncolumns*.'''
    ll = len(l)

    if ll == 0:
        return

    max_width = max([len(x) for x in l]) + 3
    n = ll // ncolumns
    if ll % 3 != 0:
        n += 1

    # build columns
    columns = [l[x * n:x * n + n] for x in range(ncolumns)]

    # add empty fields for missing columns in last row
    for x in range(ncolumns - (len(l) % ncolumns)):
        columns[-(x + 1)].append('')

    # convert to rows
    rows = list(zip(*columns))

    # build pattern for a row
    p = '%-' + str(max_width) + 's'
    pattern = ' '.join([p for x in range(ncolumns)])

    # put it all together
    return '\n'.join([pattern % row for row in rows])


def main(argv=None):

    argv = sys.argv

    # paths to look for pipelines:
    print(scpipelines.__file__)
    path = os.path.abspath(os.path.dirname(scpipelines.__file__))
    scanpy = path + "/scanpy/"
    seurat = path + "/seurat/"

    paths = [path, scanpy, seurat]

    if len(argv) == 1 or argv[1] == "--help" or argv[1] == "-h":

        print((globals()["__doc__"]))
        print("The list of available sections are:\n")
        print("main  scanpy  seurat\n")
        return

    elif argv[1] == "main":
        print((globals()["__doc__"]))

        pipelines = []
        pipelines.extend(glob.glob(os.path.join(path, "pipeline_*.py")))
        print("The list of available single cell pipelines is:\n")
        print("{}\n".format(
            printListInColumns(
                sorted([os.path.basename(x)[len("pipeline_"):-len(".py")] for x in pipelines]),
                3)))
    elif argv[1] == "seurat":
        print((globals()["__doc__"]))

        pipelines = []
        pipelines.extend(glob.glob(os.path.join(seurat, "pipeline_*.py")))
        print("The list of available single cell seurat pipelines are:\n")
        print("{}\n".format(
            printListInColumns(
                sorted([os.path.basename(x)[len("pipeline_"):-len(".py")] for x in pipelines]),
                5)))

    elif argv[1] == "scanpy":
        print((globals()["__doc__"]))

        pipelines = []
        pipelines.extend(glob.glob(os.path.join(scanpy, "pipeline_*.py")))
        print("The list of available single cell pipelines is:\n")
        print("{}\n".format(
            printListInColumns(
                sorted([os.path.basename(x)[len("pipeline_"):-len(".py")] for x in pipelines]),
                3)))

    else:
        print("please select an appropriate workflow section: main, scanpy or seurat")

    try:
        command = argv[2]
        if command == "-h" or command == "--help":
            return
        pipeline = "pipeline_{}".format(command)
    except Exception:
        print("No pipeline has been selected under the %s section" % (argv[1]))
        print("Problems?", sys.exc_info())
        return

    # remove 'scflow' from sys.argv
    del sys.argv[0]
    del sys.argv[1]

    (file, pathname, description) = imp.find_module(pipeline, paths)

    module = imp.load_module(pipeline, file, pathname, description)

    module.main(sys.argv)


if __name__ == "__main__":
    sys.exit(main())
