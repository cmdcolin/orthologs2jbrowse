"""Ortholog calls as the tables a JBrowse 2 synteny track loads."""
import importlib
import sys

__version__ = "0.1.0"

COMMANDS = {
    "orthofinder": "OrthoFinder Orthogroups.tsv to an N-genome .blocks table",
    "odp": "odp .rbh table to a .blocks table and one BED per species",
    "mcscanx": "MCScanX .collinearity to .anchors, .anchors.simple, PAF or .blocks",
    "compara": "Ensembl Compara homology TSV to pairwise .blocks with its measurements",
    "symbols": "GFF3 files joined on gene symbol to a .blocks table",
    "kaks": "dN and dS per gene pair, to color a table by selection (needs biopython)",
}


def usage():
    width = max(map(len, COMMANDS))
    return "usage: orthologs2jbrowse <command> [args]\n\n" + "\n".join(
        f"  {name.ljust(width)}  {text}" for name, text in COMMANDS.items()
    ) + "\n\nRun orthologs2jbrowse <command> -h for a command's arguments.\n"


def main():
    args = sys.argv[1:]
    if args[:1] in (["--version"], ["-V"]):
        print(__version__)
    elif not args or args[0] not in COMMANDS:
        sys.stderr.write(usage())
        sys.exit(0 if args[:1] in (["-h"], ["--help"]) else 2)
    else:
        sys.argv = [f"orthologs2jbrowse {args[0]}", *args[1:]]
        importlib.import_module(f"orthologs2jbrowse.{args[0]}").main()
