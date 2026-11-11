import contextlib
import io
import os
import subprocess
import sys
import tempfile
import unittest

from orthologs2jbrowse import compara as cmp_mod
from orthologs2jbrowse import mcscanx as mcx
from orthologs2jbrowse import orthofinder as og
from orthologs2jbrowse import symbols as sy


class Converters(unittest.TestCase):
    def setUp(self):
        self.argv = sys.argv

    def tearDown(self):
        sys.argv = self.argv

    def check(self, name, got, want):
        with self.subTest(name):
            self.assertEqual(got, want)

    def test_mcscanx_functions(self):
        check = self.check
        # `mcscanx`: the MCScan adapters throw on a gene id missing from the
        # BED and silently mis-draw a block whose columns are the wrong way round, so
        # what is pinned is the genome split and the column normalization.
        mcx_genes = {
            "a1": ("at1", 100, 200), "a2": ("at1", 300, 400), "a3": ("at1", 500, 600),
            "p1": ("pp3", 1000, 1100), "p2": ("pp3", 2000, 2100),
            "p3": ("pp3", 3000, 3100), "h1": ("at2", 10, 20),
        }
        mcx_tags = ["at", "pp"]
        converted, counts = mcx.convert_blocks(
            [("-", [("a1", "p3", 30), ("a2", "p2", 20), ("a3", "p1", 10)]),
             # MCScanX writes the same file's blocks in either column order, and reports
             # self-synteny alongside the cross-genome ones
             ("+", [("p1", "a1", 5), ("p2", "a2", 5)]),
             ("+", [("a1", "h1", 5)])],
            mcx_genes, mcx_tags)
        check("convert_blocks keeps the cross-genome blocks only", counts["kept"], 2)
        check("convert_blocks skips a same-genome block", counts["skipped"], 1)
        check("convert_blocks puts the first --species in column 1 either way round",
              [p[0] for p in converted[1][0]], ["a1", "a2"])
        check("convert_blocks keeps the block's orientation off the header",
              [c[1][5] for c in converted], ["-", "+"])
        check("convert_blocks names each side's first and last gene by coordinate",
              converted[0][1][:4], ("a1", "a3", "p1", "p3"))
        check("convert_blocks scores the simple row with the anchor count",
              converted[0][1][4], "3")
        check("convert_blocks drops a pair on an unlisted genome",
              mcx.convert_blocks([("+", [("a1", "zz9", 5)])], mcx_genes, mcx_tags)[1],
              {"kept": 0, "skipped": 1, "unknown": 1})
        # MultiGenomePAFAdapter finds each side's assembly by the PanSN prefix, so the
        # tag has to come off the refName and the genome name has to go on
        check("paf_rows spans each side's genes under a PanSN name",
              mcx.paf_rows(converted[:1], mcx_genes, mcx_tags, ["grape", "peach"],
                           lambda chrom, tag, name: chrom[len(tag):],
                           {"grape": {"1": 5000}, "peach": {"3": 9000}}),
              [["grape#1#1", 5000, 99, 600, "-", "peach#1#3", 9000, 999, 3100, 3,
                2101, 255, "genes:i:3"]])
        # one --species is MCScanX's whole-genome-duplication use, where the blocks the
        # pairwise case discards are the entire result
        self_out, self_counts = mcx.convert_blocks(
            [("+", [("a1", "h1", 5)]), ("+", [("a1", "p1", 5)])], mcx_genes, ["at"])
        check("convert_blocks keeps the same-genome blocks for one --species",
              [p[:2] for p in self_out[0][0]], [("a1", "h1")])
        check("convert_blocks skips the cross-genome block for one --species",
              (self_counts["kept"], self_counts["unknown"]), (1, 1))
        # three genomes go to a .blocks table instead, anchored on the first --species.
        # A cell holds one id and MCScanBlocksAdapter joins the non-reference columns
        # through column 0, so what is pinned is the one-ortholog-per-cell reduction and
        # that a pair missing the reference is dropped rather than written somewhere.
        mcx_genes["c1"] = ("tc5", 10, 20)
        mcx_genes["c2"] = ("tc5", 30, 40)
        rows, counts = mcx.build_table(
            [("+", [("a1", "p1", 10), ("a2", "p2", 10), ("a1", "c1", 10)]),
             ("+", [("a1", "p3", 99)]),
             ("+", [("p1", "c2", 10)])],
            mcx_genes, ["at", "pp", "tc"])
        check("build_table writes a row per reference gene, columns in --species order",
              rows, [["a1", "p3", "c1"], ["a2", "p2", "."]])
        check("build_table drops a pair between two non-reference genomes",
              counts["indirect"], 1)
        # an e_value, not a bit score; e_value=0 is MCScanX's "below what it prints"
        check("score_from_evalue is -log10", mcx.score_from_evalue("1e-77"), 77)
        check("score_from_evalue caps a zero e_value", mcx.score_from_evalue("0"), 1000)
        # BED column 1 has to match the assembly, which never carries MCScanX's tag
        check("bed_rows strips the species tag and converts to a 0-based start",
              mcx.bed_rows(mcx_genes, "pp", {"p1": "-"}, False, "Pp0")[0],
              ("Pp03", 999, 1100, "p1", 0, "-"))
        check("bed_rows defaults an unknown strand to +",
              mcx.bed_rows(mcx_genes, "at", {}, True, "")[0],
              ("at1", 99, 200, "a1", 0, "+"))
        # a refName the assembly does not have draws nothing rather than erroring, so
        # the tag-stripping is checked against the .fai before the BED is written
        fai = os.path.join(tempfile.mkdtemp(), "grape.fa.fai")
        with open(fai, "w") as fh:
            fh.write("1\t100\t0\t60\t61\n2\t100\t0\t60\t61\n")
        check("unmatched_refnames names the BED refNames the assembly lacks",
              mcx.unmatched_refnames(mcx.bed_rows(mcx_genes, "at", {}, True, ""), fai)[0],
              ["at1", "at2"])
        check("unmatched_refnames passes the stripped names",
              mcx.unmatched_refnames(mcx.bed_rows(mcx_genes, "at", {}, False, ""), fai)[0],
              [])

    def test_orthofinder_and_symbols(self):
        check = self.check
        # `orthofinder`: an OrthoFinder cell holds every gene of that genome
        # in the orthogroup, so the duplicated case is the whole design. Reducing it to
        # one gene draws a confident link the data does not support, and in a polyploid
        # it hides the duplication that is the reason for looking.
        dup = ["Os01g1", "Zm01a, Zm01b", "Sb01g"]
        check("expand emits a row per copy, every direct link kept",
              og.orthogroup_rows(dup, "expand", 4, [None] * 3),
              [["Os01g1", "Zm01a", "Sb01g"], ["Os01g1", "Zm01b", "Sb01g"]])
        check("single empties an ambiguous cell rather than choosing",
              og.orthogroup_rows(dup, "single", 4, [None] * 3),
              [["Os01g1", ".", "Sb01g"]])
        check("first takes the gene OrthoFinder happened to list first",
              og.orthogroup_rows(dup, "first", 4, [None] * 3),
              [["Os01g1", "Zm01a", "Sb01g"]])
        # expansion is index-paired, not a product: two duplicated columns cost two rows
        check("expand pairs copies by index instead of multiplying them out",
              len(og.orthogroup_rows(["A1, A2", "B1, B2", "C1"], "expand", 4, [None] * 3)), 2)
        check("a cell past --max-copies is a gene family and contributes nothing",
              og.orthogroup_rows(["Os1", "Z1, Z2, Z3", "Sb1"], "expand", 2, [None] * 3),
              [["Os1", ".", "Sb1"]])
        # ...and is COUNTED, since a threshold below the ploidy in the set empties the
        # cells the table exists to show and reads as a genome with fewer orthologs. The
        # BED filter empties a cell too, and the placement share already reports that one.
        fam = {"orthogroups": 0, "expanded": 0, "families": 0}
        og.orthogroup_rows(["Os1", "Z1, Z2, Z3", "Sb1"], "expand", 2, [None] * 3, fam)
        og.orthogroup_rows(["Os1", "Z1", "Sb1"], "expand", 2, [None] * 3, fam)
        og.orthogroup_rows(["Os1", "Z1, Z2", "Sb1"], "expand", 2, [None, {}, None], fam)
        check("a cell --max-copies emptied is counted, one the BED emptied is not",
              fam["families"], 1)
        # a row with one gene names no link, and would draw nothing while inflating the
        # table
        check("an orthogroup present in one genome only is dropped",
              og.orthogroup_rows(["Os1", "", ""], "expand", 4, [None] * 3), [])
        # OrthoFinder ids come from the protein FASTA headers; ids the BED cannot
        # resolve draw nothing, so they are dropped here where the count is reported
        check("a gene the BED does not have is not counted as an answer",
              og.orthogroup_rows(dup, "expand", 4, [None, {"Zm01b": "Zm01b"}, None]),
              [["Os01g1", "Zm01b", "Sb01g"]])
        # OrthoFinder writes `E(spl)m2-BFM` as `E_spl_m2-BFM`; the table has to carry
        # the BED's spelling, which is what MCScanBlocksAdapter looks the gene up by
        spl = os.path.join(tempfile.mkdtemp(), "spl.bed")
        with open(spl, "w") as fh:
            fh.write("3R\t100\t200\tE(spl)m2-BFM\n3R\t300\t400\tmpka1;1\n")
        check("an id OrthoFinder rewrote resolves to the BED's own spelling",
              og.orthogroup_rows(["E_spl_m2-BFM", "mpka1_1"], "expand", 4,
                                 [og.read_bed_names(spl)] * 2),
              [["E(spl)m2-BFM", "mpka1;1"]])
        # MCScanX_h picks the collinear copy from gene order, so it gets every pairing,
        # not the index-paired rows that already chose one
        check("--homology writes the product of two duplicated cells",
              og.homology_pairs(["A1, A2", "B1, B2", "C1"], [(0, 1), (1, 2)], 4, [None] * 3),
              {(0, 1): [("A1", "B1"), ("A1", "B2"), ("A2", "B1"), ("A2", "B2")],
               (1, 2): [("B1", "C1"), ("B2", "C1")]})
        check("--homology drops a gene family and an unplaced gene as the table does",
              og.homology_pairs(["A1", "B1, B2, B3", "C1, C2"], [(0, 1), (0, 2)], 2,
                                [None, None, {"C2": "C2"}]),
              {(0, 1): [], (0, 2): [("A1", "C2")]})
        # `symbols`: the same three modes over a symbol join, where a column's
        # copies are the genes sharing one symbol. Expanding them is what draws a
        # duplication as a ribbon per copy; taking the first hid the second copy of the
        # ANCHOR too, which left that gene joined to nothing in either direction.
        # orthofinder and symbols each carry pick/--max-copies/index-pairing; this
        # holds the two copies to the same answer.
        grid = [
            [["a1"], ["b1", "b2"], ["c1"]],
            [["a1", "a2"], ["b1", "b2", "b3"], []],
            [["a1", "a2", "a3"], ["b1"], ["c1", "c2"]],
            [[], ["b1"], []],
            [["a1"], [], ["c1", "c2", "c3", "c4", "c5"]],
        ]
        check("symbol_rows and orthogroup_rows agree on every pick and --max-copies",
              [sy.symbol_rows(cells, pick, mc)
               for cells in grid for pick in ("expand", "first", "single")
               for mc in (0, 1, 2, 4)],
              [og.orthogroup_rows([", ".join(c) for c in cells], pick, mc, [None] * 3)
               for cells in grid for pick in ("expand", "first", "single")
               for mc in (0, 1, 2, 4)])
        sdup = [["h1"], ["c1a", "c1b"], ["g1"]]
        check("expand emits a row per copy of the symbol",
              sy.symbol_rows(sdup, "expand", 4),
              [["h1", "c1a", "g1"], ["h1", "c1b", "g1"]])
        check("a repeated anchor symbol keeps both of its copies",
              sy.symbol_rows([["h1", "h1b"], ["c1"], []], "expand", 4),
              [["h1", "c1", "."], ["h1b", "c1", "."]])
        check("single empties a multi-copy cell rather than choosing",
              sy.symbol_rows(sdup, "single", 4), [["h1", ".", "g1"]])
        # copies pair by index, so each genome lists its own in the gene's reading
        # direction: a chromosome deposited the other way round would pair the two
        # halves of a split gene first with last
        halves = lambda strand: [
            {"ref": "c", "start": 10, "strand": strand, "label": "w"},
            {"ref": "c", "start": 50, "strand": strand, "label": "w-2"},
        ]
        check("a minus-strand symbol's copies read from the high coordinate down",
              [[g["label"] for g in sy.reading_order(halves(s))] for s in "+-"],
              [["w", "w-2"], ["w-2", "w"]])
        # no direction runs from X to Y, so copies on two sequences keep file order
        par = [
            {"ref": "X", "start": 10, "strand": "-", "label": "p"},
            {"ref": "Y", "start": 10, "strand": "-", "label": "p-2"},
        ]
        check("copies on different sequences keep the order the annotation lists",
              [g["label"] for g in sy.reading_order(par)], ["p", "p-2"])
        check("first takes the first copy in file order",
              sy.symbol_rows(sdup, "first", 4), [["h1", "c1a", "g1"]])
        check("a symbol past --max-copies is a gene family and empties its cell",
              sy.symbol_rows([["h1"], ["c1", "c2", "c3"], ["g1"]], "expand", 2),
              [["h1", ".", "g1"]])
        check("a symbol only one genome names is dropped",
              sy.symbol_rows([["h1"], [], []], "expand", 4), [])
        check("NCBI's CnHmorfK homolog spelling reads as the human CmorfK",
              [sy.human_symbol(s) for s in
               ("C1H1orf35", "C2AH2orf40", "c16h15orf48", "CXHXorf49B", "CYHYorf15a")],
              ["C1orf35", "C2orf40", "C15orf48", "CXorf49B", "CYorf15A"])
        check("a symbol outside that spelling is left as written",
              [sy.human_symbol(s) for s in ("C1orf35", "CHORF1", "sorF", "EBAG9UORF", None)],
              ["C1orf35", "CHORF1", "sorF", "EBAG9UORF", None])

        # an orthogroup is a set, not a reference-anchored row: unlike a jcvi .blocks
        # table no column anchors the others, so a pair that skips column 0 is a row
        check("a row survives without column 0, which anchors nothing",
              og.orthogroup_rows(["", "B1", "C1"], "expand", 4, [None] * 3),
              [[".", "B1", "C1"]])
        # the header row is the column order, which is what blockAssemblies must be
        check("column_names drops the proteome file's extensions",
              og.column_names(["rice.pep.fa", "maize.pep"], {}), ["rice", "maize"])
        check("column_names applies --assembly to the raw or stripped name",
              og.column_names(["rice.pep.fa", "Zm-B73.pep"], {"Zm-B73": "maize"}),
              ["rice", "maize"])
        # The four required BED columns are a valid BED, and the adapter reads one, so
        # the id check has to. Only the last field carries the line ending, which put it
        # on the name and resolved every id in that column to nothing.
        bd = tempfile.mkdtemp()
        minimal = os.path.join(bd, "minimal.bed")
        with open(minimal, "w") as fh:
            fh.write("# a comment\ntrack name=genes\nchr1\t100\t200\tOs01g1\n"
                     "chr1\t300\t400\tOs01g2\t0\t+\n")
        check("read_bed_names reads a four-column BED, and skips its header lines",
              og.read_bed_names(minimal), {"Os01g1": "Os01g1", "Os01g2": "Os01g2"})
        # The share is what makes "near zero" readable; a bare count is not, and an
        # unchecked column must not print as one that resolved.
        placed, dead = og.report_columns(
            ["rice", "maize", "sorghum"], {"rice": "rice.bed", "maize": "maize.bed"},
            [{"Os1", "Os2"}, {"Zm1", "Zm2"}, {"Sb1"}], [{"Os1"}, set(), set()])
        check("a column resolving none of its ids is named as an id mismatch",
              dead, ["maize"])
        check("the per-column report gives the share, and marks an unchecked column",
              placed.splitlines(),
              ["  rice: 1/2 ids (50%) placed by rice.bed",
               "  maize: 0/2 ids (0%) placed by maize.bed",
               "  sorghum: 1 ids, unchecked (pass --bed sorghum=FILE)"])
        # distinct ids, not output cells: `expand` repeats the single-copy gene once per
        # copy beside it, so counting cells reports more rice genes than the table holds
        _, seen, resolved, _, _ = og.build_rows(
            ["OG1\tOs01g1\tZm01a, Zm01b\n"], ["rice", "maize"], "expand", [None] * 2, 4)
        check("the placed count is distinct ids, not the rows expansion wrote",
              [len(seen[0]), len(resolved[0])], [1, 1])
        # `beds` has one entry per column, so a stray cell on a wider row indexed past
        # the end of it: IndexError, after the file had already been read.
        wide, _, _, _, _ = og.build_rows(
            ["OG1\tOs1\tZm1\tSTRAY\n"], ["rice", "maize"], "expand", [None] * 2, 4)
        check("a row wider than the header is cut to the columns, not indexed past them",
              wide, [["Os1", "Zm1"]])
        # a row SHORTER than the header is the ordinary case and still pads
        short, _, _, _, _ = og.build_rows(
            ["OG1\tOs1\tZm1\n"], ["rice", "maize", "sorghum"], "expand", [None] * 3, 4)
        check("a row shorter than the header pads to it",
              short, [["Os1", "Zm1", "."]])

    def test_compara(self):
        check = self.check
        # `compara`: Compara is the one ortholog source that publishes what
        # the inference measured, and those numbers are the whole reason for the
        # converter. Each has a unit the adapter or a ramp cares about, so a silent
        # change to one is a mis-colored track rather than a crash.
        # identity rides the [0,1] identity ramp; Compara's column is a percent
        check("identity is written as the fraction the ramp is defined on",
              cmp_mod.cell("91.5", 0.01), "0.915")
        # a missing dS is not a dS of zero: zero divides into a dN/dS of zero, which
        # reads as total purifying selection rather than as no measurement
        check("a rate Compara could not estimate is a missing cell, not a zero",
              [cmp_mod.cell(v) for v in ("NULL", "NA", "", ".", "\\N")],
              [".", ".", ".", ".", "."])
        check("a rate it did estimate survives the round trip",
              [cmp_mod.cell("0.4"), cmp_mod.cell("0")], ["0.4", "0"])
        # `copies` is the one column the converter counts rather than copies: how many
        # orthologs the reference gene has in that partner. Against a polyploid it IS the
        # ploidy, so every row of a fanned-out gene has to carry the same total — a
        # per-row count would read as "this link is one of one".
        og_src = os.path.join(tempfile.mkdtemp(), "h.tsv")
        with open(og_src, "w") as fh:
            fh.write("gene_stable_id\tprotein_stable_id\tspecies\tidentity\thomology_type\t"
                     "homology_gene_stable_id\thomology_protein_stable_id\thomology_species\t"
                     "homology_identity\tdn\tds\tgoc_score\twga_coverage\tis_high_confidence\thomology_id\n")
            for i, (g, o) in enumerate([("S1", "W1a"), ("S1", "W1b"), ("S1", "W1c"), ("S2", "W2")]):
                fh.write(f"{g}\tp\tsorghum_bicolor\t80\tortholog_one2many\t{o}\tp\t"
                         f"triticum_aestivum\t80\tNULL\tNULL\tNULL\tNULL\t1\t{i}\n")
        outdir = tempfile.mkdtemp()
        sys.argv = ["orthologs2jbrowse", og_src, "--reference", "sorghum_bicolor=sorghum",
                    "--species", "triticum_aestivum=wheat", "--outdir", outdir]
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            cmp_mod.main()
        # A gene with dozens of orthologs in one partner is a family, and `copies` has no
        # declared domain — one such gene drags the ramp so the 1-vs-3 distinction the
        # polyploid case is about lands in its bottom twentieth.
        family_src = os.path.join(tempfile.mkdtemp(), "fam.tsv")
        with open(family_src, "w") as fh:
            fh.write("gene_stable_id\tprotein_stable_id\tspecies\tidentity\thomology_type\t"
                     "homology_gene_stable_id\thomology_protein_stable_id\thomology_species\t"
                     "homology_identity\tdn\tds\tgoc_score\twga_coverage\tis_high_confidence\thomology_id\n")
            for i in range(3):
                fh.write(f"S1\tp\tsorghum_bicolor\t80\tortholog_one2many\tW{i}\tp\t"
                         f"triticum_aestivum\t80\tNULL\tNULL\tNULL\tNULL\t1\t{i}\n")
            for i in range(9):
                fh.write(f"S2\tp\tsorghum_bicolor\t80\tortholog_one2many\tF{i}\tp\t"
                         f"triticum_aestivum\t80\tNULL\tNULL\tNULL\tNULL\t1\t{100 + i}\n")
        famdir = tempfile.mkdtemp()
        sys.argv = ["orthologs2jbrowse", family_src, "--reference", "sorghum_bicolor=sorghum",
                    "--species", "triticum_aestivum=wheat", "--outdir", famdir]
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            cmp_mod.main()
        check("a gene past --max-copies is a family and contributes no rows",
              sorted({l.split("\t")[0] for l in
                      open(os.path.join(famdir, "sorghum.wheat.blocks")).read().splitlines()}),
              ["S1"])

        # `rows` is filled before --max-copies runs, so the converter's own "no rows
        # matched" guard has already passed by the time the filter can empty a table.
        # Left alone that is an exit 0 over a 0-byte .blocks, which tabix indexes and
        # JBrowse loads as a blank track -- the same silent wrong picture --max-copies
        # was added to prevent. Against a polyploid partner it is reachable rather than
        # theoretical: every reference gene there can sit above the threshold.
        allfam_src = os.path.join(tempfile.mkdtemp(), "allfam.tsv")
        with open(allfam_src, "w") as fh:
            fh.write("gene_stable_id\tprotein_stable_id\tspecies\tidentity\thomology_type\t"
                     "homology_gene_stable_id\thomology_protein_stable_id\thomology_species\t"
                     "homology_identity\tdn\tds\tgoc_score\twga_coverage\tis_high_confidence\thomology_id\n")
            for i in range(9):
                fh.write(f"S1\tp\tsorghum_bicolor\t80\tortholog_one2many\tF{i}\tp\t"
                         f"triticum_aestivum\t80\tNULL\tNULL\tNULL\tNULL\t1\t{i}\n")
        sys.argv = ["orthologs2jbrowse", allfam_src, "--reference", "sorghum_bicolor=sorghum",
                    "--species", "triticum_aestivum=wheat", "--outdir", tempfile.mkdtemp()]
        try:
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                cmp_mod.main()
            emptied_exit = "exit 0"
        except SystemExit as e:
            emptied_exit = "fatal" if str(e.code).startswith("--max-copies") else str(e.code)
        check("a run --max-copies empties entirely is fatal, not a 0-byte table",
              emptied_exit, "fatal")

        check("every row of a fanned-out gene carries the whole copy count",
              [l.split("\t")[-1] for l in
               open(os.path.join(outdir, "sorghum.wheat.blocks")).read().splitlines()],
              ["3", "3", "3", "1"])

        # the four required BED columns are a valid BED, and only the last field carries
        # the newline — the same trap `orthofinder` had
        bd2 = tempfile.mkdtemp()
        minimal2 = os.path.join(bd2, "minimal.bed")
        with open(minimal2, "w") as fh:
            fh.write("chr1\t100\t200\tENSG1\nchr1\t300\t400\tENSG2\t0\t+\n")
        check("compara read_bed_names reads a four-column BED",
              cmp_mod.read_bed_names(minimal2), {"ENSG1", "ENSG2"})

        # The one line on stdout is a contract: build_orthofinder_synteny.sh captures it
        # and hands it to three consumers, which position blockAssemblies/bedLocations
        # and index the table's columns. A comma-splitting consumer against this
        # space-separated line raised ValueError and, under `set -e`, killed the build
        # after the OrthoFinder run had already finished. Pin the separator here.
        d = tempfile.mkdtemp()
        src = os.path.join(d, "Orthogroups.tsv")
        with open(src, "w") as fh:
            fh.write("Orthogroup\tmaize.fa\trice.fa\n"
                     "OG0000000\tZm01a, Zm01b\tOs01g1\n")
        sys.argv = ["orthologs2jbrowse", src, "-o", os.path.join(d, "out.blocks")]
        stdout = io.StringIO()
        with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(io.StringIO()):
            og.main()
        check("the column order goes to stdout whitespace-separated, one line",
              stdout.getvalue().split("\n")[:2], ["maize rice", ""])
        check("the .blocks columns are in that order, not the caller's",
              open(os.path.join(d, "out.blocks")).read().splitlines(),
              ["Zm01a\tOs01g1", "Zm01b\tOs01g1"])

    def test_mcscanx_end_to_end(self):
        check = self.check
        # `mcscanx`: the one converter a reader is told to curl and run on
        # their own MCScanX output, and the only helper here with four output modes. Each
        # is driven end to end against a synthetic run, because the mistakes it can make
        # are all of the kind that writes a plausible file: the wrong genome's genes in a
        # BED, a block kept that joins the wrong pair, a coordinate off by one.
        d = tempfile.mkdtemp()
        gff = os.path.join(d, "xyz.gff")
        with open(gff, "w") as fh:
            # MCScanX's own 4-column gff: tagged chromosome, gene, start, end. g003 is
            # written end-before-start, which is how a minus-strand gene can arrive.
            fh.write("vv1\tVIT_g001\t1000\t2000\n"
                     "vv1\tVIT_g002\t3000\t4000\n"
                     "vv1\tVIT_g003\t6000\t5000\n"
                     "pp1\tPrupe_g001\t11000\t12000\n"
                     "pp1\tPrupe_g002\t13000\t14000\n"
                     "pp1\tPrupe_g003\t15000\t16000\n")
        coll = os.path.join(d, "xyz.collinearity")
        with open(coll, "w") as fh:
            fh.write("############### Parameters ###############\n"
                     "# MATCH_SCORE: 50\n"
                     "## Alignment 0: score=300.0 e_value=1e-30 N=3 vv1&pp1 plus\n"
                     "  0-  0:\tVIT_g001\tPrupe_g001\t  1e-30\n"
                     "  0-  1:\tVIT_g002\tPrupe_g002\t  1e-40\n"
                     "  0-  2:\tVIT_g003\tPrupe_g003\t  0\n"
                     # a self-synteny block, which the two-genome conversion drops and
                     # the one-genome conversion is entirely made of
                     "## Alignment 1: score=100.0 e_value=1e-10 N=2 vv1&vv1 minus\n"
                     "  1-  0:\tVIT_g001\tVIT_g003\t  1e-10\n")


        def run_mcx(*argv):
            sys.argv = ["orthologs2jbrowse", "--gff", gff, "--collinearity", coll,
                        *argv]
            with contextlib.redirect_stderr(io.StringIO()):
                mcx.main()


        def lines(*parts):
            return open(os.path.join(d, *parts)).read().splitlines()


        run_mcx("--species", "vv=grape", "--species", "pp=peach", "--outdir",
                os.path.join(d, "pair"))
        # BED is 0-based half-open where MCScanX's gff is 1-based inclusive, and the
        # reversed row is normalized rather than written as a negative-length feature
        check("the BED strips the chromosome tag and converts to 0-based",
              lines("pair", "grape.bed"),
              ["1\t999\t2000\tVIT_g001\t0\t+", "1\t2999\t4000\tVIT_g002\t0\t+",
               "1\t4999\t6000\tVIT_g003\t0\t+"])
        # the cross-species block only, with the self-synteny one dropped, and scores as
        # -log10(e_value) with e_value=0 capped
        check("only the cross-species block reaches the anchors file",
              lines("pair", "grape.peach.anchors"),
              ["###", "VIT_g001\tPrupe_g001\t30", "VIT_g002\tPrupe_g002\t40",
               "VIT_g003\tPrupe_g003\t1000"])
        # a simple row names the first and last gene of the block on each side, each by
        # its own coordinates, and takes its orientation from the block header
        check("the simple file reduces that block to one row",
              lines("pair", "grape.peach.anchors.simple"),
              ["VIT_g001\tVIT_g003\tPrupe_g001\tPrupe_g003\t3\t+"])

        # one --species keeps exactly what the two-genome case threw away
        run_mcx("--species", "vv=grape", "--outdir", os.path.join(d, "self"))
        check("a self-alignment run keeps the same-genome block and drops the rest",
              lines("self", "grape.grape.anchors"), ["###", "VIT_g001\tVIT_g003\t10"])

        # strand comes from the annotation the MCScanX input was built from, since
        # MCScanX's own gff has no strand column and every BED row would be `+`
        gff3 = os.path.join(d, "peach.gff3")
        with open(gff3, "w") as fh:
            fh.write("##gff-version 3\n"
                     "1\tx\tmRNA\t11000\t12000\t.\t-\t.\tID=mRNA:Prupe_g001\n")
        run_mcx("--species", "vv=grape", "--species", "pp=peach",
                "--strand-gff3", "peach=" + gff3, "--outdir", os.path.join(d, "strand"))
        check("--strand-gff3 recovers strand through the namespaced id too",
              lines("strand", "peach.bed")[0], "1\t10999\t12000\tPrupe_g001\t0\t-")

        # a refName the assembly does not have is the one mistake that draws an empty
        # track rather than erroring, so --fai turns it into a failure here instead
        fai = os.path.join(d, "peach.fa.fai")
        with open(fai, "w") as fh:
            fh.write("Pp01\t1000000\t10\t60\t61\n")
        try:
            run_mcx("--species", "vv=grape", "--species", "pp=peach",
                    "--fai", "peach=" + fai, "--outdir", os.path.join(d, "fai"))
            check("--fai rejects a BED refName the assembly lacks", "no exit", "SystemExit")
        except SystemExit as e:
            check("--fai rejects a BED refName the assembly lacks",
                  "not in" in str(e) and "Pp01" in str(e), True)

        # a tag that prefixes another silently files that genome's genes under the first
        try:
            run_mcx("--species", "v=grape", "--species", "vv=peach",
                    "--outdir", os.path.join(d, "ambig"))
            check("an ambiguous --species tag is rejected", "no exit", "SystemExit")
        except SystemExit as e:
            check("an ambiguous --species tag is rejected", "is a prefix of" in str(e), True)

        # --drop-overlapped-blocks: block 1 lies wholly inside the denser block 0 on vv1
        # and goes; block 2 is sparser still but only a third of it overlaps block 0
        gff = os.path.join(d, "overlap.gff")
        with open(gff, "w") as fh:
            fh.writelines(f"{chrom}\t{gene}{i}\t{1000 * i}\t{1000 * i + 500}\n"
                          for chrom, gene, n in (("vv1", "v", 10), ("pp1", "p", 4),
                                                 ("pp2", "q", 5), ("pp3", "r", 10))
                          for i in range(1, n + 1))
        coll = os.path.join(d, "overlap.collinearity")
        with open(coll, "w") as fh:
            for i, pairs in enumerate(([(1, "p1"), (2, "p2"), (3, "p3"), (4, "p4")],
                                       [(1, "q1"), (2, "q3"), (4, "q5")],
                                       [(1, "r1"), (4, "r4"), (7, "r7"), (10, "r10")])):
                fh.write(f"## Alignment {i}: score=1 e_value=1e-10 N={len(pairs)} vv1&pp1 plus\n")
                fh.writelines(f"  {i}-  {j}:\tv{v}\t{g}\t  1e-10\n"
                              for j, (v, g) in enumerate(pairs))
        run_mcx("--species", "vv=grape", "--species", "pp=peach",
                "--outdir", os.path.join(d, "overlap_off"))
        run_mcx("--species", "vv=grape", "--species", "pp=peach",
                "--drop-overlapped-blocks", "--outdir", os.path.join(d, "overlap_on"))
        check("without --drop-overlapped-blocks every block is written",
              [l.split("\t")[2] for l in lines("overlap_off", "grape.peach.anchors.simple")],
              ["p1", "q1", "r1"])
        check("--drop-overlapped-blocks drops the block inside a denser one only",
              [l.split("\t")[2] for l in lines("overlap_on", "grape.peach.anchors.simple")],
              ["p1", "r1"])

    def test_odp(self):
        check = self.check
        # `odp`: the linkage-group stack's opacity reads break_FET off the
        # hosted tables, so the column passes through as odp wrote it, a missing
        # group still '.'.
        with tempfile.TemporaryDirectory() as d:
            rbh = os.path.join(d, "A_B.rbh")
            with open(rbh, "w") as fh:
                fh.write("rbh\tA_gene\tB_gene\tgene_group\tA_scaf\tA_pos\tB_scaf\tB_pos"
                         "\tbreak_FET\tcolor\n"
                         "r1\ta1\tb1\tM\tA1\t10\tB1\t20\t0.001\t#C23D51\n"
                         "r2\ta2\tb2\tM\tA1\t30\tB2\t40\t0.4\t#C23D51\n")
            out = os.path.join(d, "A_B.blocks")
            subprocess.run([sys.executable, "-m", "orthologs2jbrowse", "odp", rbh, "-o", out,
                            "--bed-dir", d, "--species", "A", "B",
                            "--attributes", "gene_group", "color", "break_FET"],
                           check=True, capture_output=True)
            check("break_FET passes through beside the group and its color",
                  open(out).read().splitlines(),
                  ["a1\tb1\tM\t#C23D51\t0.001", "a2\tb2\tM\t#C23D51\t0.4"])

            raw = os.path.join(d, "raw.rbh")
            with open(raw, "w") as fh:
                fh.write("rbh\tA_gene\tB_gene\tA_scaf\tA_pos\tB_scaf\tB_pos\tbreak_FET\n"
                         "r1\ta1\tb1\tA1\t10\tB1\t20\t0.001\n"
                         "r2\ta9\tb9\tA1\t50\tB1\t60\t0.7\n")
            joined = os.path.join(d, "joined.blocks")
            subprocess.run([sys.executable, "-m", "orthologs2jbrowse", "odp", raw, "-o", joined,
                            "--bed-dir", d, "--species", "A", "B", "--alg", rbh,
                            "--attributes", "gene_group", "color", "break_FET"],
                           check=True, capture_output=True)
            check("--alg joins the group and its color, break_FET still off the table",
                  open(joined).read().splitlines(),
                  ["a1\tb1\tM\t#C23D51\t0.001", "a9\tb9\t.\t.\t0.7"])

            misspelt = subprocess.run(
                [sys.executable, "-m", "orthologs2jbrowse", "odp", rbh, "-o", out,
                 "--bed-dir", d, "--species", "A", "B",
                 "--attributes", "gene_group", "break_fet"],
                capture_output=True, text=True)
            check("an --attributes column the table lacks fails instead of writing '.'",
                  (misspelt.returncode != 0, "break_fet" in misspelt.stderr),
                  (True, True))

            # The odp build skips a pair whose .blocks exists, so a failed run must
            # leave none behind, or the re-run loads a blank genome.
            chrom = os.path.join(d, "B.chrom")
            with open(chrom, "w") as fh:
                fh.write("other\tB1\t+\t1\t100\n")
            dead_dir = os.path.join(d, "dead")
            dead_out = os.path.join(d, "dead.blocks")
            dead = subprocess.run(
                [sys.executable, "-m", "orthologs2jbrowse", "odp", rbh, "-o", dead_out,
                 "--bed-dir", dead_dir, "--species", "A", "B", "--chrom", f"B={chrom}"],
                capture_output=True, text=True)
            check("a genome whose ids resolve nowhere fails and writes nothing",
                  (dead.returncode != 0, "B: none of its ids resolve" in dead.stderr,
                   os.path.exists(dead_out), os.path.exists(dead_dir)),
                  (True, True, False, False))


if __name__ == '__main__':
    unittest.main()
