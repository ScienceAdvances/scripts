__author__ = "Victor Wong"
__copyright__ = "Copyright 2024, Victor Wong"
__email__ = "yehior@qq.com"
__license__ = "MIT"

import os
import sys
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

out_dir = ""
params = "--bdg "

ext = "_peaks.xls"
out_file = [o for o in snakemake.output if o.endswith(ext)][0]
out_name = os.path.basename(out_file[: -len(ext)])
out_dir = os.path.dirname(out_file)
control = snakemake.input.control
gsize=snakemake.params.gsize
pvalue=snakemake.params.pvalue
qvalue=snakemake.params.qvalue
keep_dup=snakemake.params.keep_dup
nolambda=snakemake.params.nolambda
broad=snakemake.params.broad
broad_cutoff=snakemake.params.broad_cutoff

if in_contr:
    opt_input = "-c {in_contr}"
if out_dir:
    out_dir = f"--outdir {out_dir}"
if control:
    control = f"--outdir {control} "
if gsize:
    params += f"--gsize {gsize} "
if pvalue:
    params += f"--pvalue {pvalue} "
if qvalue:
    params += f"--qvalue {qvalue} "
if nolambda:
    params += f"--nolambda "
if broad:
    params += f"--broad "
if broad_cutoff:
    params += f"--broad-cutoff {broad_cutoff} "

shell(
    "macs3 callpeak "
    "--treatment {snakemake.input.treatment} "
    "{control} "
    "--name {snakemake.params.name} "
    "--format {snakemake.params.format} "
    "--outdir {out_dir} "
    "{params} {log}"
)
